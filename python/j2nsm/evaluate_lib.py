"""評価の部品: 1 ステップ誤差、ロールアウト、基準線。"""
import numpy as np
import torch

from .data import group_of, sequences
from .physics import Physics


# ---------------------------------------------------------------- 1 ステップ誤差
def _rmse(a, mask=None):
    a = a if mask is None else a[mask]
    return float(np.sqrt(np.mean(a ** 2))) if a.size else float("nan")


@torch.no_grad()
def one_step_predictions(model, phys, d, split, bs=65536):
    """split の全遷移について、各手法の (Δq, Δdq) 予測を返す（numpy、物理量）。

    手法: model（学習済み）/ analytic（解析モデル RK4＋非弾性ストッパ）/ analytic_free（ストッパ無し）/
          zero（状態は変化しない）/ kinematic（Δq = dt·dq、Δdq = 0）
    """
    s = d["splits"][split]
    dt = d["dt"]
    xs = torch.from_numpy(s["X_spec"])
    X = torch.from_numpy(s["X"].astype(np.float64))
    model.eval()
    dq_model = []
    for i in range(0, len(xs), bs):
        dq_model.append(torch.stack(model.delta_physical(xs[i:i + bs], X[i:i + bs]), dim=1).double())
    dq_model = torch.cat(dq_model).numpy()
    q, dq, tau = X[:, 0], X[:, 1], X[:, 2]
    out = {"model": dq_model}
    for name, stop in (("analytic", True), ("analytic_free", False)):
        q1, v1 = phys.step(q, dq, tau, stop=stop)
        out[name] = torch.stack([q1 - q, v1 - dq], dim=1).numpy()
    out["zero"] = np.zeros_like(dq_model)
    out["kinematic"] = np.stack([dt * dq.numpy(), np.zeros(len(dq))], axis=1)
    return out


def one_step_table(model, phys, d, split):
    """群別・リミット拘束別の 1 ステップ RMSE（Δq [rad]、Δdq [rad/s]、加速度換算 [rad/s²]）。"""
    s = d["splits"][split]
    dt = d["dt"]
    y = s["Y"].astype(np.float64)
    preds = one_step_predictions(model, phys, d, split)
    grp = np.array([group_of(d["scenarios"]["pattern"][k]) for k in s["scenario_id"]])
    at = s["at_limit"]
    # 接触を、動いている接触（衝突・反発・離脱）と静止押し付けに分ける（MATLAB の検証と同じ定義）
    dyn = at & ((np.abs(s["X"][:, 1]) > 0.02) | (np.abs(y[:, 1]) > 0.005))
    sel = {"all": np.ones(len(y), bool), "free(接触なし)": ~at, "contact(接触あり)": at,
           "contact:dynamic(動的な接触)": dyn, "contact:static(静止押し付け)": at & ~dyn}
    for g in sorted(set(grp)):
        sel[f"group:{g}"] = grp == g
    table = {}
    for name, m in sel.items():
        if not m.any():
            continue
        row = {"n": int(m.sum()), "target_std_dq_rad_s2": _rmse(y[m, 1] - y[m, 1].mean()) / dt}
        for meth, p in preds.items():
            e = p - y
            row[meth] = {"rmse_dq_rad": _rmse(e[:, 0], m), "rmse_ddq_rad_s2": _rmse(e[:, 1], m) / dt,
                         "rmse_dq_rad_s": _rmse(e[:, 1], m)}
        table[name] = row
    return table, preds


# ---------------------------------------------------------------- ロールアウト
@torch.no_grad()
def rollout(step_fn, q0, dq0, tau):
    """q0, dq0: (B,)、tau: (B, H)。H ステップ自己回帰で進め、(B, H+1) の q, dq を返す。"""
    B, H = tau.shape
    q = torch.empty(B, H + 1, dtype=q0.dtype)
    dq = torch.empty(B, H + 1, dtype=q0.dtype)
    q[:, 0], dq[:, 0] = q0, dq0
    for k in range(H):
        q[:, k + 1], dq[:, k + 1] = step_fn(q[:, k], dq[:, k], tau[:, k])
    return q, dq


def make_step_fns(model, phys):
    def model_fn(q, dq, tau):
        return model.step_physical(q, dq, tau)

    def analytic_fn(q, dq, tau):
        return phys.step(q, dq, tau, stop=True)

    def hold_fn(q, dq, tau):  # 状態を変えない（persistence）
        return q, dq

    return {"model": model_fn, "analytic": analytic_fn, "persistence": hold_fn}


def make_windows(seqs, horizon, stride):
    """各シナリオから、開始点が stride 刻みで長さ horizon の窓を切り出す。"""
    W = {"q": [], "dq": [], "tau": [], "contact": [], "group": [], "name": []}
    for s in seqs:
        K = len(s["tau"])
        for a in range(0, K - horizon + 1, stride):
            W["q"].append(s["q"][a:a + horizon + 1]); W["dq"].append(s["dq"][a:a + horizon + 1])
            W["tau"].append(s["tau"][a:a + horizon]); W["contact"].append(bool(s["at_limit"][a:a + horizon].any()))
            W["group"].append(s["group"]); W["name"].append(s["name"])
    if not W["q"]:
        return None
    return {k: (np.array(v) if k != "name" else v) for k, v in W.items()}


def window_rollouts(model, phys, d, split, horizon=1000, stride=250):
    seqs = sequences(d, split)
    W = make_windows(seqs, horizon, stride)
    if W is None:
        return None
    fns = make_step_fns(model, phys)
    q0 = torch.from_numpy(W["q"][:, 0]); dq0 = torch.from_numpy(W["dq"][:, 0]); tau = torch.from_numpy(W["tau"])
    res = {}
    for name, fn in fns.items():
        q, dq = rollout(fn, q0, dq0, tau)
        res[name] = (q.numpy(), dq.numpy())
    return W, res


def horizon_curves(W, res, steps=(1, 10, 50, 100, 250, 500, 1000)):
    """窓全体・接触なし・接触ありについて、ステップ数ごとの RMSE（q は度、dq は rad/s）。"""
    out = {}
    sels = {"all": np.ones(len(W["contact"]), bool), "free": ~W["contact"], "contact": W["contact"]}
    for sname, m in sels.items():
        if not m.any():
            continue
        out[sname] = {"n_windows": int(m.sum())}
        for meth, (q, dq) in res.items():
            eq = q[m] - W["q"][m]
            ed = dq[m] - W["dq"][m]
            rq = np.sqrt(np.nanmean(eq ** 2, axis=0)); rd = np.sqrt(np.nanmean(ed ** 2, axis=0))
            div = float(np.mean(~np.isfinite(q[m][:, -1]) | (np.abs(q[m]).max(axis=1) > 10)))
            out[sname][meth] = {"rmse_q_deg": {str(h): float(np.degrees(rq[h])) for h in steps if h < len(rq)},
                                "rmse_dq_rad_s": {str(h): float(rd[h]) for h in steps if h < len(rd)},
                                "diverged_fraction": div,
                                "curve_q_deg": np.degrees(rq).tolist(), "curve_dq": rd.tolist()}
    return out


def horizon_by_group(W, res, h=1000):
    """群別の、ステップ h での q の RMSE（度）。"""
    out = {}
    for g in sorted(set(W["group"])):
        m = np.array(W["group"]) == g
        out[g] = {"n_windows": int(m.sum())}
        for meth, (q, dq) in res.items():
            e = q[m][:, h] - W["q"][m][:, h]
            out[g][meth] = float(np.degrees(np.sqrt(np.nanmean(e ** 2))))
    return out


# PD の設計（MATLAB の buildJ2ClosedLoop の既定: wn = 30 rad/s、ζ = 1、Kp = M·wn²、Kd = 2·ζ·M·wn、±tauPeak で飽和）
PD_WN, PD_ZETA = 30.0, 1.0


def pd_ff_reference(seq, p, dt):
    """参照角 qref（K 点）から、dqref, ddqref とフィードフォワードトルクを作る。

    τ_ff = M·ddq_ref − mgL·sin q_ref + Bv·dq_ref + Fc·tanh(dq_ref/ε)（MATLAB の j2Feedforward と同じ）
    qref は 1 kHz の記録なので、速度・加速度は数値微分（5 次多項式の滑らかな軌道なので十分な精度）。
    """
    qr = seq["qref"]
    dqr = np.gradient(qr, dt)
    ddqr = np.gradient(dqr, dt)
    tau_ff = p["M"] * ddqr - p["mgL"] * np.sin(qr) + p["Bv"] * dqr + p["Fc"] * np.tanh(dqr / p["eps"])
    return qr, dqr, tau_ff


@torch.no_grad()
def closed_loop_rollout(step_fn, seq, p, dt, tau_peak):
    """制御器（PD＋FF）を閉ループに入れ、プラントを step_fn に置き換えて全長を進める。

    開ループで記録トルクを再生する評価は、θ=0 まわりが不安定（成長率 √(mgL/M) ≈ 3.5/s）なため、数秒で
    正確なモデルでも発散して意味がない。サロゲートの本来の用途（制御器と結合）に合わせ、制御器を接続して評価する。
    """
    qr, dqr, tau_ff = pd_ff_reference(seq, p, dt)
    Kp, Kd = p["M"] * PD_WN ** 2, 2 * PD_ZETA * p["M"] * PD_WN
    K = len(qr)
    q = torch.tensor([seq["q"][0]], dtype=torch.float64)
    dq = torch.tensor([seq["dq"][0]], dtype=torch.float64)
    Q, V, U = np.empty(K + 1), np.empty(K + 1), np.empty(K)
    Q[0], V[0] = float(q), float(dq)
    for k in range(K):
        tau = tau_ff[k] + Kp * (qr[k] - float(q)) + Kd * (dqr[k] - float(dq))
        tau = min(max(tau, -tau_peak), tau_peak)
        U[k] = tau
        q, dq = step_fn(q, dq, torch.tensor([tau], dtype=torch.float64))
        Q[k + 1], V[k + 1] = float(q), float(dq)
        if not np.isfinite(Q[k + 1]) or abs(Q[k + 1]) > 10:        # 発散したら打ち切る（以降は NaN）
            Q[k + 1:], V[k + 1:], U[k + 1:] = np.nan, np.nan, np.nan
            break
    return Q, V, U, qr


def closed_loop_benchmark(model, phys, d, split="benchmark"):
    fns = {k: v for k, v in make_step_fns(model, phys).items() if k != "persistence"}
    p = d["physics"]
    out = []
    for s in sequences(d, split):
        if np.isnan(s["qref"]).any():
            continue
        r = {"name": s["name"], "pattern": s["pattern"], "n_steps": len(s["tau"]), "truth": (s["q"], s["dq"], s["tau"]),
             "pred": {}}
        for name, fn in fns.items():
            Q, V, U, qr = closed_loop_rollout(fn, s, p, d["dt"], p["tauPeak"])
            r["pred"][name] = (Q, V, U)
        r["qref"] = qr
        out.append(r)
    return out


def closed_loop_summary(rolls):
    rows = []
    for r in rolls:
        qt, dqt, ut = r["truth"]
        qr = r["qref"]
        row = {"name": r["name"], "pattern": r["pattern"], "n_steps": r["n_steps"],
               "truth_track_rmse_deg": float(np.degrees(np.sqrt(np.mean((qt[:-1] - qr) ** 2)))),
               "truth_track_max_deg": float(np.degrees(np.max(np.abs(qt[:-1] - qr))))}
        for meth, (Q, V, U) in r["pred"].items():
            e = Q - qt
            row[meth] = {"track_rmse_deg": float(np.degrees(np.sqrt(np.nanmean((Q[:-1] - qr) ** 2)))),
                         "track_max_deg": float(np.degrees(np.nanmax(np.abs(Q[:-1] - qr)))),
                         "vs_truth_rmse_q_deg": float(np.degrees(np.sqrt(np.nanmean(e ** 2)))),
                         "vs_truth_max_q_deg": float(np.degrees(np.nanmax(np.abs(e)))),
                         "vs_truth_rmse_dq_rad_s": float(np.sqrt(np.nanmean((V - dqt) ** 2))),
                         "tau_rmse_N_m": float(np.sqrt(np.nanmean((U - ut) ** 2))),
                         "diverged": bool(np.isnan(Q).any())}
        rows.append(row)
    return rows
