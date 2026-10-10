#!/usr/bin/env python3
"""学習済みモデルを評価する（1 ステップ誤差、ロールアウト、基準線との比較、図）。

    python python/evaluate.py --run baseline --report-dir reports/issue-23/baseline

基準線: 解析モデル（RK4＋非弾性ストッパ）、状態を変えない（persistence / zero / kinematic）。
ロールアウトは、(1) 1 秒窓の開ループ（真のトルク列を再生）、(2) PTP の閉ループ（PD＋FF にモデルを接続）の 2 種類。
全長の開ループ再生は、θ=0 まわりの不安定性（成長率 約 3.5/s）で正確なモデルでも発散するため行わない。
出力: <report-dir>/{metrics.json, training_curves.png, rollout_horizon.png, rollout_benchmark.png, error_maps.png}
"""
import argparse
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import torch  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from j2nsm import evaluate_lib as ev  # noqa: E402
from j2nsm.data import load  # noqa: E402
from j2nsm.model import build_model  # noqa: E402
from j2nsm.physics import Physics  # noqa: E402

plt.rcParams["font.family"] = ["Noto Sans CJK JP", "IPAexGothic", "TakaoPGothic", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False
COL = {"model": "tab:red", "analytic": "tab:blue", "persistence": "tab:gray", "truth": "k"}
LAB = {"model": "学習モデル", "analytic": "解析モデル（RK4＋ストッパ）", "persistence": "状態を変えない", "truth": "真値"}


def parse():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--data", default="data/export/j2_flat_full.h5")
    p.add_argument("--run", default="baseline")
    p.add_argument("--out-root", default="data/train")
    p.add_argument("--report-dir", default=None)
    p.add_argument("--horizon", type=int, default=1000)
    p.add_argument("--stride", type=int, default=250)
    p.add_argument("--threads", type=int, default=4)
    return p.parse_args()


def load_model(path):
    ck = torch.load(path, map_location="cpu")
    m = build_model(ck["config"])
    m.load_state_dict(ck["state_dict"])
    m.eval()
    return m, ck


def plot_training(run_dir, out):
    h = json.load(open(os.path.join(run_dir, "history.json")))
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.semilogy([r["epoch"] for r in h], [r["train_loss"] for r in h], label="train")
    ax.semilogy([r["epoch"] for r in h], [r["val_loss"] for r in h], label="val")
    ax.set_xlabel("エポック"); ax.set_ylabel("損失（標準化した出力の MSE）"); ax.grid(True, which="both", alpha=0.3); ax.legend()
    ax.set_title("学習曲線")
    fig.tight_layout(); fig.savefig(out, dpi=110); plt.close(fig)


def plot_horizon(curves, out, horizon):
    fig, axes = plt.subplots(2, 3, figsize=(16, 8))
    t = np.arange(horizon + 1)
    for j, sel in enumerate(("all", "free", "contact")):
        if sel not in curves:
            for i in range(2):
                axes[i, j].axis("off")
            continue
        for meth in ("persistence", "analytic", "model"):
            c = curves[sel][meth]
            axes[0, j].semilogy(t, np.maximum(c["curve_q_deg"], 1e-6), color=COL[meth], label=LAB[meth])
            axes[1, j].semilogy(t, np.maximum(c["curve_dq"], 1e-6), color=COL[meth], label=LAB[meth])
        ttl = {"all": "全窓", "free": "接触なしの窓", "contact": "接触を含む窓"}[sel]
        axes[0, j].set_title(f"{ttl}（{curves[sel]['n_windows']} 窓）")
        axes[0, j].set_ylabel("θ の RMSE [deg]"); axes[1, j].set_ylabel("θ̇ の RMSE [rad/s]")
        for i in range(2):
            axes[i, j].set_xlabel("ロールアウトのステップ数（1 ステップ = 1 ms）"); axes[i, j].grid(True, which="both", alpha=0.3)
    axes[0, 0].legend()
    fig.suptitle("ロールアウト誤差（真のトルク列を入力、初期状態だけ与えて自己回帰）")
    fig.tight_layout(); fig.savefig(out, dpi=100); plt.close(fig)


def plot_benchmark(rolls, out):
    """PTP ベンチマークの閉ループロールアウト（最も遅い軌道と最も速い軌道）。"""
    pick = [rolls[0], rolls[-1]]
    fig, axes = plt.subplots(4, 2, figsize=(15, 12), sharex="col")
    for j, r in enumerate(pick):
        qt, dqt, ut = r["truth"]
        t = np.arange(len(qt)) * 1e-3
        for meth in ("analytic", "model"):
            Q, V, U = r["pred"][meth]
            axes[0, j].plot(t, np.degrees(Q), color=COL[meth], lw=1, label=LAB[meth])
            axes[1, j].plot(t, V, color=COL[meth], lw=1)
            axes[2, j].plot(t[:-1], U, color=COL[meth], lw=1)
            axes[3, j].plot(t, np.degrees(Q - qt), color=COL[meth], lw=1)
        axes[0, j].plot(t, np.degrees(qt), "k--", lw=1.2, label=LAB["truth"])
        axes[1, j].plot(t, dqt, "k--", lw=1.2); axes[2, j].plot(t[:-1], ut, "k--", lw=1.2)
        axes[0, j].set_title(f"{r['name']}（{r['n_steps']} ステップ）")
        axes[0, j].set_ylabel("θ [deg]"); axes[1, j].set_ylabel("θ̇ [rad/s]"); axes[2, j].set_ylabel("τ [N·m]")
        axes[3, j].set_ylabel("θ − 真値 [deg]"); axes[3, j].set_xlabel("t [s]")
        for i in range(4):
            axes[i, j].grid(True, alpha=0.3)
    axes[0, 0].legend()
    fig.suptitle("PTP ベンチマークの閉ループロールアウト（PD＋FF の制御器にプラントとして接続。最も遅い軌道と最も速い軌道）")
    fig.tight_layout(); fig.savefig(out, dpi=100); plt.close(fig)


def plot_error_maps(d, preds, out, split="test"):
    s = d["splits"][split]
    dt, y = d["dt"], s["Y"].astype(np.float64)
    q, dq, tau = np.degrees(s["X"][:, 0]), s["X"][:, 1], s["X"][:, 2]
    fig, axes = plt.subplots(2, 2, figsize=(14, 9))
    qmin, qmax = np.degrees(d["physics"]["qMin"]), np.degrees(d["physics"]["qMax"])
    for i, meth in enumerate(("model", "analytic")):
        e = np.abs(preds[meth][:, 1] - y[:, 1]) / dt
        for j, (yy, ylab, yr) in enumerate(((dq, "θ̇ [rad/s]", (-5.3, 5.3)), (tau, "τ [N·m]", (-870, 870)))):
            H, xe, ye = np.histogram2d(q, yy, bins=(60, 60), range=((qmin, qmax), yr))
            S, _, _ = np.histogram2d(q, yy, bins=(xe, ye), weights=e)
            with np.errstate(invalid="ignore", divide="ignore"):
                M = np.where(H > 0, S / H, np.nan)
            im = axes[i, j].pcolormesh(xe, ye, M.T, shading="auto", norm=matplotlib.colors.LogNorm(vmin=1e-2, vmax=1e3), cmap="viridis")
            fig.colorbar(im, ax=axes[i, j], label="平均 |加速度誤差| [rad/s²]")
            axes[i, j].set_xlabel("θ [deg]"); axes[i, j].set_ylabel(ylab)
            axes[i, j].set_title(f"{LAB[meth]}: θ–{ylab.split()[0]}（{split}）")
    fig.suptitle("1 ステップ誤差（加速度換算）の分布")
    fig.tight_layout(); fig.savefig(out, dpi=100); plt.close(fig)


def main():
    a = parse()
    torch.set_num_threads(a.threads)
    rep = a.report_dir or os.path.join("reports", "eval", a.run)
    os.makedirs(rep, exist_ok=True)
    run_dir = os.path.join(a.out_root, a.run)
    model, ck = load_model(os.path.join(run_dir, "model.pt"))
    d = load(a.data)
    phys = Physics({k: d["physics"][k] for k in ("M", "mgL", "Fc", "Bv", "eps", "qMin", "qMax")}, d["dt"])
    metrics = {"run": a.run, "kind": model.kind, "epoch_best": ck["epoch"], "structured": model.structured,
               "n_param": sum(p.numel() for p in model.parameters())}

    print("1 ステップ誤差 ...", flush=True)
    metrics["one_step"] = {}
    preds_test = None
    for sp in ("val", "test", "benchmark"):
        tab, preds = ev.one_step_table(model, phys, d, sp)
        metrics["one_step"][sp] = tab
        if sp == "test":
            preds_test = preds
    print("ウィンドウのロールアウト ...", flush=True)
    metrics["rollout_windows"] = {}
    for sp in ("test",):
        r = ev.window_rollouts(model, phys, d, sp, a.horizon, a.stride)
        if r is None:
            continue
        W, res = r
        cur = ev.horizon_curves(W, res)
        metrics["rollout_windows"][sp] = {"horizon": a.horizon, "stride": a.stride, "n_windows": int(len(W["q"])),
                                          "by_group_rmse_q_deg_at_horizon": ev.horizon_by_group(W, res, a.horizon),
                                          "curves": {k: {m: {kk: vv for kk, vv in v.items() if not kk.startswith("curve_")}
                                                         if isinstance(v, dict) else v for m, v in c.items()}
                                                     for k, c in cur.items()}}
        plot_horizon(cur, os.path.join(rep, "rollout_horizon.png"), a.horizon)
    print("PTP ベンチマークの閉ループロールアウト ...", flush=True)
    rolls = ev.closed_loop_benchmark(model, phys, d, "benchmark")
    metrics["closed_loop_benchmark"] = ev.closed_loop_summary(rolls)
    plot_benchmark(rolls, os.path.join(rep, "rollout_benchmark.png"))
    plot_training(run_dir, os.path.join(rep, "training_curves.png"))
    plot_error_maps(d, preds_test, os.path.join(rep, "error_maps.png"))
    with open(os.path.join(rep, "metrics.json"), "w") as f:
        json.dump(metrics, f, indent=1, ensure_ascii=False)
    print_summary(metrics)
    print("保存:", rep)


def print_summary(m):
    print(f"\n== {m['run']}（種類={m.get('kind', 'nsm')}、構造付き={m['structured']}、パラメータ {m['n_param']:,}、best epoch {m['epoch_best']}）")
    for sp in ("test", "benchmark"):
        t = m["one_step"][sp]
        print(f"[1 ステップ {sp}] 加速度換算 RMSE [rad/s²]（目標の標準偏差 {t['all']['target_std_dq_rad_s2']:.1f}）")
        for k in ("all", "free(接触なし)", "contact(接触あり)", "contact:dynamic(動的な接触)", "contact:static(静止押し付け)"):
            if k in t:
                r = t[k]
                print(f"   {k:26s} n={r['n']:>8,}  model {r['model']['rmse_ddq_rad_s2']:8.3f}  analytic {r['analytic']['rmse_ddq_rad_s2']:8.3f}  zero {r['zero']['rmse_ddq_rad_s2']:8.3f}")
    w = m["rollout_windows"].get("test")
    if w:
        print(f"[ロールアウト test {w['n_windows']} 窓、{w['horizon']} ステップ] θ RMSE [deg]")
        for sel, c in w["curves"].items():
            print(f"   {sel:8s} n={c['n_windows']:4d}  " + "  ".join(
                f"{meth}: {c[meth]['rmse_q_deg'][str(w['horizon'])]:.3g}" for meth in ("model", "analytic", "persistence")
                if str(w["horizon"]) in c[meth]["rmse_q_deg"]) + f"  発散率 model {c['model']['diverged_fraction']:.3f}")
    print("[PTP ベンチマーク 閉ループロールアウト] θ の真値との RMSE [deg]（model / analytic）、追従誤差 RMSE [deg]（真値 / model / analytic）")
    for r in m["closed_loop_benchmark"]:
        print(f"   {r['name']:14s} {r['n_steps']:6d} 歩  vs真値 {r['model']['vs_truth_rmse_q_deg']:8.4f} / {r['analytic']['vs_truth_rmse_q_deg']:8.4f}"
              f"   追従 {r['truth_track_rmse_deg']:7.4f} / {r['model']['track_rmse_deg']:7.4f} / {r['analytic']['track_rmse_deg']:7.4f}"
              f"{'  [model 発散]' if r['model']['diverged'] else ''}")


if __name__ == "__main__":
    main()
