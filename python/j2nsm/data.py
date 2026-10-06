"""エクスポート済み HDF5 の読み込みと、シナリオ系列の復元。"""
import h5py
import numpy as np

SPLITS = ("train", "val", "test", "benchmark")


def _s(x):
    return x.decode() if isinstance(x, bytes) else str(x)


def load(path, splits=SPLITS):
    """HDF5 を読み、辞書で返す。

    d['splits'][split] : X (N,3), Y (N,2), X_spec (N,4), Y_spec (N,2), scenario_id, step, at_limit(bool), qref (N,)
    d['scale'], d['physics'] : スカラー定数の辞書（SI 単位）
    d['stats'] : x_mean, x_std, y_mean, y_std（train の全遷移）
    d['scenarios'] : name, pattern, type, split, n_transition
    d['dt'] : 1 ステップの時間 [s]
    """
    d = {"splits": {}}
    with h5py.File(path, "r") as f:
        for sp in splits:
            g = f[sp]
            e = {k: g[k][:] for k in ("X", "Y", "X_spec", "Y_spec", "scenario_id", "step")}
            e["at_limit"] = g["at_limit"][:].astype(bool)
            e["qref"] = g["qref"][:].astype(np.float64) if "qref" in g else np.full(len(e["step"]), np.nan)
            d["splits"][sp] = e
        # MATLAB が書く数値スカラーは長さ 1 の配列に見える
        d["scale"] = {k: float(np.ravel(f["scale"][k][()])[0]) for k in f["scale"]}
        d["physics"] = {k: float(np.ravel(f["physics"][k][()])[0]) for k in f["physics"]}
        d["stats"] = {k: f["stats"][k][:] for k in f["stats"]}
        sc = f["scenarios"]
        d["scenarios"] = {
            "name": [_s(x) for x in sc["name"][:]],
            "pattern": [_s(x) for x in sc["pattern"][:]],
            "type": [_s(x) for x in sc["type"][:]],
            "split": [_s(x) for x in sc["split"][:]],
            "n_transition": sc["n_transition"][:].astype(int),
        }
        d["dt"] = 1.0 / float(np.ravel(f.attrs["fs_data_hz"])[0])
        d["exclude_at_limit"] = bool(np.ravel(f.attrs.get("exclude_at_limit", 0))[0])
    return d


def group_of(pattern):
    """パターン名を評価用の群にまとめる。"""
    for prefix, g in (("bln_", "bln"), ("hold_", "hold"), ("contact_", "contact"), ("ptp_", "ptp")):
        if pattern.startswith(prefix):
            return g
    return pattern  # freefall, chirp, step, micro, nearlimit


def sequences(d, split):
    """split 内のシナリオごとの時系列（q, dq: K+1 点、tau, at_limit: K 点）を復元する。

    エクスポートの行はシナリオごとに連続し、遷移番号の昇順に並ぶ（リミット拘束を除外していない場合）。
    """
    if d["exclude_at_limit"]:
        raise ValueError("リミット拘束を除外したエクスポートからは時系列を復元できません")
    s = d["splits"][split]
    sid, step, X, Y = s["scenario_id"], s["step"], s["X"].astype(np.float64), s["Y"].astype(np.float64)
    out = []
    bounds = np.flatnonzero(np.diff(sid)) + 1
    for lo, hi in zip(np.r_[0, bounds], np.r_[bounds, len(sid)]):
        k = int(sid[lo])
        assert np.array_equal(step[lo:hi], np.arange(hi - lo)), "遷移番号が連続していません"
        assert hi - lo == d["scenarios"]["n_transition"][k], "遷移数がシナリオ表と一致しません"
        q = np.r_[X[lo:hi, 0], X[hi - 1, 0] + Y[hi - 1, 0]]
        dq = np.r_[X[lo:hi, 1], X[hi - 1, 1] + Y[hi - 1, 1]]
        out.append({
            "sid": k, "name": d["scenarios"]["name"][k], "pattern": d["scenarios"]["pattern"][k],
            "group": group_of(d["scenarios"]["pattern"][k]),
            "q": q, "dq": dq, "tau": X[lo:hi, 2], "at_limit": s["at_limit"][lo:hi],
            "qref": s["qref"][lo:hi],   # 閉ループ（姿勢保持・PTP）のみ。励振は NaN
        })
    return out
