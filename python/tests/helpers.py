"""テスト用の共通部品: エクスポートと同じ構成の小さな HDF5 を作る。"""
import h5py
import numpy as np

PHYS = {"M": 6.0857, "mgL": 73.25, "Fc": 23.82, "Bv": 1.42, "eps": 0.0130,
        "qMin": -2.7576, "qMax": 1.1345, "qdMax": 5.236, "tauRated": 286.8, "tauPeak": 860.4}
SCALE = {"DTHETA_MAX": 5.236, "TAU_MAX": 860.4, "D_THETA_MAX": 5.236e-3, "D_DTHETA_MAX": 860.4 / 6.0857 * 1e-3}


def write_tiny_h5(path, exclude_at_limit=False, K=60, seed=0):
    """train 2 本・test 1 本・benchmark 1 本（各 K 遷移）。X=[q dq tau]、Y=[Δq Δdq]。"""
    rs = np.random.RandomState(seed)
    names, patterns, splits = [], [], []
    spec = [("train", "bln_normal"), ("train", "contact_fall"), ("test", "hold_-090"), ("benchmark", "ptp_v40")]
    seqs = []
    for k, (sp, pat) in enumerate(spec):
        q = -1.0 + np.cumsum(0.002 * rs.randn(K + 1))
        dq = 0.1 * rs.randn(K + 1)
        tau = 50 * rs.randn(K)
        qref = q[:K] + 0.001 if pat.startswith(("hold", "ptp")) else np.full(K, np.nan)
        at = np.zeros(K, bool)
        if pat.startswith("contact"):
            at[5:15] = True
        seqs.append((sp, q, dq, tau, qref, at))
        names.append(f"s{k}"); patterns.append(pat); splits.append(sp)
    with h5py.File(path, "w") as f:
        for sp in ("train", "val", "test", "benchmark"):
            rows = [(i, s) for i, s in enumerate(seqs) if s[0] == sp]
            X, Y, sid, st, at, qr = [], [], [], [], [], []
            for i, (_, q, dq, tau, qref, atl) in rows:
                X.append(np.c_[q[:-1], dq[:-1], tau]); Y.append(np.c_[np.diff(q), np.diff(dq)])
                sid += [i] * len(tau); st += list(range(len(tau))); at.append(atl); qr.append(qref)
            if not rows:  # val は空にできないので 1 行だけのダミー
                X, Y, sid, st, at, qr = [np.zeros((1, 3))], [np.zeros((1, 2))], [0], [0], [np.zeros(1, bool)], [np.zeros(1)]
            X, Y = np.vstack(X).astype(np.float32), np.vstack(Y).astype(np.float32)
            g = f.create_group(sp)
            g["X"], g["Y"] = X, Y
            g["X_spec"] = np.c_[np.sin(X[:, 0]), np.cos(X[:, 0]), X[:, 1] / SCALE["DTHETA_MAX"], X[:, 2] / SCALE["TAU_MAX"]].astype(np.float32)
            g["Y_spec"] = np.c_[Y[:, 0] / SCALE["D_THETA_MAX"], Y[:, 1] / SCALE["D_DTHETA_MAX"]].astype(np.float32)
            g["scenario_id"], g["step"] = np.array(sid, np.int32), np.array(st, np.int32)
            g["at_limit"], g["qref"] = np.concatenate(at).astype(np.uint8), np.concatenate(qr)
        for grp, dct in (("scale", SCALE), ("physics", PHYS)):
            for k, v in dct.items():
                f[f"{grp}/{k}"] = np.array([v])          # MATLAB は長さ 1 の配列として書く
        for k in ("x_mean", "x_std", "y_mean", "y_std"):
            f[f"stats/{k}"] = np.zeros(3 if k.startswith("x") else 2)
        sc = f.create_group("scenarios")
        dt = h5py.string_dtype()
        sc.create_dataset("name", data=names, dtype=dt); sc.create_dataset("pattern", data=patterns, dtype=dt)
        sc.create_dataset("type", data=["excitation"] * 4, dtype=dt); sc.create_dataset("split", data=splits, dtype=dt)
        sc["n_transition"] = np.full(4, K, np.int32)
        f.attrs["fs_data_hz"] = np.array([1000.0])
        f.attrs["exclude_at_limit"] = np.array([1.0 if exclude_at_limit else 0.0])
    return seqs
