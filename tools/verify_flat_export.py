#!/usr/bin/env python3
"""j2ExportFlat で書き出したファイルを Python（h5py / scipy）から読み、整合を検証する。

使い方:
    python tools/verify_flat_export.py data/export/j2_flat_full.h5 [data/export/j2_flat_full.mat]

検証内容:
  - 各 split の形状 (N,3) (N,2) (N,4) (N,2)、型、有限性
  - 遷移数の合計が /scenarios/n_transition の合計と一致（at_limit を除外していない場合）
  - 正規化統計 /stats が train の X, Y の平均・標準偏差と一致（train の全遷移から算出しているため）
  - X_spec / Y_spec が X / Y と /scale から再計算した値と一致（sin²+cos²=1 を含む）
  - scenario_id / step / split の整合（train の行は split == 'train' のシナリオに属する）
  - .mat（v7）も与えた場合、scipy.io.loadmat で読めて h5 と同じ値になる
依存: numpy, h5py, scipy（.mat を検証する場合）
"""
import sys
import numpy as np
import h5py


def check(cond, msg):
    print(('  [OK] ' if cond else '  [NG] ') + msg)
    return bool(cond)


def main(h5path, matpath=None):
    ok = True
    with h5py.File(h5path, 'r') as f:
        print(f'ファイル: {h5path}')
        print('ルート属性:', {k: (v.decode() if isinstance(v, bytes) else (float(np.ravel(v)[0]) if np.ndim(v) else v)) for k, v in f.attrs.items()
                              if k in ('preset', 'created', 'fs_data_hz', 'dtype', 'parent_commit', 'exclude_at_limit')})
        splits = [s for s in ('train', 'val', 'test', 'benchmark') if s in f]
        scale = {k: float(f['scale'][k][()]) if f['scale'][k].shape == () else float(f['scale'][k][0]) for k in f['scale']}
        names = [x.decode() if isinstance(x, bytes) else x for x in f['scenarios/name'][:]]
        sc_split = np.array([x.decode() if isinstance(x, bytes) else x for x in f['scenarios/split'][:]])
        n_tr = f['scenarios/n_transition'][:]
        total = 0
        for sp in splits:
            g = f[sp]
            X, Y, Xs, Ys = g['X'], g['Y'], g['X_spec'], g['Y_spec']
            N = X.shape[0]
            total += N
            print(f'[{sp}] N = {N:,}')
            ok &= check(X.shape == (N, 3) and Y.shape == (N, 2) and Xs.shape == (N, 4) and Ys.shape == (N, 2),
                        f'形状 X{X.shape} Y{Y.shape} X_spec{Xs.shape} Y_spec{Ys.shape}')
            ok &= check(X.dtype == np.float32 or X.dtype == np.float64, f'型 {X.dtype}')
            Xa, Ya, Xsa, Ysa = X[:], Y[:], Xs[:], Ys[:]
            ok &= check(all(np.isfinite(a).all() for a in (Xa, Ya, Xsa, Ysa)), '有限（NaN/Inf なし）')
            # spec 形式の再計算
            ok &= check(np.allclose(Xsa[:, 0]**2 + Xsa[:, 1]**2, 1.0, atol=1e-5), 'sin²+cos² = 1')
            ok &= check(np.allclose(Xsa[:, 2], Xa[:, 1] / scale['DTHETA_MAX'], atol=1e-5), 'X_spec[:,2] = dq / DTHETA_MAX')
            ok &= check(np.allclose(Xsa[:, 3], Xa[:, 2] / scale['TAU_MAX'], atol=1e-5), 'X_spec[:,3] = tau / TAU_MAX')
            ok &= check(np.allclose(Ysa[:, 0], Ya[:, 0] / scale['D_THETA_MAX'], rtol=1e-4, atol=1e-4), 'Y_spec[:,0] = Δq / D_THETA_MAX')
            ok &= check(np.allclose(Ysa[:, 1], Ya[:, 1] / scale['D_DTHETA_MAX'], rtol=1e-4, atol=1e-4), 'Y_spec[:,1] = Δdq / D_DTHETA_MAX')
            # メタ情報
            sid, stp, at = g['scenario_id'][:], g['step'][:], g['at_limit'][:]
            ok &= check(sid.min() >= 0 and sid.max() < len(names), 'scenario_id は 0 始まりで範囲内')
            ok &= check(bool(np.all(sc_split[sid] == sp)), f"全行が split == '{sp}' のシナリオに属する")
            ok &= check(bool(np.all(stp >= 0) and np.all(stp < n_tr[sid])), 'step は 0 ≤ step < n_transition')
            print(f'     at_limit の割合 {at.mean()*100:.2f}%，|dq| 最大 {np.abs(Xa[:,1]).max():.3f} rad/s，|tau| 最大 {np.abs(Xa[:,2]).max():.1f} N·m')
            if sp == 'train':
                st = {k: f['stats'][k][:] for k in f['stats']}
                ok &= check(np.allclose(st['x_mean'], Xa.astype(np.float64).mean(0), rtol=1e-4, atol=1e-5), '/stats/x_mean = train の X の平均')
                ok &= check(np.allclose(st['x_std'], Xa.astype(np.float64).std(0, ddof=1), rtol=1e-4, atol=1e-5), '/stats/x_std = train の X の標準偏差')
                ok &= check(np.allclose(st['y_mean'], Ya.astype(np.float64).mean(0), rtol=1e-3, atol=1e-7), '/stats/y_mean = train の Y の平均')
                ok &= check(np.allclose(st['y_std'], Ya.astype(np.float64).std(0, ddof=1), rtol=1e-3, atol=1e-7), '/stats/y_std = train の Y の標準偏差')
        # MATLAB が書く数値属性は h5py から長さ 1 の配列に見える（スカラーではない）
        if not int(np.ravel(f.attrs.get('exclude_at_limit', 0))[0]):
            ok &= check(total == int(n_tr.sum()), f'遷移数の合計 {total:,} = シナリオ別 n_transition の合計 {int(n_tr.sum()):,}')
        if matpath:
            from scipy.io import loadmat
            m = loadmat(matpath)
            Xm = m['train_X']
            ok &= check(Xm.shape == f['train/X'].shape, f'.mat の train_X 形状 {Xm.shape}')
            ok &= check(np.array_equal(Xm, f['train/X'][:]), '.mat と .h5 の train_X が完全に一致')
            ok &= check(list(m['scenario_name'].ravel()) == names or [str(x).strip() for x in m['scenario_name'].ravel()] == names,
                        '.mat のシナリオ名が一致')
            ok &= check(int(m['scenario_n_transition'].sum()) == int(n_tr.sum()), '.mat の n_transition の合計が一致')
    print('結果:', '合格' if ok else '不合格')
    return 0 if ok else 1


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(2)
    sys.exit(main(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else None))
