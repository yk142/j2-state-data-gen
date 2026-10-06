"""j2nsm: J2 単軸の Neural State Model（1 ステップ遷移モデル）の学習・評価。

エクスポート済みの HDF5（MATLAB の j2ExportFlat）を読み、PyTorch で学習し、
1 ステップ誤差とロールアウトで評価する。解析モデル（physics）を基準線とする。
"""
