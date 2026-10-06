# 他言語向けエクスポート レポート（issue #15）

実施日: 2026-10-06 / MATLAB R2025a、Python 3.11（h5py 3.16、numpy 2.4、scipy 1.17。検証は使い捨ての仮想環境で実施）

## 成果物
| ファイル | 内容 |
|---|---|
| `src/data/j2ExportFlat.m` | データセットを HDF5（`.h5`）または v7 mat（`.mat`）に書き出す |
| `test/test_j2_export.m` | 書き出しと読み戻しのテスト 7 件（合成データ） |
| `tools/verify_flat_export.py` | 書き出したファイルを Python から読んで整合を確認する |
| `docs/dataset.md` | エクスポートの構成と Python の読み込み例を追記 |

## 設計
- **形式は HDF5 と v7 mat の 2 つ**。HDF5 は h5py（Python）・MATLAB・他言語から読め、圧縮もできる。v7 mat は `scipy.io.loadmat` で読める（構造体配列を持つ MATLAB v7.3 は HDF5 の参照の入れ子になって読みにくいため、フラットな変数で書く）
- **配列の向き**: MATLAB は列優先なので、(N, d) の行列を転置して書き、**h5py から (N, d) に見える**ようにした。MATLAB の `h5read` は (d, N) で返す
- **インデックスは 0 始まり**（`scenario_id`、`step`）。他言語向けのため。MATLAB で使うときは +1 する
- 内容: 入力・目標（SI 単位）、引継ぎ資料の正規化形式（`X_spec`/`Y_spec`）、遷移ごとのメタ情報（シナリオ番号、遷移番号、リミット拘束フラグ）、シナリオ表、正規化統計、スケール、ルート属性（preset、固定姿勢、親資産のコミット等）
- **既定は単精度**（約 103 MB。倍精度は約 165 MB）。Δq は倍精度で計算してから丸めるので、単精度でも差分の精度は保たれる（テストで確認）
- 既定は `excludeAtLimit = false`（`j2BuildFlat` と同じ）。`at_limit` を持つので、読み込み側で除外もできる

## 検証
### MATLAB（`test_j2_export`、7 件すべて合格）
`j2BuildFlat` の結果との一致（倍精度で完全一致）、メタ情報からシナリオ表への逆引き（`scenario_id`/`step` から元の `q, dq, tau` を復元して一致）、シナリオ名・シード・遷移数・`phase`、統計・スケール・属性、単精度の許容誤差、`excludeAtLimit` と split の選択、v7 mat の読み戻し、未対応の拡張子のエラー。

### Python（`tools/verify_flat_export.py`、full 2,856,681 遷移で全項目合格）
| 出力 | サイズ | 遷移数 |
|---|---|---|
| `j2_flat_full.h5`（single） | 103.3 MB | train 1,958,000 / val 386,000 / test 446,000 / benchmark 66,681 |
| `j2_flat_full.mat`（v7、single） | 99.7 MB | 同上 |
| `j2_flat_full_double.h5`（double） | 165.0 MB | 同上 |

- 形状 (N,3) (N,2) (N,4) (N,2)、型 float32、NaN・Inf なし
- `X_spec`/`Y_spec` を `X`/`Y` と `/scale` から再計算して一致（sin²+cos² = 1 を含む）
- `scenario_id`/`step`/split の整合（各 split の全行が、その split のシナリオに属する）
- `/stats` が train の `X`, `Y` の平均・標準偏差と一致
- 遷移数の合計 2,856,681 = シナリオ別 `n_transition` の合計
- `.mat` を scipy で読み、`.h5` と `train_X` が完全に一致
- 3 ファイルすべてで `tools/verify_flat_export.py` が合格（倍精度は型 float64 で確認）
- 生成時間: h5 23 s、mat 17 s

## 注意点
- MATLAB が書く数値のスカラー（`/scale/*`、数値のルート属性）は、h5py からは**長さ 1 の配列**に見える。`[0]` で取り出す（`dataset.md` に明記）
- 文字列はバイト列（vlen）で返るので `.decode()` が必要
- 書き出したファイルは `data/export/`（git 管理外）。`j2ExportFlat` で再生成できる
- numpy の `.npz` は MATLAB から直接書けないため対象外（HDF5 と v7 mat で代替できる）

## 完了条件の確認（issue #15）
- Python から split ごとの学習行列（`X`/`Y`、正規化形式も）とメタ情報（シナリオ、パターン、リミット拘束）が読める ✔
- 書き出して読み戻し、MATLAB 側の `j2BuildFlat` と一致することを確認 ✔
- 読み込み例を `docs/dataset.md` に追加 ✔
