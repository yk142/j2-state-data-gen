# データセット仕様

`genJ2Dataset` が生成する J2 単軸の学習データセット（full: 394 シナリオ、2,856,681 遷移、1 kHz）の仕様。

## 物理系
EPSON C8-A901S の J2。J1, J4〜J6 = 0、**J3 = 75.0684°** で固定した単軸モデル。

- θ = 0 が倒立の不安定平衡点。重力トルクは `+mgL·sin θ`（θ を増やす向き）、mgL = 73.25 N·m
- 可動範囲 −158° 〜 +65°（−2.758 〜 +1.134 rad）、最大角速度 5.236 rad/s、トルク容量: 定格 286.8 N·m・瞬時最大 860.4 N·m
- 解析モデル: `M·ddq = τ + mgL·sin θ − Bv·dq − Fc·tanh(dq/ε)`（M = 6.09 kg·m²（armature 込み）、Fc = 23.8 N·m、Bv = 1.42 N·m/(rad/s)）。1 kHz データの加速度を残差比 3.5% で説明する。リミット接触は説明できない（残差比 1.08）
- シミュレーション: 8 kHz 固定ステップ（`ode14x`）、単位は SI（rad、rad/s、N·m）

## ファイル
| ファイル | 内容 |
|---|---|
| `data/j2_dataset_<preset>.mat` | 1 kHz のデータセット（MATLAB v7.3）。`ds.scen`, `ds.stats`, `ds.scale`, `ds.manifest` |
| `data/cache/<preset>_<name>.mat` | シナリオごとの 8 kHz 生データ（変数 `e`: `q dq tau qref`（single）、`q0 dq0 contact8k fs wall`） |

どちらも git 管理外。`tools/j2_run_all.sh` / `genJ2Dataset` で再生成できる（シード固定）。

## `ds.scen`（構造体配列、シナリオ 1 本 = 1 要素）
| フィールド | 内容 |
|---|---|
| `name`, `pattern`, `type` | シナリオ名、パターン名、種別（`excitation` / `hold` / `ptp`） |
| `phase` | 1, 2（収集計画のフェーズ）、3（リミット接触）、`'benchmark'` |
| `split` | `train` / `val` / `test` / `benchmark`（シナリオ単位。パターンごとに 70/15/15、PTP は benchmark） |
| `seed` | シナリオの乱数シード |
| `q`, `dq` | [K+1 × 1] 瞬時値（1 kHz） |
| `tau` | [K × 1] 区間 [t_k, t_{k+1}) の**平均トルク**（実印加トルク。バリア・飽和後） |
| `tauInst` | [K × 1] 区間先頭の瞬時トルク |
| `qref` | [K+1 × 1] 参照角（閉ループのみ。励振は空） |
| `atLimit` | [K × 1] logical。遷移の開始・終了状態のどちらかが、リミットから 0.1° 以内または外側 |
| `nTransition` | K |

**遷移 k**: 入力 `x_k = (q(k), dq(k), tau(k))` → 次状態 `(q(k+1), dq(k+1))`（Δt = 1 ms）。

## シナリオ（full）
| 種別 | パターン | 本数 | 長さ | 内容 |
|---|---|---|---|---|
| 励振（フェーズ 1） | `bln_normal` / `bln_high` / `bln_low` | 60 / 30 / 30 | 10 / 5 / 15 s | BLN 0.5 / 1.0 / 0.2 τ_ref、20 / 20 / 5 Hz。初期状態は Sobol 列 |
| 励振（フェーズ 1） | `freefall` | 6 | 3 s | 微小 BLN（0.05 τ_ref）での自由落下 |
| 励振（フェーズ 2） | `chirp` / `step` / `micro` | 10 / 20 / 20 | 20 / 2 / 10 s | 0.1→20 Hz 掃引 / 最大トルク（860 N·m、飽和）ステップ / 保持トルク＋微小 BLN |
| 姿勢保持 | `hold_+000` `+030` `+060` `-045` `-090` `-120` `-150` | 7 × 20 | 5 s | PD＋FF で姿勢まわりの滑らかな参照に追従し、BLN 外乱（0.3 τ_ref）を加える（倒立点近傍を取る） |
| リミット近傍 | `nearlimit` | 10 | 10 s | 参照をリミットの 4° 内側に置いた閉ループ保持、外乱 0.5 τ_ref |
| リミット接触（#11） | `contact_fall` / `contact_drive` / `contact_bln` | 12 / 24 / 20 | 3 / 4 / 10 s | バリアをリミットの外側に移し、衝突・反発・押し付け・離脱を取る |
| ベンチマーク | `ptp_v20` `v40` `v60` `v80` | 4 × 3 | 最大 60 s | PTP 連続 GO（6 点、停留 0.2 s、速度スケール 0.2〜0.8）。**評価用で学習には使わない** |

τ_ref = 定格トルク 286.8 N·m。励振モデルのソフトバリアは位置が可動域の 85%・速度が最大角速度の 60%（接触シナリオは位置 110%・速度 80%）。トルクは ±860.4 N·m で飽和する。

## 正規化・スケール
- `ds.stats`: `xMean xStd`（[q dq tau]）、`yMean yStd`（[Δq Δdq]）。**train の全遷移**（リミット拘束を含む）から算出
- `ds.scale`: 引継ぎ資料の定数を実機値で保存。`DTHETA_MAX` = 5.236 rad/s、`TAU_MAX` = 860.4 N·m、`D_THETA_MAX` = 5.236e-3 rad、`D_DTHETA_MAX` = 0.1414 rad/s
- `j2BuildFlat(ds, split, struct('features','spec'))`: `X = [sin q, cos q, dq/DTHETA_MAX, tau/TAU_MAX]`、`Y = [Δq/D_THETA_MAX, Δdq/D_DTHETA_MAX]`

## 使い方（MATLAB）
```matlab
ds = load('data/j2_dataset_full.mat');
[Xtr, Ytr] = j2BuildFlat(ds, 'train');                 % 学習
[Xva, Yva] = j2BuildFlat(ds, 'val');                   % 検証
[Xte, Yte] = j2BuildFlat(ds, 'test');                  % テスト
[Xb,  Yb ] = j2BuildFlat(ds, 'benchmark');             % PTP ベンチマーク（評価のみ）
[~, ~] = j2BuildFlat(ds, 'train', struct('excludeAtLimit', true));   % リミット拘束を除外して比較
```

## 他言語向けエクスポート（Python 等）
`ds.scen` は構造体配列で MATLAB 以外では読みにくいため、`j2ExportFlat` でフラットな形式に書き出す。

```matlab
ds = load('data/j2_dataset_full.mat');
j2ExportFlat(ds, 'data/export/j2_flat_full.h5');                       % HDF5（既定は single、gzip 圧縮、約 103 MB）
j2ExportFlat(ds, 'data/export/j2_flat_full.mat');                      % v7 mat（scipy.io.loadmat で読める、約 100 MB）
j2ExportFlat(ds, 'data/export/j2_flat_full.h5', struct('dtype','double'));   % 倍精度（約 165 MB）
```
オプション: `splits`（既定は 4 つすべて）、`excludeAtLimit`（既定 false）、`dtype`（`single`|`double`）、`deflate`（0〜9、既定 4）。

### HDF5 の構成（Python/h5py から見た形）
| パス | 形状 | 内容 |
|---|---|---|
| `/<split>/X` | (N, 3) | `[q rad, dq rad/s, tau N·m]`（SI 単位） |
| `/<split>/Y` | (N, 2) | `[Δq rad, Δdq rad/s]`（1 ms 後 − 現在） |
| `/<split>/X_spec` | (N, 4) | `[sin q, cos q, dq/DTHETA_MAX, tau/TAU_MAX]`（引継ぎ資料の正規化形式） |
| `/<split>/Y_spec` | (N, 2) | `[Δq/D_THETA_MAX, Δdq/D_DTHETA_MAX]` |
| `/<split>/scenario_id` | (N,) int32 | `/scenarios/*` の添字（**0 始まり**） |
| `/<split>/step` | (N,) int32 | シナリオ内の遷移番号（**0 始まり**） |
| `/<split>/at_limit` | (N,) uint8 | リミット拘束の遷移なら 1 |
| `/<split>/qref` | (N,) | 参照角 [rad]（姿勢保持・PTP の閉ループのみ。励振は NaN）。遷移 k の開始時刻の値。制御器を閉ループに入れた評価（学習試行 #23）に使う |
| `/scenarios/{name,pattern,type,phase,split}` | (S,) 文字列 | シナリオ表（S = 394）。`phase` は `'1'` `'2'` `'3'` `'benchmark'` |
| `/scenarios/{seed,n_transition}` | (S,) | シード、遷移数 |
| `/stats/{x_mean,x_std,y_mean,y_std}` | (3,) / (2,) | 正規化統計（train の全遷移） |
| `/scale/*` | (1,) | `DTHETA_MAX`, `TAU_MAX`, `D_THETA_MAX`, `D_DTHETA_MAX` |
| `/physics/*` | (1,) | 解析モデルの係数 `M`, `mgL`, `Fc`, `Bv`, `eps`, `qMin`, `qMax`, `qdMax`, `tauRated`, `tauPeak`（SI 単位。基準線の計算用） |

ルート属性: `preset`, `created`, `fs_data_hz`, `q_fixed_deg`, `parent_commit`, `matlab`, `exclude_at_limit`, `dtype`, `layout`, `X_columns`, `Y_columns`。
**MATLAB が書く数値のスカラーは、h5py からは長さ 1 の配列に見える**（例: `f['scale/TAU_MAX'][0]`、`f.attrs['fs_data_hz'][0]`）。

### Python での読み込み例
```python
import h5py, numpy as np

with h5py.File("data/export/j2_flat_full.h5", "r") as f:
    Xtr, Ytr = f["train/X"][:], f["train/Y"][:]            # (1958000, 3), (1958000, 2) float32
    Xva, Yva = f["val/X"][:],   f["val/Y"][:]
    Xs,  Ys  = f["train/X_spec"][:], f["train/Y_spec"][:]  # 正規化済み（4 入力 → 2 出力）
    scale    = {k: float(f["scale"][k][0]) for k in f["scale"]}
    names    = [x.decode() for x in f["scenarios/name"][:]]
    pattern  = np.array([x.decode() for x in f["scenarios/pattern"][:]])
    pat_of_row = pattern[f["train/scenario_id"][:]]         # 各行のパターン名（接触だけを取り出す等に使う）
    at_limit = f["train/at_limit"][:].astype(bool)
    contact_rows = at_limit                                  # リミット拘束の遷移

# 例: 接触を除いた学習データ
Xtr_free, Ytr_free = Xtr[~at_limit], Ytr[~at_limit]
```
v7 mat の場合:
```python
from scipy.io import loadmat
m = loadmat("data/export/j2_flat_full.mat")
Xtr, Ytr = m["train_X"], m["train_Y"]                       # (N, 3), (N, 2)
names = [str(x).strip() for x in m["scenario_name"].ravel()]
```

### 検証
`tools/verify_flat_export.py` が、書き出したファイルを Python（h5py / scipy）から読んで整合を確認する（形状・型・有限性、`X_spec`/`Y_spec` の再計算、split とシナリオの対応、`/stats` が train の平均・標準偏差と一致、遷移数の合計、`.h5` と `.mat` の一致）。
```bash
python tools/verify_flat_export.py data/export/j2_flat_full.h5 data/export/j2_flat_full.mat
```
full（2,856,681 遷移）で全項目合格（[export_report.md](export_report.md)）。

## 品質（検証 15 項目すべて合格）
NaN なし、|τ| ≤ 860.4 N·m、|dq| ≤ 4.54 rad/s、q は −158.7° 〜 +65.6°（接触シナリオはリミットを最大 0.65° 超える）、統計は train のみ由来、リミット拘束は全体の 6.6%（接触シナリオ内で 56.9%）。カバレッジ（20 ビン格子）: 3 次元占有率 25.6%、PTP が学習系の占有セルに入る割合 94.7%（[phase3_report.md](phase3_report.md)、[limit_contact_report.md](limit_contact_report.md)）。

## 注意点（学習に使う前に）
1. **リミット接触遷移の 93.6% は静止押し付け**（次状態が現在状態と同じ）。動的な接触は −側 7,047、+側 5,092 件。サンプル重みや間引きの検討が必要になりうる
2. **大トルクはバリア付近に偏る**: θ–τ 平面で |τ| > 300 N·m はバリア端（約 −141°、+48°）に集中し、内部の θ では薄い。τ の PSD の 9 Hz のピークはバリア（wn = 60 rad/s）の振動による
3. **倒立点近傍（|θ| < 約 19°）は摩擦が重力を上回り、トルク無しでは動かない**。この領域はトルク励振と姿勢保持シナリオで取っている
4. **固定姿勢は 1 つ（J3 = 75.0684°）**。J3 への依存は学習できない
5. 姿勢保持・PTP の PD は連続形（wn = 30 rad/s）で、ソルバ刻み 8 kHz で評価される（離散制御器ではない）
6. 初期角速度が非ゼロの実行で、Simscape の初期条件の警告が 1 行出る（結果には影響しない。`simJ2` / `runJ2Sim` で詳細警告のみ抑制）
7. `benchmark` は評価専用。学習に混ぜると PTP の被覆評価（学習系の占有セルとの比較）の意味がなくなる
