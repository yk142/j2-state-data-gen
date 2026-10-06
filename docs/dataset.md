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

## 使い方
```matlab
ds = load('data/j2_dataset_full.mat');
[Xtr, Ytr] = j2BuildFlat(ds, 'train');                 % 学習
[Xva, Yva] = j2BuildFlat(ds, 'val');                   % 検証
[Xte, Yte] = j2BuildFlat(ds, 'test');                  % テスト
[Xb,  Yb ] = j2BuildFlat(ds, 'benchmark');             % PTP ベンチマーク（評価のみ）
[~, ~] = j2BuildFlat(ds, 'train', struct('excludeAtLimit', true));   % リミット拘束を除外して比較
```
MATLAB 以外で使う場合は、行列を v7 形式で保存すると `scipy.io.loadmat` で読める: `save('flat.mat','Xtr','Ytr','Xva','Yva','Xte','Yte','-v7')`（`.mat` v7.3 は HDF5 だが、構造体配列は参照の入れ子になり読みにくい）。

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
