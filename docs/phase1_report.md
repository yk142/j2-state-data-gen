# Phase 1 レポート — J2 単軸プラントと物理検証（issue #5）

実施日: 2026-10-06 / MATLAB R2025a

## 成果物
| ファイル | 内容 |
|---|---|
| `src/plant/buildJ2Plant.m` | J2 単軸プラント生成（`spike/phase0_fixA.m` の本実装）。I/F: `tau`[スカラ] → `q`, `dq` |
| `src/plant/setJ2InitialState.m` | 初期 (θ₀, θ̇₀) の設定（Position/Velocity Target） |
| `src/plant/simJ2.m` | トルク列（ゼロ次保持）を与えて実行し、`t, q, dq, tau` を返す |
| `test/test_j2_plant.m` | 物理検証 13 件（`matlab.unittest`） |

## 検証結果（13 件すべて合格、約 168 s）
| テスト | 確認内容 |
|---|---|
| interfaceIsScalar | revolute 1・weld 6、`ode14x`、刻み 125 µs、I/F がスカラ |
| gravityIsPendulum | J2 重力トルクが `−mgL·sin θ`（mgL=73.25 N·m）に残差 1e-4 N·m 以内で一致、θ=0 で 0 |
| holdsWithGravityCompensation ×5 | J2=−150, −90, −45, +30, +60° で RNEA 重力補償を与えると静止（ドリフト < 1e-3°、\|dq\| < 1e-3 rad/s）。Simscape の重力 = `c8RNEA` |
| uprightIsUnstable | θ=0 から ±5° の点は tau=0 で原点から遠ざかる。保持トルクの θ 微分が負 |
| frictionMasksNearUpright | θ0=0.5° では 3 s 後も 0.5° 以内（摩擦が重力を上回る） |
| freeFallReachesNegativeLimit | −90° から落下し −158° 付近に到達、めり込み 3° 以内、静止、速度上限内、NaN なし |
| freeFallReachesPositiveLimit | +60° から +65° 付近に到達して静止 |
| initialVelocityIsApplied | 初期角速度 1 rad/s が反映される（dq(1)=1.00000） |
| randomTorqueStaysFinite | ±150 N·m のランダム入力で NaN・発散なし |

## 実装上の判断
- **固定角の丸め**: J3=75.0684°（4 桁）のため、θ=0 の重力トルクは 3.8e-6 N·m 残る。許容誤差 1e-4 N·m とした（約 1.4e-6 rad の傾きに相当し、データ品質への影響はない）
- **初期角速度の警告**: 初期角速度が非ゼロのとき、Simscape が armature 回転子速度の初期条件を緩和して解き直し、`physmod:simscape:engine:core:dae_errors:IcDiagnostics` 警告を出す。関節速度ターゲット自体は満たされる（dq(1)=1.00000）ため、`simJ2` で非ゼロ時のみ抑制した。ターゲット優先度（High/Low）では回避できなかった。短い警告 1 行（初期条件の最初の求解が収束しませんでした）は残る
- `simJ2` のトルク入力は区間定数（ゼロ次保持）。8 kHz のソルバ刻みで保持される
- 親資産のコードは変更していない（`nJoints=1` の `p` と block 改名で `addC8*` を再利用）

## ゲート G1
FR-2 の全チェックが通った ✔。θ=0 が不安定平衡、重力トルクが −mgL·sin θ に一致する ✔。
