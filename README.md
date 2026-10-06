# j2-state-data-gen

EPSON C8-A901S の **J2 軸（肩ピッチ）** の Neural State Model（次時刻状態推論モデル）用の学習データを、Simscape Multibody の真値プラントから生成する MATLAB プロジェクト。

- 入力 `(θ, θ̇, τ)` → 出力 `(θ_{t+1}, θ̇_{t+1})`（1 kHz、Δt = 1 ms）
- J2 以外を固定した **J2 単軸**の物理モデル（J3 = 75.0684° で固定し、θ = 0 を倒立の不安定平衡点にしている）
- 励振（BLN・チャープ・ステップ・微小トルク）、安定化制御下の姿勢保持、リミット接触、PTP 連続 GO（ベンチマーク）の 394 シナリオ、2,856,681 遷移

## 前提
- MATLAB R2025a。使用するのは Simulink、Simscape、Simscape Multibody、Signal Processing Toolbox のみ（Deep Learning・System Identification・Statistics and Machine Learning Toolbox は使わない）
- 真値の出典として、親資産 `robotArmSurrogate-matlab`（`c8_params`、`buildC8Plant` 周辺、`c8RNEA`）を **path 参照のみ**で使う（コピー・改変しない）。既定の場所は `~/matlab-projects/robotArmSurrogate-matlab`、環境変数 `J2_PARENT_DIR` で変更できる
- Parallel Computing Toolbox は使えない前提（長時間の生成は別プロセスで分担する）
- 実行環境は MATLAB Home ライセンス（personal use only）

## 使い方
```matlab
addpath('config'); addpath(genpath('src'));
c = j2_setup_path();                         % 親資産と本リポジトリを path に追加

% 1. 動作確認: J2 単軸プラントを作って 1 シナリオ実行
mE  = buildJ2Excitation();                   % 開ループ励振モデル（バリア＋トルク飽和）
mC  = buildJ2ClosedLoop();                   % PD＋FF 閉ループモデル
S   = j2ScenarioSet('small');                % 動作確認用の小さなシナリオ表（25 本）
raw = runJ2Scenario(S(1), struct('excite',mE,'closed',mC));   % 8 kHz の生データ

% 2. データセット生成（シナリオ単位でキャッシュ、中断しても再開できる）
ds = genJ2Dataset(struct('preset','small')); % 約 8 分。'full' は約 4〜5 時間（394 本）
R  = validateJ2Data(ds);                     % 検証（15 項目）

% 3. 学習用の行列を取り出す
[X, Y] = j2BuildFlat(ds, 'train');                                   % X=[q dq tau], Y=[Δq Δdq]（SI 単位）
[Xs, Ys] = j2BuildFlat(ds, 'train', struct('features','spec'));      % 引継ぎ資料の正規化形式

% 4. カバレッジ評価（PNG 8 枚と指標）
M = j2Coverage(ds, 'reports/my_run');
```

長時間の `full` 生成は別プロセスで連続実行できる:
```bash
tools/j2_run_all.sh full reports/my_run      # 未生成シナリオの生成 → 組み立て・検証・評価
```
メモリの都合で 1 プロセス約 3.7 GB を使うため、並列にはしない。

## テスト
```matlab
runtests('test/test_j2_signals.m')       % 信号・Sobol・導出パラメータ（約 3 s）
runtests('test/test_j2_scenarioset.m')   % シナリオ表・参照軌道（約 15 s）
runtests('test/test_j2_dataset.m')       % 8k→1k 変換・フラット化・検証（約 6 s）
runtests('test/test_j2_coverage.m')      % カバレッジ評価（約 2 分）
runtests('test/test_j2_plant.m')         % 物理検証 13 件（約 3 分）
runtests('test/test_j2_scenarios.m')     % シナリオ実行 17 件（約 5.5 分）
```
再現確認（別プロセス）: `matlab -batch "addpath('tools'); j2_repro_check()"`

## 構成
| ディレクトリ | 内容 |
|---|---|
| `config/` | `j2_config.m`（親資産パス、固定姿勢、サンプリング、シード） |
| `src/model/` | `j2Params`（J2 の解析パラメータ: 実効慣性、mgL、摩擦、制限） |
| `src/plant/` | プラント・励振・閉ループモデルの生成、初期状態、実行 |
| `src/control/` | PTP 軌道、フィードフォワード、滑らかな参照 |
| `src/data/` | 励振信号、Sobol 列、シナリオ表と実行、1 kHz 変換、データセット生成、検証 |
| `src/analysis/` | カバレッジ評価（時系列・位相平面・周波数・接触） |
| `tools/` | 分担ワーカー、完了後の組み立て・評価、連続実行、再現確認 |
| `test/` | `matlab.unittest` のテスト |
| `spike/` | Phase 0 の検証用コード（記録として保存） |
| `docs/` | 要件・計画・各フェーズのレポート・データセット仕様 |
| `reports/` | 評価結果の PNG と指標（`issue-<N>/`） |
| `data/` | 生成物（git 管理外）: `cache/`（8 kHz 生データ）、`j2_dataset_<preset>.mat` |

## ドキュメント
- [docs/dataset.md](docs/dataset.md): データセットの仕様・使い方・注意点
- [docs/requirements.md](docs/requirements.md): 要件定義 / [docs/plan.md](docs/plan.md): 計画と進捗
- フェーズごとのレポート: [Phase 0](docs/phase0_report.md)、[Phase 1](docs/phase1_report.md)、[Phase 2](docs/phase2_report.md)、[Phase 3](docs/phase3_report.md)、[リミット接触](docs/limit_contact_report.md)、[Phase 4](docs/phase4_report.md)
- 作業は「issue → ブランチ → 作業 → 報告 → コミット・issue コメント → PR」の順で行う（`CLAUDE.md`）

## 主な設計判断
- **J2 単軸化**: 生成 URDF の J2 以外の joint を `fixed` に後処理して `smimport`。他軸 PD 保持と J2 応答が最大 0.0008° で一致し、約 5 倍速い（Phase 0）
- **固定姿勢 J3 = 75.0684°**: θ = 0 が不安定平衡になり、J2 の重力トルクが `−mgL·sin θ`（mgL = 73.25 N·m）に厳密一致する
- **1 kHz 変換**: 状態は瞬時値を間引き、トルクは 1 ms 区間の平均（遷移データのため状態にフィルタをかけない）
- **解析モデルとの整合**: `M·ddq = τ + mgL·sin θ − Bv·dq − Fc·tanh(dq/ε)`（M = 6.09 kg·m²）が 1 kHz データの加速度を残差 3.5% で説明する。リミット接触はこのモデルでは説明できない（残差比 1.08）
