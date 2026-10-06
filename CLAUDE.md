# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## 言語

ユーザーとのやり取りは日本語で行うこと。

## 参照先

- `README.md`（使い方）、`docs/dataset.md`（データセット仕様と注意点）、`docs/plan.md`（進捗）、`docs/phase*_report.md`（各フェーズの結果）
- 再現確認: `matlab -batch "addpath('tools'); j2_repro_check()"`（キャッシュなしで再生成し、保存済み 8 kHz データとのビット一致を確認）

## 作業サイクル（必須）

1. issue 作成 → ブランチ作成 → 作業
2. 作業終了後、**ユーザーに報告**（ここで止まる）
3. **ユーザー指示があって初めて**コミットし、結果を issue にコメント投稿
4. **ユーザー指示があって初めて** PR を作成 → ユーザーが承認

- **main への直接 push は禁止**。ユーザーの明示指示なしにコミット・push・PR 作成・issue 投稿をしない。
- GitHub リポジトリは現時点で未作成（ユーザーが用意する）。この間は git 操作・issue 作成をしない。
- 実装は「要件定義 → 計画 → 実装」の順（`~/matlab-projects/robotArmSurrogate-matlab` の流儀）。いきなり実装しない。

## 決定事項

- 実装は **MATLAB/Simulink**（このリポジトリ内）。`robotArmSurrogate-matlab` の `c8_params` やプラント関連は path 参照のみ（コピーしない）。Deep Learning / System Identification / Statistics and ML Toolbox は使用禁止。
- 対象は EPSON C8-A901S の J2。**他軸を固定した J2 単軸の別モデル**を作る。固定姿勢は J1, J4〜J6 = 0、**J3 = 75.0684°**（θ=0 を倒立平衡点にする値。J3=0 だと平衡点が θ≈+11° にずれる）。この姿勢で J2 の重力トルクは −mgL·sin θ（mgL=73.25 N·m）。
- 固定軸の実現は **URDF の該当 joint を fixed に後処理して smimport する方式（A）**に決定（`spike/phase0_fixA.m`、固定角は origin の rpy に畳み込む）。
- リミット接触: 基本の 338 本はバリアでリミットに接触しない。接触は `contact_fall/drive/bln`（56 本、バリアをリミットの外側に移し、速度バリアは 80%）で取る。接触遷移（`atLimit`）は既定で学習対象に含め、`j2BuildFlat` の `excludeAtLimit=true` で除外できる。
- 引継ぎ資料の `simulate_joint2` と Python + matlab.engine 構成は仮のもので、採用しない。

## コマンド（MATLAB R2025a）

```matlab
addpath('config'); addpath(genpath('src'));
c   = j2_setup_path();                      % 親資産と本リポジトリを path に追加し、設定を返す
mdl = buildJ2Plant();                       % J2 単軸プラントを生成（data/models/j2_plant.slx、約 15 s）
out = simJ2(mdl, tVec, tau, q0, dq0);       % トルク列を与えて実行 → out.t, q, dq, tau
results = runtests('test/test_j2_plant.m'); % 物理検証（13 件、約 3 分）
S   = j2ScenarioSet('full');                % シナリオ表（338 本）。'small' は動作確認用（19 本）
raw = runJ2Scenario(S(k), struct('excite',mE,'closed',mC));   % 1 シナリオを 8 kHz で実行
%   mE = buildJ2Excitation(); mC = buildJ2ClosedLoop();  （開ループ励振 / PD+FF 閉ループ）
```

- テスト: `test_j2_signals`（約 3 s）、`test_j2_scenarioset`（約 10 s）、`test_j2_scenarios`（約 250 s、実行を伴う）
- J2 の解析モデル: `M·ddq = τ + mgL·sin θ − Bv·dq − Fc·tanh(dq/ε)`（M=6.09 kg·m²（armature 込み）、mgL=73.25 N·m）。係数は `j2Params()`

- データセット生成: `ds = genJ2Dataset(struct('preset','small'))`（約 6 分）。`'full'` は約 6〜7 時間かかる。シナリオ単位で `data/cache/` にキャッシュされ、中断しても同じ呼び出しで再開できる
- 長時間生成は別プロセスで: `matlab -batch "addpath('tools'); j2_gen_worker('full', i, n)"`（n 分割の i 番目を担当）。ワーカー 1 本で約 3.7 GB 使うため、メモリが足りなければ並列にしない。生成後は `j2_finish('full')` で組み立て・検証・カバレッジ評価
- 検証: `validateJ2Data(ds)`、フラット化: `[X,Y] = j2BuildFlat(ds,'train')`（`features='spec'` で引継ぎ資料の正規化形式）、カバレッジ: `j2Coverage(ds, outDir)`（PNG 7 枚と指標）
- 1 kHz 変換は、状態は瞬時値の間引き、トルクは区間平均（フィルタなし。遷移データのため）。リミット拘束（リミットから 0.1° 以内または外側）の遷移にはフラグを付けるが、既定では除外しない
- 単一テストは `runtests('test/test_j2_plant.m','Name','test_j2_plant/gravityIsPendulum')`
- 親資産のパスは環境変数 `J2_PARENT_DIR` で上書き可（既定 `~/matlab-projects/robotArmSurrogate-matlab`）
- MATLAB MCP サーバー経由の実行は 120 s でバックグラウンド化される。長い実行は完了通知を待つ
- 8 kHz のシミュレーションは実時間の約 5〜6 倍かかる（Parallel Computing Toolbox 未許諾のため直列）

## 現状

Phase 4 まで完了（生成・検証・カバレッジ評価・再現確認・ドキュメント）。Phase 2 まで（プラント `src/plant/`、評価 `src/analysis/`、制御・参照軌道 `src/control/`、励振・シナリオ `src/data/`、解析パラメータ `src/model/`、テスト `test/`）。要件・計画は `docs/`。リポジトリ直下には日本語の引継ぎ資料（`Neural State Model 学習データ生成スクリプト 引継ぎ資料…md`、ステータス: 未着手）のみがあり、コード・ビルド・lint・テストの設定はまだ存在しない。git リポジトリではない。実装前に必ず資料を読むこと（各モジュールのコード例と収集計画が載っている）。以下は特に間違えやすい点をまとめたもの。

## 目的

6軸ロボットの **2軸目** の **Neural State Model**（次時刻状態推論モデル）用の学習データを、Python から MATLAB シミュレーター（`matlab.engine`。この環境では MATLAB MCP サーバーも利用可能）を呼び出して生成する。

- モデル入力: `(sin θ, cos θ, θ̇/DTHETA_MAX, τ/TAU_MAX)`、出力: 正規化した差分残差 `(θ_{t+1}-θ)/D_THETA_MAX, (θ̇_{t+1}-θ̇)/D_DTHETA_MAX`。
- 予定構成: `data_collection/` 配下に `config.py`、`excitation.py`、`initial_states.py`、`simulator_interface.py`、`collectors/phase{1_broadband,2_structured,3_active}.py`、`coverage.py`、`dataset.py`、`run_collection.py`。

## 物理系の規約（間違えやすい点）

- 垂直平面。**θ=0 は倒立位置（不安定平衡点）**、安定平衡点は θ=±π（垂下）。重力トルクは `τ_g = −mgL·sin θ`。
- 可動範囲は非対称: **−158° 〜 +65°**（−2.757 〜 1.134 rad）。安定点（−180°）は範囲外なので、自由落下させると必ずリミットに到達する。
- シミュレーションは **8 kHz**（Δt=125 µs）、学習データは **1 kHz**（Δt=1 ms、1/8 間引き）。8kHz の生データを保存してから間引くこと。
- MATLAB モデルがリミット処理を内包しているか、単位系（deg/rad、rpm/rad/s）、固定/可変ステップかは**未確認**。依存する前に確認し、内包していなければ Python 側でクリップ処理を追加する。リミット接触時のサンプルを除外するか含めるかも未決定。
- J、m、L、b、TAU_MAX、DTHETA_MAX（仮値 3.0 rad/s、TAU_MAX 300 N·m）は MATLAB モデルから取得すること。`D_DTHETA_MAX = (TAU_MAX/J)·0.001`。
- 資料中の `simulate_joint2(theta0, dtheta0, tau_seq, dt)` は仮のインターフェース。実際の関数名・引数を確認すること。

## 収集計画（約360軌跡、1kHz換算で約340万サンプル）

- フェーズ1（140軌跡）: 振幅・帯域を変えたバンドリミテッドノイズ＋自由落下。初期状態 (θ₀, θ̇₀) は Sobol 列で生成。
- フェーズ2（160軌跡）: 重要7姿勢（0°, +30°, +60°, −45°, −90°, −120°, −150°）× 20軌跡、加えてチャープ／最大トルクステップ／微小トルク／リミット付近の各パターン。θ=0 近傍は自由落下だと通過するだけなので、安定化制御＋外乱で収集する。
- フェーズ3（60軌跡）: 能動学習。フェーズ1+2でアンサンブルを学習し、予測分散の大きい領域を初期状態として3イテレーション追加収集。
- 保存形式: HDF5。軌跡ごとに `traj_XXXX` グループを作り、`theta`・`dtheta`・`tau`（float32）と `pattern`/`phase` 属性を持たせる。`build_flat_dataset` でフラットな (X, Y) 配列を生成。カバレッジ評価は (θ, θ̇, τ) の3次元ヒストグラム占有率と位相平面プロット。
- 初期状態はリミットから5°内側の安全マージンを取る。
