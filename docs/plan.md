# 実装計画 — C8-A901S J2 単軸 学習データ生成

要件は [requirements.md](requirements.md)。各フェーズを 1 つ以上の issue にし、「issue → ブランチ → 作業 → 報告 → 指示でコミット → 指示で PR」で回す。各フェーズ末尾に**検証ゲート**を置く。

## Phase 0: 基盤整備とスパイク

目的: 未決事項（requirements §7 の #1〜#3, #5）を実測で潰す。

- 0.1 リポジトリ構成と path 設定（親資産への参照を 1 箇所に集約）
  - `config/j2_config.m`（親資産パス、固定姿勢、シナリオ設定）
  - `src/`、`test/`、`data/`（`.gitignore` 対象）
- 0.2 親資産の動作確認: `c8_params`、`buildC8Plant` が実行できるか
- 0.3 固定軸の実現方式を 2 案で試作して比較
  - A: 生成 URDF の J1, J3〜J6 を `fixed` に後処理して `smimport`
  - B: 6 軸のまま他軸を高剛性で位置保持
  - 比較観点: 数値安定性、実行時間、J2 の応答が 6 軸（他軸ほぼ静止）と一致するか
- 0.4 実行時間の実測（8 kHz × 10 s）、リミット接触時の挙動観察

**ゲート G0**: 方式を決定し、J2 単軸が 8 kHz で発散なく動く。実行時間から総規模を見積もる。

## Phase 1: J2 単軸プラントと物理検証

- 1.1 `buildJ2Plant`: 方式確定版（Variant Subsystem は不要。真値のみ）
- 1.2 `setJ2InitialPose`: 初期 (θ₀, θ̇₀) の設定。θ̇₀ の与え方も確認する
- 1.3 物理検証テスト（FR-2）: 倒立点近傍の不安定性、重力トルク vs `c8RNEA`、自由落下、リミット到達

**ゲート G1**: FR-2 の全チェックが通る。q=0 の姿勢が資料の想定と一致する。

## Phase 2: 励振とシナリオ

- 2.1 `j2ExcitationModel`: 励振トルク入力＋ソフトバリア＋トルク飽和（`buildC8Excitation` の流儀。比率を設定化）
- 2.2 励振信号: バンドリミテッドノイズ、チャープ、ステップ、微小トルク（Signal Processing Toolbox は許諾確認済み）
- 2.3 初期状態サンプリング: Sobol（Statistics Toolbox は禁止のため `sobolset` 不可。**自前実装**または低食い違い列で代替）
- 2.4 `j2ScenarioSet`: フェーズ1・2 のシナリオ表（シード・種別・split）
- 2.5 倒立点近傍の安定化制御＋外乱シナリオ
- 2.6 PTP 連続 GO ベンチマーク（FR-4b）: J2 閉ループモデル（制御器＋FF）、ウェイポイント生成、速度スケール違いの複数本

**ゲート G2**: 全シナリオ種別で 1 本ずつ実行でき、バリア・飽和が意図どおり働く。

## Phase 3: データ生成とカバレッジ

- 3.1 `genJ2Dataset`: シナリオ単位キャッシュ・再開可能・8 kHz 生データ保存
- 3.2 1 kHz への間引き（`c8Downsample` 利用または同等実装）
- 3.3 split 分割、train のみの正規化統計、マニフェスト
- 3.4 `validateJ2Data`
- 3.5 `j2Coverage`（FR-8）:
  - `plotTimeSeries`（θ, θ̇, τ の 3 段）
  - `plotPhasePlanes`（θ–θ̇、θ–τ、θ̇–τ の 3 枚、密度表示）
  - `plotFreqCoverage`（PSD とスペクトログラム。`pwelch`/`spectrogram` 使用。不可の環境では Python フォールバック、requirements FR-8）
  - 励振系とベンチマークを分けて出力し、比較図も作る
- 3.6 レポート出力: PNG を `reports/issue-<N>/` に保存し、issue コメント用 Markdown（SHA 固定 raw URL の画像埋め込み付き）を生成する

**ゲート G3**: 検証関数エラー 0、カバレッジが目標値以上、PTP ベンチマークと励振系のカバレッジ比較が出せる。

## Phase 4: まとめ

- 4.1 全体の再現確認（シード固定で同一データ）
- 4.2 ドキュメント整備、後続課題（フェーズ3 能動学習、J3 拡張）の issue 化

## issue コメントへの画像埋め込み

GitHub の issue コメントに gh CLI から直接画像はアップロードできないため、次の手順にする。
1. PNG を `reports/issue-<N>/` に保存（ユーザー指示でコミット・push）
2. コミット SHA 固定の URL（`https://raw.githubusercontent.com/yk142/j2-state-data-gen/<SHA>/reports/issue-<N>/<file>.png`）を `![](...)` で埋め込む
3. issue コメントを投稿（これもユーザー指示後）

## 注意

- 8 kHz 真値の実行は長時間かかる。実行はバックグラウンドで行い、キャッシュで再開する
- 親資産は変更しない。必要な変更はこのリポジトリ側で完結させる
- Toolbox の使用可否は親資産 `docs/requirements.md` §3 に従い、`license('checkout',…)` で実機確認してから使う（Parallel Computing Toolbox は未許諾のため直列実行）
