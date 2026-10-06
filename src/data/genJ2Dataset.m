function ds = genJ2Dataset(opts)
%GENJ2DATASET シナリオを実行して J2 単軸の学習データセットを生成する。
%   ds = GENJ2DATASET(opts)
%     opts.preset    'full'（既定）| 'small'
%     opts.outDir    保存先（既定 <repo>/data）。cache/ と データセット .mat を置く
%     opts.shard     [i n]: n 分割したうちの i 番目だけ新規生成する（複数の MATLAB プロセスで分担）
%     opts.maxScen   1 回の呼び出しで新規生成するシナリオ数の上限（既定 Inf）
%     opts.cacheOnly true なら生成だけ行い、データセットの組み立てはしない（分担実行用）
%     opts.save      データセットを .mat に保存するか（既定 true）
%   戻り値 ds: 全シナリオのキャッシュが揃っていれば組み立てたデータセット、無ければ []。
%
%   シナリオ単位でキャッシュする（8 kHz 生データを single で保存）ので、中断しても再開できる。
%   Parallel Computing Toolbox は使えないため、別プロセスで shard を指定して分担する
%   （tools/j2_gen_worker.m）。各プロセスは自分の作業フォルダでモデルを生成する。
%
%   組み立て:
%     1. 8 kHz → 1 kHz（j2Downsample）。生データは cache に残す
%     2. リミット拘束フラグ（遷移の開始・終了状態のどちらかが、リミットから 0.1° 以内または外側）を付ける。
%        リミット接触シナリオ（#11）を追加したため、拘束サンプルは既定で学習対象に含める
%        （除外して比較したい場合は j2BuildFlat の excludeAtLimit を使う）
%     3. 正規化統計は train の全遷移から算出（val/test のリーク防止。拘束を含む）
%     4. マニフェスト（生成条件・シード・固定姿勢・親資産のコミット）を保存
c = j2_setup_path();
if nargin < 1, opts = struct(); end
def = struct('preset','full','outDir',c.dataDir,'shard',[1 1],'maxScen',Inf,'cacheOnly',false,'save',true);
f = fieldnames(def);
for i = 1:numel(f)
    if ~isfield(opts,f{i}) || isempty(opts.(f{i})), opts.(f{i}) = def.(f{i}); end
end
jp = j2Params();
S = j2ScenarioSet(opts.preset);
cacheDir = fullfile(opts.outDir,'cache');
if ~isfolder(cacheDir), mkdir(cacheDir); end
cacheFn = @(s) fullfile(cacheDir, sprintf('%s_%s.mat', opts.preset, s.name));

% ---- 1. 未生成のシナリオを（自分の分担ぶんだけ）実行 ----
todo = find(arrayfun(@(s) ~isfile(cacheFn(s)), S));
mine = todo(mod(todo-1, opts.shard(2)) == opts.shard(1)-1);
fprintf('=== J2 データ生成 preset=%s: 全 %d 本、未生成 %d 本、このプロセス(shard %d/%d)の担当 %d 本 ===\n', ...
    opts.preset, numel(S), numel(todo), opts.shard(1), opts.shard(2), numel(mine));
if ~isempty(mine)
    workDir = tempname;  mkdir(workDir);                 % Simscape のキャッシュをプロセスごとに分離
    oldDir = cd(workDir);
    cleanWork = onCleanup(@() cleanupWork(oldDir, workDir));   % 終了時に cwd を戻し、作業フォルダを削除
    mE = buildJ2Excitation(struct('modelName','j2_excite','slxPath',fullfile(workDir,'j2_excite.slx')));
    mC = buildJ2ClosedLoop(struct('modelName','j2_closedloop','slxPath',fullfile(workDir,'j2_closedloop.slx')));
    models = struct('excite',mE,'closed',mC);
    nNew = 0;
    for k = mine
        if nNew >= opts.maxScen, break; end
        s = S(k);  t0 = tic;
        raw = runJ2Scenario(s, models);
        e = struct('name',s.name,'q',single(raw.q),'dq',single(raw.dq),'tau',single(raw.tau), ...
            'qref',single(raw.qref),'q0',raw.q0,'dq0',raw.dq0,'contact8k',raw.contact,'fs',jp.fsSim, ...
            'wall',toc(t0)); %#ok<NASGU>
        tmp = [cacheFn(s) '.tmp'];  save(tmp, 'e');  movefile(tmp, cacheFn(s));   % 書き込み途中のファイルを残さない
        nNew = nNew + 1;
        fprintf('  [%3d/%3d] %-18s %-10s %6.1f s 実行 %5.1f s 経過\n', k, numel(S), s.name, s.type, raw.t(end), toc(t0));
    end
    for m = {mE, mC}, if bdIsLoaded(m{1}), close_system(m{1},0); end, end
end

% ---- 2. すべて揃っていれば組み立て ----
ds = [];
if opts.cacheOnly, return; end
if any(arrayfun(@(s) ~isfile(cacheFn(s)), S))
    fprintf('未生成のシナリオが残っています。同じ引数で再実行すると続きから生成します。\n');
    return;
end
ds = assembleJ2Dataset(S, cacheFn, opts, jp);
end

function ds = assembleJ2Dataset(S, cacheFn, opts, jp)
tol = deg2rad(0.1);
scen = repmat(struct('name','','pattern','','type','','phase',[],'split','','seed',0, ...
    'q',[],'dq',[],'tau',[],'tauInst',[],'qref',[],'atLimit',[],'nTransition',0), numel(S), 1);
for k = 1:numel(S)
    s = S(k);  C = load(cacheFn(s));  e = C.e;
    D = j2Downsample(e, jp.fsSim, jp.fsData);
    % 遷移 k: 開始・終了状態のどちらかがリミット拘束
    lim = @(q) q <= jp.qMin + tol | q >= jp.qMax - tol;
    atLimit = lim(D.q(1:end-1)) | lim(D.q(2:end));
    scen(k) = struct('name',s.name,'pattern',s.pattern,'type',s.type,'phase',s.phase,'split',s.split, ...
        'seed',s.seed,'q',D.q,'dq',D.dq,'tau',D.tau,'tauInst',D.tauInst,'qref',D.qref, ...
        'atLimit',atLimit,'nTransition',numel(D.tau));
end
% ---- 正規化統計: train かつ非拘束のみ ----
isTrain = strcmp({scen.split},'train');
X = [];  Y = [];
for k = find(isTrain)
    sc = scen(k);
    X = [X; sc.q(1:end-1), sc.dq(1:end-1), sc.tau]; %#ok<AGROW>
    Y = [Y; diff(sc.q), diff(sc.dq)];                %#ok<AGROW>
end
stats.xMean = mean(X,1);  stats.xStd = std(X,0,1);
stats.yMean = mean(Y,1);  stats.yStd = std(Y,0,1);
stats.xStd(stats.xStd < eps) = 1;  stats.yStd(stats.yStd < eps) = 1;
stats.note = '正規化統計は train の全遷移（リミット拘束を含む）から算出（val/test のリーク防止）。X=[q dq tau], Y=[Δq Δdq]';
% 引継ぎ資料のスケール（仕様どおり）: 実機値で再定義
dt = 1/jp.fsData;
scale.DTHETA_MAX = jp.qdMax;                 scale.TAU_MAX = jp.tauPeak;
scale.D_THETA_MAX = jp.qdMax*dt;             scale.D_DTHETA_MAX = (jp.tauPeak/jp.M)*dt;

manifest.created = char(datetime('now'));
manifest.preset = opts.preset;  manifest.fsRaw = jp.fsSim;  manifest.fsData = jp.fsData;
manifest.nScenario = numel(scen);
manifest.nTrain = nnz(isTrain);  manifest.nVal = nnz(strcmp({scen.split},'val'));
manifest.nTest = nnz(strcmp({scen.split},'test'));  manifest.nBenchmark = nnz(strcmp({scen.split},'benchmark'));
manifest.nTransition = sum([scen.nTransition]);
manifest.qFixedDeg = rad2deg(jp.qFixed);
manifest.limitTolDeg = 0.1;
manifest.downsample = '状態は瞬時値を間引き、トルクは区間平均（j2Downsample）';
manifest.params = jp;
manifest.parentCommit = parentCommit();
manifest.matlab = version;
ds.scen = scen;  ds.stats = stats;  ds.scale = scale;  ds.manifest = manifest;
if opts.save
    if ~isfolder(opts.outDir), mkdir(opts.outDir); end
    fn = fullfile(opts.outDir, sprintf('j2_dataset_%s.mat', opts.preset));
    save(fn, '-struct', 'ds', '-v7.3');
    fprintf('保存: %s（%d シナリオ、%d 遷移、train/val/test/bench = %d/%d/%d/%d）\n', fn, numel(scen), ...
        manifest.nTransition, manifest.nTrain, manifest.nVal, manifest.nTest, manifest.nBenchmark);
end
end

function cleanupWork(oldDir, workDir)
cd(oldDir);
if isfolder(workDir), rmdir(workDir, 's'); end
end

function h = parentCommit()
c = j2_config();
[st, out] = system(sprintf('git -C "%s" rev-parse HEAD 2>/dev/null', c.parentDir));
if st == 0, h = strtrim(out); else, h = 'unknown'; end
end
