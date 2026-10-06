function info = j2ExportFlat(ds, outFile, opts)
%J2EXPORTFLAT データセットを、MATLAB 以外（Python 等）から読めるフラットな形式で書き出す。
%   info = J2EXPORTFLAT(ds, outFile, opts)
%     ds       genJ2Dataset の出力（または data/j2_dataset_<preset>.mat を load したもの）
%     outFile  出力先。拡張子で形式が決まる: .h5（HDF5、既定）/ .mat（v7。scipy.io.loadmat で読める）
%     opts.splits          書き出す split（既定 {'train','val','test','benchmark'}）
%     opts.excludeAtLimit  リミット拘束の遷移を除く（既定 false。j2BuildFlat と同じ）
%     opts.dtype           行列の型 'single'（既定）| 'double'。単精度でも Δq は倍精度で計算してから丸める
%     opts.deflate         HDF5 の圧縮レベル 0〜9（既定 4）
%   info: 出力ファイル、サイズ、split ごとの遷移数
%
%   HDF5 の構成（Python/h5py から見た形。MATLAB は列優先なので転置して書いている）:
%     /<split>/X          (N, 3)  [q rad, dq rad/s, tau N·m]            ← 入力（SI 単位）
%     /<split>/Y          (N, 2)  [Δq rad, Δdq rad/s]（1 ステップ = 1 ms）← 目標
%     /<split>/X_spec     (N, 4)  [sin q, cos q, dq/DTHETA_MAX, tau/TAU_MAX]  ← 引継ぎ資料の正規化形式
%     /<split>/Y_spec     (N, 2)  [Δq/D_THETA_MAX, Δdq/D_DTHETA_MAX]
%     /<split>/scenario_id (N,)   int32、/scenarios/* の添字（0 始まり）
%     /<split>/step        (N,)   int32、シナリオ内の遷移番号（0 始まり）
%     /<split>/at_limit    (N,)   uint8、リミット拘束フラグ（遷移の開始・終了状態がリミットから 0.1° 以内または外側）
%     /scenarios/{name,pattern,type,phase,split}  文字列、/scenarios/{seed,n_transition}  整数
%     /stats/{x_mean,x_std,y_mean,y_std}   正規化統計（train の全遷移）、/scale/*  引継ぎ資料のスケール
%     ルート属性: preset, created, fs_data, q_fixed_deg, parent_commit, matlab, layout など
%   インデックスは他言語向けに 0 始まり。MATLAB で使うときは +1 する。
if nargin < 3, opts = struct(); end
def = struct('splits',{{'train','val','test','benchmark'}}, 'excludeAtLimit',false, 'dtype','single', 'deflate',4);
f = fieldnames(def);
for i = 1:numel(f)
    if ~isfield(opts,f{i}) || isempty(opts.(f{i})), opts.(f{i}) = def.(f{i}); end
end
[outDir,~,ext] = fileparts(outFile);
if ~isempty(outDir) && ~isfolder(outDir), mkdir(outDir); end
fo = struct('excludeAtLimit', opts.excludeAtLimit);
S = ds.scen;
off = [0; cumsum([S.nTransition]')];            % 全シナリオの at_limit を連結したときの先頭位置
atAll = vertcat(S.atLimit);

% ---- split ごとの行列を作る ----
D = struct();
for i = 1:numel(opts.splits)
    sp = opts.splits{i};
    [X, Y, inf1] = j2BuildFlat(ds, sp, fo);
    [Xs, Ys] = j2BuildFlat(ds, sp, struct('excludeAtLimit',opts.excludeAtLimit,'features','spec'));
    D.(sp).X = X;  D.(sp).Y = Y;  D.(sp).Xs = Xs;  D.(sp).Ys = Ys;
    D.(sp).scenario_id = int32(inf1.scenario - 1);
    D.(sp).step        = int32(inf1.k - 1);
    D.(sp).at_limit    = uint8(atAll(off(inf1.scenario) + inf1.k));
end
cast = @(A) cast_(A, opts.dtype);

switch lower(ext)
    case '.h5'
        if isfile(outFile), delete(outFile); end
        for i = 1:numel(opts.splits)
            sp = opts.splits{i};  d = D.(sp);
            writeMat(outFile, ['/' sp '/X'],      cast(d.X),  opts, 'q [rad], dq [rad/s], tau [N*m]');
            writeMat(outFile, ['/' sp '/Y'],      cast(d.Y),  opts, 'dq_q [rad] (q_next - q), d_dq [rad/s] (dq_next - dq), per 1 ms step');
            writeMat(outFile, ['/' sp '/X_spec'], cast(d.Xs), opts, 'sin q, cos q, dq/DTHETA_MAX, tau/TAU_MAX');
            writeMat(outFile, ['/' sp '/Y_spec'], cast(d.Ys), opts, 'd_q/D_THETA_MAX, d_dq/D_DTHETA_MAX');
            writeVec(outFile, ['/' sp '/scenario_id'], d.scenario_id, opts, '0-based index into /scenarios');
            writeVec(outFile, ['/' sp '/step'],        d.step,        opts, '0-based transition index within the scenario');
            writeVec(outFile, ['/' sp '/at_limit'],    d.at_limit,    opts, '1 = limit-bound transition');
        end
        writeStr(outFile, '/scenarios/name',    string({S.name}));
        writeStr(outFile, '/scenarios/pattern', string({S.pattern}));
        writeStr(outFile, '/scenarios/type',    string({S.type}));
        writeStr(outFile, '/scenarios/phase',   string(cellfun(@num2str, {S.phase}, 'UniformOutput', false)));
        writeStr(outFile, '/scenarios/split',   string({S.split}));
        writeVec(outFile, '/scenarios/seed',         int64([S.seed]'), opts, '');
        writeVec(outFile, '/scenarios/n_transition', int32([S.nTransition]'), opts, '');
        st = ds.stats;
        writeVec(outFile, '/stats/x_mean', double(st.xMean(:)), opts, 'q, dq, tau (train, all transitions)');
        writeVec(outFile, '/stats/x_std',  double(st.xStd(:)),  opts, '');
        writeVec(outFile, '/stats/y_mean', double(st.yMean(:)), opts, 'd_q, d_dq');
        writeVec(outFile, '/stats/y_std',  double(st.yStd(:)),  opts, '');
        sc = ds.scale;
        for nm = fieldnames(sc)'
            writeVec(outFile, ['/scale/' nm{1}], double(sc.(nm{1})), opts, '');
        end
        writeAttrs(outFile, ds, opts);
    case '.mat'
        out = struct();
        for i = 1:numel(opts.splits)
            sp = opts.splits{i};  d = D.(sp);
            out.([sp '_X']) = cast(d.X);   out.([sp '_Y']) = cast(d.Y);
            out.([sp '_X_spec']) = cast(d.Xs);  out.([sp '_Y_spec']) = cast(d.Ys);
            out.([sp '_scenario_id']) = d.scenario_id;  out.([sp '_step']) = d.step;  out.([sp '_at_limit']) = d.at_limit;
        end
        out.scenario_name = {S.name}';  out.scenario_pattern = {S.pattern}';  out.scenario_type = {S.type}';
        out.scenario_phase = cellfun(@num2str, {S.phase}', 'UniformOutput', false);  out.scenario_split = {S.split}';
        out.scenario_seed = int64([S.seed]');  out.scenario_n_transition = int32([S.nTransition]');
        out.stats_x_mean = ds.stats.xMean(:);  out.stats_x_std = ds.stats.xStd(:);
        out.stats_y_mean = ds.stats.yMean(:);  out.stats_y_std = ds.stats.yStd(:);
        for nm = fieldnames(ds.scale)'
            out.(['scale_' nm{1}]) = ds.scale.(nm{1});
        end
        if isfield(ds,'manifest') && isfield(ds.manifest,'preset')
            out.preset = ds.manifest.preset;  out.fs_data = ds.manifest.fsData;
            out.q_fixed_deg = ds.manifest.qFixedDeg(:);  out.parent_commit = ds.manifest.parentCommit;
        end
        out.layout = 'rows = transitions. indices are 0-based. X=[q dq tau], Y=[dq_q d_dq]';
        save(outFile, '-struct', 'out', '-v7');
    otherwise
        error('j2ExportFlat:ext', '未対応の拡張子: %s（.h5 または .mat）', ext);
end
d = dir(outFile);
info.file = outFile;  info.bytes = d.bytes;  info.dtype = opts.dtype;
for i = 1:numel(opts.splits), info.n.(opts.splits{i}) = size(D.(opts.splits{i}).X, 1); end
end

% ======================= ヘルパ =======================
function A = cast_(A, dtype)
if strcmp(dtype,'single'), A = single(A); else, A = double(A); end
end

function writeMat(file, path, A, opts, columns)
% A は (N, d)。MATLAB は列優先なので (d, N) に転置して書く → h5py からは (N, d) に見える
[N, d] = size(A);
if N == 0, error('j2ExportFlat:empty', '%s は空です', path); end
h5create(file, path, [d N], 'Datatype', class(A), 'ChunkSize', [d min(N, 65536)], 'Deflate', opts.deflate);
h5write(file, path, A');
if ~isempty(columns), h5writeatt(file, path, 'columns', columns); end
end

function writeVec(file, path, v, opts, note)
v = v(:);
if isempty(v), return; end
if numel(v) > 1 && opts.deflate > 0
    h5create(file, path, numel(v), 'Datatype', class(v), 'ChunkSize', min(numel(v), 65536), 'Deflate', opts.deflate);
else
    h5create(file, path, numel(v), 'Datatype', class(v));
end
h5write(file, path, v);
if ~isempty(note), h5writeatt(file, path, 'note', note); end
end

function writeStr(file, path, s)
h5create(file, path, numel(s), 'Datatype', 'string');
h5write(file, path, s(:));
end

function writeAttrs(file, ds, opts)
m = ds.manifest;
put = @(k, v) h5writeatt(file, '/', k, v);
if isfield(m,'preset'), put('preset', char(m.preset)); end
if isfield(m,'created'), put('created', char(m.created)); end
if isfield(m,'fsData'), put('fs_data_hz', double(m.fsData)); end
if isfield(m,'fsRaw'), put('fs_raw_hz', double(m.fsRaw)); end
if isfield(m,'qFixedDeg'), put('q_fixed_deg', double(m.qFixedDeg(:))); end
if isfield(m,'parentCommit'), put('parent_commit', char(m.parentCommit)); end
if isfield(m,'matlab'), put('matlab', char(m.matlab)); end
put('exclude_at_limit', double(opts.excludeAtLimit));
put('dtype', opts.dtype);
put('layout', 'row = transition; HDF5 datasets are (N, d) as seen from h5py. indices are 0-based.');
put('X_columns', 'q [rad], dq [rad/s], tau [N*m]');
put('Y_columns', 'q_next-q [rad], dq_next-dq [rad/s] (1 ms step)');
put('description', 'EPSON C8-A901S J2 single-axis (J3 fixed at 75.0684 deg) one-step transition dataset, 1 kHz');
end
