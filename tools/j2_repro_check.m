function R = j2_repro_check(preset, names, outDir)
%J2_REPRO_CHECK キャッシュを使わずにシナリオを再生成し、保存済みの生データと一致するか確認する。
%   R = J2_REPRO_CHECK(preset, names, outDir)
%     preset  'full'（既定）、names 再生成するシナリオ名の cell、outDir 結果の保存先（既定 reports/issue-13）
%   別プロセス（新しい MATLAB）で実行すると、セッションの状態に依存しないことも確認できる:
%     matlab -batch "addpath('tools'); j2_repro_check()"
%   比較: 再生成した 8 kHz の q, dq, tau を single に丸めたものと、キャッシュ（single）の差の最大値。
%   決定的な実行なら差は 0（ビット一致）になる。
repo = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(repo,'config'), genpath(fullfile(repo,'src')));
c = j2_config();
if nargin < 1 || isempty(preset), preset = 'full'; end
if nargin < 2 || isempty(names)
    names = {'bln_high_001','freefall_001','step_001','micro_001','hold_-090_001','nearlimit_001', ...
        'ptp_v40_001','contact_fall_001','contact_drive_001','contact_bln_001'};
end
if nargin < 3 || isempty(outDir), outDir = fullfile(repo,'reports','issue-13'); end
if ~isfolder(outDir), mkdir(outDir); end
S = j2ScenarioSet(preset);
workDir = tempname;  mkdir(workDir);
oldDir = cd(workDir);
cleanWork = onCleanup(@() cleanup(oldDir, workDir));
mE = buildJ2Excitation(struct('modelName','j2_excite','slxPath',fullfile(workDir,'j2_excite.slx')));
mC = buildJ2ClosedLoop(struct('modelName','j2_closedloop','slxPath',fullfile(workDir,'j2_closedloop.slx')));
models = struct('excite',mE,'closed',mC);
R = struct('name',{},'type',{},'nSample',{},'maxAbsDiff_q',{},'maxAbsDiff_dq',{},'maxAbsDiff_tau',{},'bitIdentical',{},'wall',{});
for i = 1:numel(names)
    s = S(strcmp({S.name}, names{i}));
    assert(isscalar(s), 'j2_repro_check:name', 'シナリオが見つかりません: %s', names{i});
    t0 = tic;  raw = runJ2Scenario(s, models);  wall = toc(t0);
    C = load(fullfile(c.dataDir,'cache',sprintf('%s_%s.mat', preset, s.name)));  e = C.e;
    assert(numel(e.q) == numel(raw.q), 'j2_repro_check:len', '%s: 長さが違います（%d vs %d）', names{i}, numel(e.q), numel(raw.q));
    dq_ = max(abs(double(e.q)   - double(single(raw.q))));
    dd_ = max(abs(double(e.dq)  - double(single(raw.dq))));
    dt_ = max(abs(double(e.tau) - double(single(raw.tau))));
    R(end+1) = struct('name',s.name,'type',s.type,'nSample',numel(raw.q),'maxAbsDiff_q',dq_, ...
        'maxAbsDiff_dq',dd_,'maxAbsDiff_tau',dt_,'bitIdentical',(dq_ == 0 && dd_ == 0 && dt_ == 0),'wall',wall); %#ok<AGROW>
    fprintf('  %-20s %-10s %7d 点  max|Δq|=%.3g  max|Δdq|=%.3g  max|Δtau|=%.3g  %s  (%.0f s)\n', s.name, s.type, ...
        numel(raw.q), dq_, dd_, dt_, ternary(R(end).bitIdentical,'ビット一致','不一致'), wall);
end
for m = {mE, mC}, if bdIsLoaded(m{1}), close_system(m{1},0); end, end
out = struct('created', char(datetime('now')), 'preset', preset, 'matlab', version, 'results', R);
fid = fopen(fullfile(outDir,'repro.json'),'w');  fwrite(fid, jsonencode(out,'PrettyPrint',true));  fclose(fid);
fprintf('再現確認: %d 本中 %d 本がビット一致\n', numel(R), nnz([R.bitIdentical]));
end

function cleanup(oldDir, workDir)
cd(oldDir);
if isfolder(workDir), rmdir(workDir, 's'); end
end

function s = ternary(c, a, b)
if c, s = a; else, s = b; end
end
