function j2_finish(preset, outDir)
%J2_FINISH 全シナリオの生成後に、データセットの組み立て・検証・カバレッジ評価を行う。
%   別プロセスから実行: matlab -batch "addpath('tools'); j2_finish('full')"
%   outDir の既定は reports/issue-9。結果の指標は outDir/metrics.json と data/logs に残す。
repo = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(repo,'config'), genpath(fullfile(repo,'src')));
if nargin < 2, outDir = fullfile(repo,'reports','issue-9'); end
ds = genJ2Dataset(struct('preset',preset));
assert(~isempty(ds), 'j2_finish:incomplete', '未生成のシナリオが残っています');
R = validateJ2Data(ds, true);
M = j2Coverage(ds, outDir);
out.validation = struct('nErrors',R.nErrors,'metrics',R.metrics, ...
    'checks',arrayfun(@(c) struct('name',c.name,'pass',c.pass,'value',c.value), R.checks, 'UniformOutput', false));
out.coverage = M;  out.manifest = ds.manifest;
fid = fopen(fullfile(outDir,'metrics.json'),'w');  fwrite(fid, jsonencode(out, 'PrettyPrint', true));  fclose(fid);
fprintf('完了: 検証 不合格 %d 件、3 次元占有率 %.1f%%、ベンチマークが学習系の占有セルに入る割合 %.1f%%\n', ...
    R.nErrors, 100*M.occupancy3D, 100*M.benchmarkInOccupied3D);
end
