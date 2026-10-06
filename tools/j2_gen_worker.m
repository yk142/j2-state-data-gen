function j2_gen_worker(preset, shardI, shardN)
%J2_GEN_WORKER データ生成を分担する作業プロセス用のエントリポイント。
%   別プロセスから実行: matlab -batch "addpath('tools'); j2_gen_worker('full', 1, 2)"
%   shardI/shardN で担当を分ける（Parallel Computing Toolbox が使えないための代替）。
%   生成だけ行い（cacheOnly）、組み立ては全員の生成後に genJ2Dataset を 1 回呼んで行う。
repo = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(repo,'config'), genpath(fullfile(repo,'src')));
if ischar(shardI), shardI = str2double(shardI); end
if ischar(shardN), shardN = str2double(shardN); end
genJ2Dataset(struct('preset',preset,'shard',[shardI shardN],'cacheOnly',true));
end
