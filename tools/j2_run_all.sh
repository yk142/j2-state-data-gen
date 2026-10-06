#!/bin/bash
# 未生成のシナリオを 1 プロセスで生成し、続けて組み立て・検証・カバレッジ評価を行う。
# 使い方: tools/j2_run_all.sh <preset> <評価の出力先（リポジトリからの相対パス）>
cd /home/yk/matlab-projects/j2-state-data-gen || exit 1
M=/usr/local/MATLAB/R2025a/bin/matlab
PRESET=${1:-full}; OUT=${2:-reports/issue-9}
echo "[$(date +%T)] 生成開始 preset=$PRESET"
$M -batch "cd('/home/yk/matlab-projects/j2-state-data-gen'); addpath('tools'); j2_gen_worker('$PRESET',1,1)" > data/logs/gen_$PRESET.log 2>&1
echo "[$(date +%T)] 生成終了。組み立て・検証・評価（出力先 $OUT）"
$M -batch "cd('/home/yk/matlab-projects/j2-state-data-gen'); addpath('tools'); j2_finish('$PRESET', fullfile(pwd,'$OUT'))" > data/logs/finish_$PRESET.log 2>&1
echo "[$(date +%T)] すべて終了"
