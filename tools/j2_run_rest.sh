#!/bin/bash
# 1 本目のワーカーの終了を待ち、残り（2 本目の担当分）の生成 → 組み立て・検証・カバレッジ評価を順に実行する。
# 使い方: tools/j2_run_rest.sh <1 本目のワーカーの PID>
cd /home/yk/matlab-projects/j2-state-data-gen || exit 1
M=/usr/local/MATLAB/R2025a/bin/matlab
while kill -0 "$1" 2>/dev/null; do sleep 30; done
echo "[$(date +%T)] worker1 終了。worker2 開始"
$M -batch "cd('/home/yk/matlab-projects/j2-state-data-gen'); addpath('tools'); j2_gen_worker('full',2,2)" > data/logs/worker2.log 2>&1
echo "[$(date +%T)] worker2 終了。組み立て・検証・カバレッジ評価を開始"
$M -batch "cd('/home/yk/matlab-projects/j2-state-data-gen'); addpath('tools'); j2_finish('full')" > data/logs/finish.log 2>&1
echo "[$(date +%T)] すべて終了"
