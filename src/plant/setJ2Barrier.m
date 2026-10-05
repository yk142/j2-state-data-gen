function [qSoftMin, qSoftMax, dqSoft] = setJ2Barrier(mdl, qSoftFrac, dqSoftFrac)
%SETJ2BARRIER 励振モデルのソフトリミットバリアの安全帯を設定する。
%   [qSoftMin, qSoftMax, dqSoft] = SETJ2BARRIER(mdl, qSoftFrac, dqSoftFrac)
%     qSoftFrac   可動域に対する安全帯の比率（可動域の中心から ±frac·半幅。1.0 でリミット位置）
%     dqSoftFrac  最大角速度に対する速度安全帯の比率
jp = j2Params();
mid  = (jp.qMax + jp.qMin)/2;   half = (jp.qMax - jp.qMin)/2 * qSoftFrac;
qSoftMin = mid - half;  qSoftMax = mid + half;  dqSoft = jp.qdMax * dqSoftFrac;
set_param([mdl '/qDead'], 'LowerValue',sprintf('%.12g',qSoftMin), 'UpperValue',sprintf('%.12g',qSoftMax));
set_param([mdl '/dqDead'],'LowerValue',sprintf('%.12g',-dqSoft),  'UpperValue',sprintf('%.12g',dqSoft));
end
