function T = j2PtpSegTime(WP, speedScale, accScale)
%J2PTPSEGTIME 速度・加速度制限から各区間の所要時間を決める。
%   T = J2PTPSEGTIME(WP, speedScale, accScale)
%     speedScale  最大角速度に対する使用率（既定 0.3）
%     accScale    加速度上限に対する使用率（既定 0.5）。加速度上限は (tauPeak − mgL − Fc)/M。
%   5 次多項式のピークは max|dq| = 1.875·Δq/T、max|ddq| = 5.7735·Δq/T²。
%   両方を満たす最小の T を採る（下限 0.05 s）。
if nargin < 2 || isempty(speedScale), speedScale = 0.3; end
if nargin < 3 || isempty(accScale),   accScale   = 0.5; end
jp = j2Params();
ddqMax = (jp.tauPeak - jp.mgL - jp.Fc) / jp.M;
d = abs(diff(WP(:)));
Tv = 1.875 * d / (speedScale * jp.qdMax);
Ta = sqrt(5.7735 * d / (accScale * ddqMax));
T = max([Tv, Ta, 0.05*ones(size(d))], [], 2);
end
