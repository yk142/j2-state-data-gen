function [qr, dqr, ddqr] = j2SmoothRef(tVec, center, amp, fMax, seed)
%J2SMOOTHREF 中心まわりに滑らかにゆらぐ参照軌道（姿勢保持＋ゆっくりした設定値変動）。
%   [qr, dqr, ddqr] = J2SMOOTHREF(tVec, center, amp, fMax, seed)
%     center 中心角 [rad]、amp 振れ幅の最大値 [rad]、fMax 帯域 [Hz]（1〜2 Hz 程度）、seed 乱数シード
%   BLN を 200 Hz で生成して 3 次スプライン（C2）で 8 kHz に補間し、数値微分で dq, ddq を得る。
%   C2 連続なので ddq_ref が不連続にならず、フィードフォワードが滑らかになる。
tVec = tVec(:);
fsG = 200;
tG = (tVec(1):1/fsG:tVec(end)+1/fsG)';
rs = RandStream('twister','Seed',seed);
[b,a] = butter(4, fMax/(fsG/2));
y = filtfilt(b,a,randn(rs,numel(tG),1));
y = y / max(abs(y));
qr = center + amp*interp1(tG, y, tVec, 'spline');
dt = tVec(2) - tVec(1);
dqr  = gradient(qr, dt);
ddqr = gradient(dqr, dt);
end
