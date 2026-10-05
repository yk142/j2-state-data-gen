function [qr, dqr, ddqr] = j2PtpTrajectory(t, WP, T)
%J2PTPTRAJECTORY J2 の PTP 目標軌道（5 次多項式・マルチウェイポイント、ベクトル化）。
%   [qr, dqr, ddqr] = J2PTPTRAJECTORY(t, WP, T)
%     t   時刻 [N x 1] [s]、WP ウェイポイント [(nSeg+1) x 1] [rad]、T 各区間の所要時間 [nSeg x 1] [s]
%   区間内は s(τ) = 10τ³ − 15τ⁴ + 6τ⁵（τ=(t−t0)/T）で補間し、各点で速度・加速度が 0 になる
%   （親資産 c8TrajGen と同じ。J2 単軸用にベクトル化）。t<0 は開始点、t>ΣT は最終点を保持する。
t = t(:);  WP = WP(:);  T = T(:);
tEnd = [0; cumsum(T)];                    % 各区間の開始時刻 + 終了時刻
qr = repmat(WP(1), size(t));  dqr = zeros(size(t));  ddqr = zeros(size(t));
for k = 1:numel(T)
    idx = t > tEnd(k) & t <= tEnd(k+1);
    tau = (t(idx) - tEnd(k)) / T(k);
    d = WP(k+1) - WP(k);
    qr(idx)   = WP(k) + d*(10*tau.^3 - 15*tau.^4 + 6*tau.^5);
    dqr(idx)  = d*(30*tau.^2 - 60*tau.^3 + 30*tau.^4) / T(k);
    ddqr(idx) = d*(60*tau - 180*tau.^2 + 120*tau.^3) / T(k)^2;
end
qr(t > tEnd(end)) = WP(end);
end
