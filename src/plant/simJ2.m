function out = simJ2(mdl, tVec, tau, q0, dq0)
%SIMJ2 J2 単軸プラントを指定トルク列で実行し、状態の時系列を返す。
%   out = SIMJ2(mdl, tVec, tau, q0, dq0)
%     mdl   buildJ2Plant で生成済み（ロード済み）のモデル名
%     tVec  時刻列 [N x 1] [s]（先頭は 0）
%     tau   トルク列 [N x 1] [N*m]、またはスカラ（一定トルク）
%     q0    初期角度 [rad]、dq0 初期角速度 [rad/s]（既定 0）
%   out: .t [M x 1], .q, .dq, .tau（実印加トルクのログではなく入力列を 8 kHz に保持補間したもの）
%
%   トルクは tVec 上の区間定数入力（ゼロ次保持）として与える。tVec が粗くても
%   ソルバの 8 kHz 刻みで入力が保持される。
if nargin < 5, dq0 = 0; end
tVec = tVec(:);
if isscalar(tau), tau = repmat(tau, size(tVec)); end
tau = tau(:);
assert(numel(tau) == numel(tVec), 'simJ2:size', 'tVec と tau の長さが一致しません');

% 初期角速度が非ゼロのとき、Simscape が armature の回転子速度（RMI で関節速度に拘束）の
% 初期条件を緩和して解き直し、IcDiagnostics 警告を出す。関節速度ターゲット自体は
% 満たされる（test_j2_plant の initialVelocityIsApplied で確認）ため、この警告だけ抑制する。
if dq0 ~= 0
    wst = warning('off','physmod:simscape:engine:core:dae_errors:IcDiagnostics');
    restoreWarn = onCleanup(@() warning(wst));
end
setJ2InitialState(mdl, q0, dq0);
% ゼロ次保持を厳密にするため、各区間の終端直前に同値の点を挿入する
tt = reshape([tVec(1:end-1) tVec(2:end)-eps(tVec(2:end))]', [], 1);
uu = reshape([tau(1:end-1) tau(1:end-1)]', [], 1);
tt = [tt; tVec(end)];  uu = [uu; tau(end)];
assignin('base','J2_tau_in',[tt uu]);

o = sim(mdl, 'StopTime', num2str(tVec(end),'%.12g'), 'LoadExternalInput','on', ...
    'ExternalInput','J2_tau_in', 'SaveOutput','on', 'SaveFormat','Dataset');
out.t   = o.yout.getElement(1).Values.Time;
out.q   = o.yout.getElement(1).Values.Data(:,1);
out.dq  = o.yout.getElement(2).Values.Data(:,1);
out.tau = interp1(tVec, tau, out.t, 'previous', 'extrap');
end
