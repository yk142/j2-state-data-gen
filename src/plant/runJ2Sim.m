function out = runJ2Sim(mdl, tVec, U, q0, dq0)
%RUNJ2SIM 外部入力を与えて J2 モデル（プラント／励振／閉ループ）を実行する。
%   out = RUNJ2SIM(mdl, tVec, U, q0, dq0)
%     mdl   ロード済みモデル名（buildJ2Plant / buildJ2Excitation / buildJ2ClosedLoop の生成物）
%     tVec  時刻列 [N x 1] [s]（先頭 0。ソルバ刻み 8 kHz の格子に載せること）
%     U     入力列 [N x nIn]（モデルの Inport 幅に一致）
%     q0, dq0  初期角度 [rad]・角速度 [rad/s]
%   out: .t .q .dq（スカラ状態）、.tau（3 番目の Outport があれば実印加トルク、無ければ空）
if nargin < 5, dq0 = 0; end
tVec = tVec(:);
assert(size(U,1) == numel(tVec), 'runJ2Sim:size', 'tVec と U の行数が一致しません');
if dq0 ~= 0
    % 初期角速度が非ゼロのときの armature 回転子初期条件の警告（simJ2 と同じ理由で抑制）
    wst = warning('off','physmod:simscape:engine:core:dae_errors:IcDiagnostics');
    restoreWarn = onCleanup(@() warning(wst));
end
setJ2InitialState(mdl, q0, dq0);
assignin('base','J2_in',[tVec U]);
o = sim(mdl, 'StopTime', num2str(tVec(end),'%.12g'), 'LoadExternalInput','on', ...
    'ExternalInput','J2_in', 'SaveOutput','on', 'SaveFormat','Dataset');
y = o.yout;
out.t  = y.getElement(1).Values.Time;
out.q  = y.getElement(1).Values.Data(:,1);
out.dq = y.getElement(2).Values.Data(:,1);
out.tau = [];
if y.numElements >= 3
    tauTS = y.getElement(3).Values;
    out.tau = tauTS.Data(:,1);
    if numel(out.tau) ~= numel(out.t)
        out.tau = interp1(tauTS.Time, out.tau, out.t, 'previous', 'extrap');
    end
end
end
