function R = phase0_holdB(srcSlx, wnHold, q2deg, T)
%PHASE0_HOLDB 方式B: 6軸プラントの J2 以外を PD で位置保持し、J2 単軸として使えるか試す。
%   R = PHASE0_HOLDB(srcMdl, wnHold, q2deg, T)
%     srcSlx  buildC8Plant 生成済みモデルの .slx パス（tau[6] -> q[6], dq[6]）
%     wnHold  保持ゲイン設計の固有角周波数 [rad/s]（c8GainDesign）
%     q2deg   J2 初期角 [deg]、T 実行時間 [s]
%   J2 には重力補償トルク c8RNEA(q,0,0) の J2 成分を与え、他軸の変位と J2 のドリフトを測る。
c = j2_setup_path();
p = c8_params(c.meshDir);
g = c8GainDesign(p, struct('wn', wnHold));

[d,srcName] = fileparts(srcSlx);
mdl = [srcName '_holdB'];
if bdIsLoaded(mdl), close_system(mdl,0); end
if bdIsLoaded(srcName), close_system(srcName,0); end
load_system(srcSlx);
save_system(srcName, fullfile(d,[mdl '.slx']));   % 別名保存（親モデルは変更しない）

% 既存の tau 入力 → demux の直結を外し、(保持PD + J2 トルク) に置換
delete_line(mdl,'tau/1','tau_demux/1');
set_param([mdl '/tau'],'PortDimensions','1');     % J2 トルク（スカラ）
Kp = g.Kp(:); Kd = g.Kd(:); Kp(2) = 0; Kd(2) = 0;      % J2 は保持しない（他軸のみ PD で保持）
assignin('base','KpHold_j2', Kp);  assignin('base','KdHold_j2', Kd);
assignin('base','e2_j2', [0;1;0;0;0;0]);
assignin('base','qFix_j2', zeros(6,1));
add_block('simulink/Math Operations/Gain',[mdl '/e2'],'Gain','e2_j2','Multiplication','Matrix(K*u)');
add_block('simulink/Math Operations/Sum',[mdl '/eq'],'Inputs','+-');
add_block('simulink/Sources/Constant',[mdl '/qFix'],'Value','qFix_j2');
add_block('simulink/Math Operations/Gain',[mdl '/KpH'],'Gain','-KpHold_j2','Multiplication','Element-wise(K.*u)');
add_block('simulink/Math Operations/Gain',[mdl '/KdH'],'Gain','-KdHold_j2','Multiplication','Element-wise(K.*u)');
add_block('simulink/Math Operations/Sum',[mdl '/tauSum'],'Inputs','+++');
add_line(mdl,'q_mux/1','eq/1'); add_line(mdl,'qFix/1','eq/2');
add_line(mdl,'eq/1','KpH/1'); add_line(mdl,'dq_mux/1','KdH/1');
add_line(mdl,'tau/1','e2/1');
add_line(mdl,'e2/1','tauSum/1'); add_line(mdl,'KpH/1','tauSum/2'); add_line(mdl,'KdH/1','tauSum/3');
add_line(mdl,'tauSum/1','tau_demux/1');

% 初期姿勢（J2 のみ q2deg、他は 0）
for i = 1:6
    qi = 0; if i==2, qi = deg2rad(q2deg); end
    set_param(sprintf('%s/joint%d',mdl,i),'PositionTargetSpecify','on', ...
        'PositionTargetValue',sprintf('%.12g',qi),'PositionTargetValueUnits','rad');
end
% 重力補償（J2、他軸 0）。q2 固定で静止させる
qv = zeros(1,6); qv(2) = deg2rad(q2deg);
tauG = c8RNEA(qv, zeros(1,6), zeros(1,6), p);
tauG2 = tauG(2);
assignin('base','tauJ2_in', [0 tauG2; T tauG2]);

t0 = tic;
o = sim(mdl,'StopTime',num2str(T),'LoadExternalInput','on','ExternalInput','tauJ2_in', ...
    'SaveOutput','on','SaveFormat','Dataset');
R.wall = toc(t0);
t = o.yout.getElement(1).Values.Time;  q = o.yout.getElement(1).Values.Data;  dq = o.yout.getElement(2).Values.Data;
R.t = t; R.q = q; R.dq = dq; R.tauG2 = tauG2;
R.maxOtherDev_deg = rad2deg(max(abs(q(:,[1 3:6])),[],1));
R.j2Drift_deg = rad2deg(q(end,2) - q(1,2));
R.maxAbsDq = max(abs(dq),[],1);
R.finite = all(isfinite(q(:))) && all(isfinite(dq(:)));
fprintf('wall=%.1fs (sim %.2fs) finite=%d J2 drift=%.3f deg tauG2=%.1f N*m\n', R.wall, T, R.finite, R.j2Drift_deg, tauG2);
fprintf('他軸最大変位[deg] J1,J3..J6 = %s\n', mat2str(R.maxOtherDev_deg,3));
end
