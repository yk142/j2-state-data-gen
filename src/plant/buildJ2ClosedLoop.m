function mdl = buildJ2ClosedLoop(opts)
%BUILDJ2CLOSEDLOOP J2 の PD＋フィードフォワード閉ループモデル（トルク飽和付き）を生成する。
%   mdl = BUILDJ2CLOSEDLOOP(opts)
%     opts は buildJ2Plant のものに加え:
%       opts.modelName  既定 'j2_closedloop'
%       opts.wn         PD の固有角周波数 [rad/s]（既定 30）、opts.zeta 減衰比（既定 1）
%   モデル I/F:
%     Inport 'ref'[4] = [q_ref, dq_ref, tau_ff, tau_dist]
%     Outport 'q', 'dq', 'tau'（飽和後の実印加トルク）
%     τ = sat( tau_ff + tau_dist + Kp·(q_ref − q) + Kd·(dq_ref − dq) ),  ±tauPeak で飽和
%   Kp = M·wn², Kd = 2·ζ·M·wn（M は実効慣性）。PD は連続形で、ソルバ刻み 8 kHz で評価される。
%   PTP 軌道追従（tau_dist=0）と、姿勢保持＋外乱トルク（tau_dist≠0）の両方に使う。
if nargin < 1, opts = struct(); end
if ~isfield(opts,'modelName') || isempty(opts.modelName), opts.modelName = 'j2_closedloop'; end
wn = 30;   if isfield(opts,'wn') && ~isempty(opts.wn), wn = opts.wn; end
zeta = 1;  if isfield(opts,'zeta') && ~isempty(opts.zeta), zeta = opts.zeta; end
opts = rmfield_safe(opts, {'wn','zeta'});
mdl = buildJ2Plant(opts);
jp = j2Params();

set_param([mdl '/tau'],'Name','ref');
set_param([mdl '/ref'],'PortDimensions','4');
delete_line(mdl,'ref/1','tau_demux/1');
add_block('simulink/Signal Routing/Demux',[mdl '/refDemux'],'Outputs','4');
add_block('simulink/Math Operations/Sum',[mdl '/eq'],'Inputs','+-');
add_block('simulink/Math Operations/Sum',[mdl '/edq'],'Inputs','+-');
add_block('simulink/Math Operations/Gain',[mdl '/Kp'],'Gain',sprintf('%.10g',jp.M*wn^2));
add_block('simulink/Math Operations/Gain',[mdl '/Kd'],'Gain',sprintf('%.10g',2*zeta*jp.M*wn));
add_block('simulink/Math Operations/Sum',[mdl '/tauSum'],'Inputs','++++');
add_block('simulink/Discontinuities/Saturation',[mdl '/tauSat'], ...
    'UpperLimit',sprintf('%.10g',jp.tauPeak),'LowerLimit',sprintf('%.10g',-jp.tauPeak));
add_block('simulink/Sinks/Out1',[mdl '/tau']);  set_param([mdl '/tau'],'Port','3');
add_line(mdl,'ref/1','refDemux/1');
add_line(mdl,'refDemux/1','eq/1');   add_line(mdl,'q_mux/1','eq/2');
add_line(mdl,'refDemux/2','edq/1');  add_line(mdl,'dq_mux/1','edq/2');
add_line(mdl,'eq/1','Kp/1');         add_line(mdl,'edq/1','Kd/1');
add_line(mdl,'Kp/1','tauSum/1');     add_line(mdl,'Kd/1','tauSum/2');
add_line(mdl,'refDemux/3','tauSum/3'); add_line(mdl,'refDemux/4','tauSum/4');
add_line(mdl,'tauSum/1','tauSat/1'); add_line(mdl,'tauSat/1','tau_demux/1'); add_line(mdl,'tauSat/1','tau/1');
save_system(mdl, get_param(mdl,'FileName'));
end

function s = rmfield_safe(s, names)
for i = 1:numel(names)
    if isfield(s, names{i}), s = rmfield(s, names{i}); end
end
end
