function mdl = buildJ2Excitation(opts)
%BUILDJ2EXCITATION 開ループ励振モデル（励振トルク＋ソフトリミットバリア＋トルク飽和）を生成する。
%   mdl = BUILDJ2EXCITATION(opts)
%     opts は buildJ2Plant のものに加え:
%       opts.modelName  既定 'j2_excite'
%       opts.barrierWn  バリアの固有角周波数 [rad/s]（既定 60）
%   モデル I/F:  Inport 'tauEx'[スカラ] → Outport 'q', 'dq', 'tau'（飽和後の実印加トルク）
%
%   ソフトリミットバリア（親資産 buildC8Excitation と同じ考え方）:
%       τ = sat( τ_ex − Kp_bar·deadzone(q; qSoftMin, qSoftMax) − Kd_bar·deadzone(dq; ±dqSoft) )
%   Dead Zone は安全帯の内側で厳密に 0 を返すため、帯の内側では励振の分布を歪めない。
%   安全帯は setJ2Barrier(mdl, qSoftFrac, dqSoftFrac) で実行前に変更できる。
%   ゲインは実効慣性 M から決める: Kp = M·wn², Kd = 2·M·wn（臨界減衰）。
if nargin < 1, opts = struct(); end
if ~isfield(opts,'modelName') || isempty(opts.modelName), opts.modelName = 'j2_excite'; end
if ~isfield(opts,'barrierWn') || isempty(opts.barrierWn), opts.barrierWn = 60; end
bw = opts.barrierWn;  opts = rmfield(opts,'barrierWn');
mdl = buildJ2Plant(opts);
jp = j2Params();

set_param([mdl '/tau'],'Name','tauEx');
delete_line(mdl,'tauEx/1','tau_demux/1');
add_block('simulink/Discontinuities/Dead Zone',[mdl '/qDead'],'LowerValue','-1','UpperValue','1');
add_block('simulink/Discontinuities/Dead Zone',[mdl '/dqDead'],'LowerValue','-1','UpperValue','1');
add_block('simulink/Math Operations/Gain',[mdl '/KpBar'],'Gain',sprintf('-%.10g',jp.M*bw^2));
add_block('simulink/Math Operations/Gain',[mdl '/KdBar'],'Gain',sprintf('-%.10g',2*jp.M*bw));
add_block('simulink/Math Operations/Sum',[mdl '/tauSum'],'Inputs','+++');
add_block('simulink/Discontinuities/Saturation',[mdl '/tauSat'], ...
    'UpperLimit',sprintf('%.10g',jp.tauPeak),'LowerLimit',sprintf('%.10g',-jp.tauPeak));
add_block('simulink/Sinks/Out1',[mdl '/tau']);  set_param([mdl '/tau'],'Port','3');
add_line(mdl,'q_mux/1','qDead/1');   add_line(mdl,'qDead/1','KpBar/1');
add_line(mdl,'dq_mux/1','dqDead/1'); add_line(mdl,'dqDead/1','KdBar/1');
add_line(mdl,'tauEx/1','tauSum/1');  add_line(mdl,'KpBar/1','tauSum/2'); add_line(mdl,'KdBar/1','tauSum/3');
add_line(mdl,'tauSum/1','tauSat/1'); add_line(mdl,'tauSat/1','tau_demux/1'); add_line(mdl,'tauSat/1','tau/1');
setJ2Barrier(mdl, 0.85, 0.60);
save_system(mdl, get_param(mdl,'FileName'));
end
