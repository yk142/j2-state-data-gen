function setJ2InitialState(mdl, q0, dq0)
%SETJ2INITIALSTATE J2 単軸プラントの初期角度・初期角速度を設定する。
%   SETJ2INITIALSTATE(mdl, q0)        初期角度 q0 [rad]（初期角速度は 0）
%   SETJ2INITIALSTATE(mdl, q0, dq0)   初期角速度 dq0 [rad/s] も指定
%
%   可動関節 'joint1'（buildJ2Plant が J2 を改名したもの）の Position/Velocity Target を設定する。
if nargin < 3, dq0 = 0; end
blk = [mdl '/joint1'];
set_param(blk, ...
    'PositionTargetSpecify','on', 'PositionTargetValue',sprintf('%.12g',q0), 'PositionTargetValueUnits','rad', ...
    'VelocityTargetSpecify','on', 'VelocityTargetValue',sprintf('%.12g',dq0), 'VelocityTargetValueUnits','rad/s');
end
