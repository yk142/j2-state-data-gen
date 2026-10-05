function mdl = phase0_fixA(outDir, name)
%PHASE0_FIXA 方式A: 生成 URDF の J2 以外を fixed joint にして J2 単軸プラントを作る。
%   mdl = PHASE0_FIXA(outDir, name)
%   親資産の writeC8URDF で 6 軸 URDF を作り、joint1,3..6 を fixed に後処理して smimport する。
%   固定姿勢は j2_config の qFixed（URDF の origin に Rz(qFixed(i)) を合成した rpy で表現）。
%   親の addC8Actuation/Armature/Friction は joint1..nJoints を前提とするため、
%   nJoints=1 の p を作り、J2 のブロックを 'joint1' に改名して適用する（親コードは変更しない）。
c = j2_setup_path();
p = c8_params(c.meshDir);
if nargin < 2, name = 'j2_fixA'; end
urdf6 = fullfile(outDir,[name '_6ax.urdf']);
urdf1 = fullfile(outDir,[name '.urdf']);
writeC8URDF(p, urdf6);

% ---- URDF 後処理: joint1,3..6 を fixed に（axis/limit を削除）----
txt = fileread(urdf6);
for i = [1 3 4 5 6]
    pat = sprintf('(<joint name="joint%d" type=")revolute(">.*?</joint>)', i);
    tok = regexp(txt, pat, 'tokens', 'once', 'dotexceptnewline');
    blk = regexp(txt, sprintf('<joint name="joint%d" type="revolute">.*?</joint>', i), 'match', 'once');
    assert(~isempty(blk), 'joint%d が見つかりません', i);
    nb = strrep(blk, 'type="revolute"', 'type="fixed"');
    nb = regexprep(nb, '\s*<axis[^>]*/>', '');
    nb = regexprep(nb, '\s*<limit[^>]*/>', '');
    if c.qFixed(i) ~= 0                         % 固定角を origin の回転に畳み込む
        rpy = rotm2rpy(rpy2rotm(p.joint(i).rpy) * rotz(c.qFixed(i)));
        nb = regexprep(nb, 'rpy="[^"]*"', sprintf('rpy="%.12g %.12g %.12g"', rpy), 'once');
    end
    txt = strrep(txt, blk, nb);
end
fid = fopen(urdf1,'w'); fwrite(fid, txt); fclose(fid);

% ---- smimport ----
mdl = name;
if bdIsLoaded(mdl), close_system(mdl,0); end
smimport(urdf1, 'ModelName', mdl);
mc = find_system(mdl,'MaskType','Mechanism Configuration');
set_param(mc{1}, 'GravityVector', '[0 0 -9.80665]');
set_param(mdl,'SolverType','Fixed-step','Solver','ode14x','FixedStep',num2str(1/c.fsSim,'%.12g'), ...
    'StopTime','0.5','SimscapeLogType','none');
try, set_param(mdl,'SimMechanicsOpenEditorOnUpdate','off'); catch, end

% ---- J2 を 'joint1' に改名（親の add* が joint1..n を参照するため）----
weld1 = find_system(mdl,'SearchDepth',1,'Name','joint1');
if ~isempty(weld1), set_param([mdl '/joint1'],'Name','weldJ1'); end
set_param([mdl '/joint2'],'Name','joint1');

% ---- nJoints=1 の p（J2 の値）----
p1 = p;  p1.nJoints = 1;
p1.armature = p.armature(2);   p1.frictionFc = p.frictionFc(2);
p1.frictionBv = p.frictionBv(2); p1.frictionEps = p.frictionEps(2);
addC8Actuation(mdl, p1);
addC8Armature(mdl, p1);
addC8Friction(mdl, p1);
save_system(mdl, fullfile(outDir,[name '.slx']));
fprintf('phase0_fixA: 保存 %s (revolute=%d, weld=%d)\n', fullfile(outDir,[name '.slx']), ...
    numel(find_system(mdl,'MaskType','Revolute Joint')), numel(find_system(mdl,'MaskType','Weld Joint')));
end

% ---- 回転ヘルパ（URDF の rpy は R = Rz(y)*Ry(p)*Rx(r)）----
function R = rotz(a), R = [cos(a) -sin(a) 0; sin(a) cos(a) 0; 0 0 1]; end
function R = rpy2rotm(v)
r = v(1); p = v(2); y = v(3);
Rx = [1 0 0; 0 cos(r) -sin(r); 0 sin(r) cos(r)];
Ry = [cos(p) 0 sin(p); 0 1 0; -sin(p) 0 cos(p)];
R = rotz(y) * Ry * Rx;
end
function rpy = rotm2rpy(R)
p = -asin(max(-1,min(1,R(3,1))));
if abs(cos(p)) > 1e-9
    r = atan2(R(3,2), R(3,3));  y = atan2(R(2,1), R(1,1));
else                              % ジンバルロック
    r = 0;  y = atan2(-R(1,2), R(2,2));
end
rpy = [r p y];
end
