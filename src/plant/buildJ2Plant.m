function mdl = buildJ2Plant(opts)
%BUILDJ2PLANT J2 以外を固定した C8-A901S の J2 単軸プラント（Simscape Multibody）を生成する。
%   mdl = BUILDJ2PLANT()       既定設定で生成し data/models/j2_plant.slx に保存する。
%   mdl = BUILDJ2PLANT(opts)   構造体 opts で設定を上書きする。
%     opts.modelName  モデル名                 (既定 'j2_plant')
%     opts.slxPath    保存先                   (既定 <dataDir>/models/j2_plant.slx)
%     opts.urdfPath   中間 URDF の保存先       (既定 slxPath と同じフォルダ)
%     opts.qFixed     固定軸の角度 [1x6 rad]   (既定 j2_config の qFixed。J2 の値は無視)
%
%   共通 I/F:
%     Inport  'tau' [スカラ] (N*m) → Outport 'q' [rad], 'dq' [rad/s]（いずれもスカラ）
%   可動関節のブロック名は 'joint1'（親資産の addC8* が joint1..nJoints を前提とするため）。
%
%   方式（Phase 0 で決定、docs/phase0_report.md）:
%     1. 親資産 writeC8URDF で 6 軸 URDF を生成
%     2. J2 以外の joint を fixed に後処理（axis/limit を削除し、固定角は origin の rpy へ畳み込む）
%     3. smimport → 重力・固定ステップ ode14x(8 kHz) を設定
%     4. J2 を 'joint1' に改名し、nJoints=1 の p で addC8Actuation/Armature/Friction を適用
%   親資産のコードは変更しない。

c = j2_setup_path();
if nargin < 1, opts = struct(); end
def = struct('modelName','j2_plant', ...
    'slxPath',fullfile(c.dataDir,'models','j2_plant.slx'), ...
    'urdfPath','', 'qFixed',c.qFixed);
f = fieldnames(def);
for i = 1:numel(f)
    if ~isfield(opts,f{i}) || isempty(opts.(f{i})), opts.(f{i}) = def.(f{i}); end
end
[slxDir,~,~] = fileparts(opts.slxPath);
if isempty(opts.urdfPath), opts.urdfPath = fullfile(slxDir,[opts.modelName '.urdf']); end
if ~isfolder(slxDir), mkdir(slxDir); end

p = c8_params(c.meshDir);
J = c.jointIdx;                                   % 可動関節 (=2)
fixedIdx = setdiff(1:p.nJoints, J);

% ---- 1-2. 6 軸 URDF を生成し、J2 以外を fixed に ----
urdf6 = [tempname '.urdf'];
cleanUrdf = onCleanup(@() delete_if_exists(urdf6));
writeC8URDF(p, urdf6);
txt = fileread(urdf6);
for i = fixedIdx
    blk = regexp(txt, sprintf('<joint name="joint%d" type="revolute">.*?</joint>', i), 'match', 'once');
    assert(~isempty(blk), 'buildJ2Plant:noJoint', 'URDF に joint%d が見つかりません', i);
    nb = strrep(blk, 'type="revolute"', 'type="fixed"');
    nb = regexprep(nb, '\s*<axis[^>]*/>', '');
    nb = regexprep(nb, '\s*<limit[^>]*/>', '');
    if opts.qFixed(i) ~= 0      % 固定角 q は子フレームの z 軸回りの回転: R = R_origin * Rz(q)
        rpy = rotm2rpy(rpy2rotm(p.joint(i).rpy) * rotz(opts.qFixed(i)));
        nb = regexprep(nb, 'rpy="[^"]*"', sprintf('rpy="%.12g %.12g %.12g"', rpy), 'once');
    end
    txt = strrep(txt, blk, nb);
end
fid = fopen(opts.urdfPath,'w');  assert(fid > 0, 'buildJ2Plant:open', '書き込めません: %s', opts.urdfPath);
fwrite(fid, txt);  fclose(fid);

% ---- 3. smimport と重力・ソルバ設定 ----
mdl = opts.modelName;
if bdIsLoaded(mdl), close_system(mdl,0); end
smimport(opts.urdfPath, 'ModelName', mdl);
mc = find_system(mdl,'MaskType','Mechanism Configuration');
set_param(mc{1}, 'GravityVector', sprintf('[%.9g %.9g %.9g]', p.gravity));
% armature(RMI) が DAE を導入するため陰的ソルバ ode14x が必要（親資産 buildC8Plant と同じ理由）
set_param(mdl,'SolverType','Fixed-step','Solver','ode14x', ...
    'FixedStep',num2str(1/c.fsSim,'%.12g'),'StopTime','0.5','SimscapeLogType','none');
try, set_param(mdl,'SimMechanicsOpenEditorOnUpdate','off'); catch, end

% ---- 4. J2 を joint1 に改名して親の add* を適用 ----
if ~isempty(find_system(mdl,'SearchDepth',1,'Name','joint1'))
    set_param([mdl '/joint1'],'Name','weldJ1');   % 固定化された J1（溶接）と名前が衝突する
end
set_param(sprintf('%s/joint%d',mdl,J),'Name','joint1');
p1 = p;  p1.nJoints = 1;
p1.armature   = p.armature(J);     p1.frictionFc  = p.frictionFc(J);
p1.frictionBv = p.frictionBv(J);   p1.frictionEps = p.frictionEps(J);
addC8Actuation(mdl, p1);
addC8Armature(mdl, p1);
addC8Friction(mdl, p1);

save_system(mdl, opts.slxPath);
fprintf('buildJ2Plant: 保存 %s (revolute=%d, weld=%d, step=%.4g s, 固定角[deg]=%s)\n', ...
    opts.slxPath, numel(find_system(mdl,'MaskType','Revolute Joint')), ...
    numel(find_system(mdl,'MaskType','Weld Joint')), 1/c.fsSim, mat2str(rad2deg(opts.qFixed),6));
end

% ================= ヘルパ =================
function delete_if_exists(f)
if exist(f,'file'), delete(f); end
end

function R = rotz(a)
R = [cos(a) -sin(a) 0; sin(a) cos(a) 0; 0 0 1];
end

function R = rpy2rotm(v)
% URDF の rpy: R = Rz(y)*Ry(p)*Rx(r)
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
