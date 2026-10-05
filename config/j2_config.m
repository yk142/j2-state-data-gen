function c = j2_config()
%J2_CONFIG J2 単軸データ生成の設定（親資産パス、固定姿勢、サンプリング）を返す。
%   c = J2_CONFIG()
%
%   親資産 robotArmSurrogate-matlab は参照専用（コピー・改変しない）。
%   パスは環境変数 J2_PARENT_DIR で上書きできる。addpath はここだけで行う。

c.parentDir = getenv('J2_PARENT_DIR');
if isempty(c.parentDir)
    c.parentDir = fullfile(getenv('HOME'),'matlab-projects','robotArmSurrogate-matlab');
end
c.meshDir = fullfile(c.parentDir,'urdf','mesh');

c.jointIdx = 2;                  % 可動にする関節（J2）
% 固定軸の角度 [rad]（J2 の値は無視される）。
% J3 = 75.0684 deg は、J2 の重力トルクが θ=0 で 0 になる値（c8RNEA を fzero で求解）。
% J3=0 のままだと J2 の重力平衡点が +11 deg にずれるため、θ=0 を倒立平衡点にする目的で選んだ。
c.qFixed   = [0 0 deg2rad(75.0684) 0 0 0];
c.fsSim    = 8000;               % シミュレーション周波数 [Hz]
c.fsData   = 1000;               % 学習データ周波数 [Hz]
c.seedBase = 20261006;           % 乱数シードの基準

root = fileparts(fileparts(mfilename('fullpath')));
c.rootDir = root;
c.dataDir = fullfile(root,'data');
end
