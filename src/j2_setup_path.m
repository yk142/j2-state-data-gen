function c = j2_setup_path()
%J2_SETUP_PATH 親資産と本リポジトリを path に追加し、設定を返す。
%   c = J2_SETUP_PATH()   親資産の config/src と本リポジトリの config/src を addpath する。

c = j2_config();
assert(isfolder(c.parentDir), 'j2:noParent', '親資産が見つかりません: %s', c.parentDir);
addpath(fullfile(c.parentDir,'config'));
addpath(genpath(fullfile(c.parentDir,'src')));
addpath(fullfile(c.rootDir,'config'), fullfile(c.rootDir,'src'));
end
