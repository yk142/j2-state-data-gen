function [X, Y, info] = j2BuildFlat(ds, split, opts)
%J2BUILDFLAT データセットから (入力, 目標) のフラットな行列を作る。
%   [X, Y, info] = J2BUILDFLAT(ds, split, opts)
%     split  'train' | 'val' | 'test' | 'benchmark'（文字列の cell で複数指定可）
%     opts.excludeAtLimit  リミット拘束の遷移を除く（既定 true）
%     opts.features        'raw'（既定）: X=[q dq tau]、Y=[Δq Δdq]（SI 単位）
%                          'spec'      : 引継ぎ資料の形（正規化済み）
%                              X=[sin q, cos q, dq/DTHETA_MAX, tau/TAU_MAX]
%                              Y=[Δq/D_THETA_MAX, Δdq/D_DTHETA_MAX]
%   info.scenario  各行が属するシナリオの添字、info.k 遷移の添字
if nargin < 3, opts = struct(); end
if ~isfield(opts,'excludeAtLimit'), opts.excludeAtLimit = true; end
if ~isfield(opts,'features'), opts.features = 'raw'; end
if ischar(split), split = {split}; end
X = [];  Y = [];  sIdx = [];  kIdx = [];
for k = find(ismember({ds.scen.split}, split))
    sc = ds.scen(k);
    Xk = [sc.q(1:end-1), sc.dq(1:end-1), sc.tau];
    Yk = [diff(sc.q), diff(sc.dq)];
    ok = true(size(sc.tau));
    if opts.excludeAtLimit, ok = ~sc.atLimit; end
    X = [X; Xk(ok,:)];  Y = [Y; Yk(ok,:)];                       %#ok<AGROW>
    sIdx = [sIdx; repmat(k, nnz(ok), 1)];  kIdx = [kIdx; find(ok)]; %#ok<AGROW>
end
if strcmp(opts.features, 'spec')
    sc = ds.scale;
    X = [sin(X(:,1)), cos(X(:,1)), X(:,2)/sc.DTHETA_MAX, X(:,3)/sc.TAU_MAX];
    Y = [Y(:,1)/sc.D_THETA_MAX, Y(:,2)/sc.D_DTHETA_MAX];
end
info.scenario = sIdx;  info.k = kIdx;
end
