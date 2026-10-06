function R = validateJ2Data(ds, verbose)
%VALIDATEJ2DATA 生成したデータセットの妥当性を検証する（要件 FR-6、ゲート G3）。
%   R = VALIDATEJ2DATA(ds, verbose)
%   R.checks: 検査名・合否・値の構造体配列、R.nErrors: 不合格数、R.metrics: 参考指標。
%   検査:
%     有限性 / トルク飽和内 / 速度制限内 / 可動範囲（±3° のめり込みまで）/ 配列長の整合 /
%     split の整合（重複・ベンチマークの分離）/ 正規化統計が train の全遷移由来 /
%     リミット拘束の割合 / 解析モデルとの 1 ステップ整合（RMS 残差が加速度 RMS の半分未満）/
%     接触シナリオ（contact_*）が両側のリミットに実際に接触し、接触遷移が 1000 件以上
if nargin < 2, verbose = true; end
jp = j2Params();
S = ds.scen;
R.checks = struct('name',{},'pass',{},'value',{});
add = @(name, pass, value) addCheck(name, pass, value, verbose);
chk = {};

allQ = vertcat(S.q);  allDq = vertcat(S.dq);  allTau = vertcat(S.tau);
chk{end+1} = add('すべて有限', all(isfinite([allQ; allDq; allTau])), nnz(~isfinite([allQ; allDq; allTau])));
% キャッシュは single 保存のため相対 1e-6 の丸めを許容
chk{end+1} = add('|tau| <= tauPeak', max(abs(allTau)) <= jp.tauPeak*(1 + 1e-6), max(abs(allTau)));
chk{end+1} = add('|dq| < qdMax', max(abs(allDq)) < jp.qdMax, max(abs(allDq)));
chk{end+1} = add('q が可動範囲 ±3° 以内', min(allQ) > jp.qMin - deg2rad(3) && max(allQ) < jp.qMax + deg2rad(3), ...
    rad2deg([min(allQ) max(allQ)]));
okLen = all(arrayfun(@(s) numel(s.q) == numel(s.tau)+1 && numel(s.dq) == numel(s.tau)+1 && numel(s.atLimit) == numel(s.tau), S));
chk{end+1} = add('配列長の整合 (q,dq = K+1, tau,atLimit = K)', okLen, nnz(~okLen));

names = {S.name};
chk{end+1} = add('シナリオ名が一意', numel(unique(names)) == numel(names), numel(names) - numel(unique(names)));
isBench = strcmp({S.phase},'benchmark');
chk{end+1} = add('benchmark は split=benchmark のみ', all(strcmp({S(isBench).split},'benchmark')) && ...
    all(ismember({S(~isBench).split},{'train','val','test'})), nnz(isBench));
% 小さい preset（本数の少ないパターンは train に寄せる）では val/test が無くてよい。100 本以上で必須
chk{end+1} = add('train/val/test がそろっている（100 本以上のとき）', numel(S) < 100 || all(ismember({'train','val','test'}, {S.split})), 0);

% 正規化統計が train の全遷移（拘束を含む）由来であること（組み立て時と同じ計算で再現）
X = [];  Y = [];
for k = find(strcmp({S.split},'train'))
    sc = S(k);
    X = [X; sc.q(1:end-1), sc.dq(1:end-1), sc.tau];  Y = [Y; diff(sc.q), diff(sc.dq)]; %#ok<AGROW>
end
chk{end+1} = add('正規化統計 = train 全遷移の統計', ...
    max(abs(ds.stats.xMean - mean(X,1))./max(abs(mean(X,1)),1e-9)) < 1e-9 && ...
    max(abs(ds.stats.yStd - std(Y,0,1))./std(Y,0,1)) < 1e-9, 0);

nAll = sum([S.nTransition]);  nLim = sum(arrayfun(@(s) nnz(s.atLimit), S));
R.metrics.atLimitFraction = nLim / nAll;
chk{end+1} = add('リミット拘束の割合 < 30%', R.metrics.atLimitFraction < 0.30, R.metrics.atLimitFraction);

% 解析モデルとの 1 ステップ整合: 区間平均トルクでの加速度の予測と実測の差
dt = 1/jp.fsData;  num = 0;  den = 0;
for k = 1:numel(S)
    sc = S(k);  ok = ~sc.atLimit;
    qm  = (sc.q(1:end-1) + sc.q(2:end))/2;  dqm = (sc.dq(1:end-1) + sc.dq(2:end))/2;
    aMeas  = diff(sc.dq)/dt;
    aModel = (sc.tau + jp.mgL*sin(qm) - jp.Bv*dqm - jp.Fc*tanh(dqm/jp.eps)) / jp.M;
    num = num + sum((aMeas(ok) - aModel(ok)).^2);  den = den + sum(aMeas(ok).^2);
end
R.metrics.modelResidualRatio = sqrt(num/den);
chk{end+1} = add('解析モデルの 1 ステップ加速度残差 / 加速度 RMS < 0.5', R.metrics.modelResidualRatio < 0.5, R.metrics.modelResidualRatio);

% リミット接触シナリオ（pattern が contact_ で始まる）が、群として両側のリミットに実際に接触していること
isC = startsWith({S.pattern}, 'contact_');
if any(isC)
    qc = vertcat(S(isC).q);  tol = deg2rad(0.1);
    chk{end+1} = add('接触シナリオが −側リミットに接触', min(qc) <= jp.qMin + tol, rad2deg(min(qc)));
    chk{end+1} = add('接触シナリオが +側リミットに接触', max(qc) >= jp.qMax - tol, rad2deg(max(qc)));
    nC = sum(arrayfun(@(s) nnz(s.atLimit), S(isC)));
    R.metrics.contactTransitions = nC;
    chk{end+1} = add('接触の遷移が 1000 件以上', nC >= 1000, nC);
    % 接触遷移の大半は静止したまま押し付けている（情報が重複）。動いている接触（衝突・反発・離脱）の件数も側ごとに確認
    nDynNeg = 0;  nDynPos = 0;
    for k = find(isC)
        sc = S(k);
        dyn = sc.atLimit & (abs(sc.dq(1:end-1)) > 0.02 | abs(diff(sc.dq)) > 0.005);
        low = sc.q(1:end-1) <= jp.qMin + tol;
        nDynNeg = nDynNeg + nnz(dyn & low);  nDynPos = nDynPos + nnz(dyn & ~low);
    end
    R.metrics.contactDynamicTransitions = struct('neg', nDynNeg, 'pos', nDynPos);
    chk{end+1} = add('動的な接触の遷移が各側 1000 件以上', nDynNeg >= 1000 && nDynPos >= 1000, [nDynNeg nDynPos]);
end

R.checks = [chk{:}];
R.nErrors = nnz(~[R.checks.pass]);
R.metrics.nTransition = nAll;
if verbose, fprintf('検証: %d 項目中 %d 件不合格（遷移 %d、リミット拘束 %.2f%%、モデル残差比 %.3f）\n', ...
    numel(R.checks), R.nErrors, nAll, 100*R.metrics.atLimitFraction, R.metrics.modelResidualRatio); end
end

function c = addCheck(name, pass, value, verbose)
c = struct('name',name,'pass',logical(pass),'value',value);
if verbose, fprintf('  [%s] %s  (%s)\n', ternary(pass,'OK','NG'), name, mat2str(value,4)); end
end

function s = ternary(c,a,b)
if c, s = a; else, s = b; end
end
