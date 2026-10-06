classdef test_j2_dataset < matlab.unittest.TestCase
%TEST_J2_DATASET 8k→1k 変換・フラット化・検証関数の単体テスト（合成データ、シミュレーション不要）。

    methods (TestClassSetup)
        function setup(tc) %#ok<MANU>
            repo = fileparts(fileparts(mfilename('fullpath')));
            addpath(fullfile(repo,'config'), genpath(fullfile(repo,'src')));
            j2_setup_path();
        end
    end

    methods (Test)
        function downsampleKeepsStatesAndAveragesTorque(tc)
            fs = 8000;  N = 8*50 + 1;  t = (0:N-1)'/fs;
            raw.q = sin(2*pi*3*t);  raw.dq = cos(2*pi*3*t);
            raw.tau = 100*double(mod(floor(t*1000), 2) == 0);      % 1 kHz の矩形波（半周期 1 ms）
            D = j2Downsample(raw, 8000, 1000);
            tc.verifyEqual(numel(D.q), 51);  tc.verifyEqual(numel(D.tau), 50);
            tc.verifyEqual(D.q, raw.q(1:8:end), 'AbsTol', 0);        % 状態は瞬時値のまま（フィルタなし）
            % 区間 [t_k, t_k+1ms) はちょうど 1 周期の半分: 偶数区間は 100、奇数区間は 0
            tc.verifyEqual(D.tau(1:2:end), 100*ones(25,1), 'AbsTol', 1e-12);
            tc.verifyEqual(D.tau(2:2:end), zeros(25,1), 'AbsTol', 1e-12);
            tc.verifyEqual(D.tauInst, raw.tau(1:8:end-8), 'AbsTol', 0);
        end

        function downsampleAveragesStepWithinInterval(tc)
            % 区間の途中でステップした場合は、平均が部分的な値になる
            raw.q = zeros(17,1);  raw.dq = zeros(17,1);
            raw.tau = [zeros(12,1); 80*ones(5,1)];                  % 13 点目以降 80
            D = j2Downsample(raw, 8000, 1000);
            % 区間 2 は 9〜16 点目の 8 点。うち 13〜16 点目の 4 点が 80 なので平均は 40
            tc.verifyEqual(D.tau, [0; 40], 'AbsTol', 1e-12);
        end

        function flatAndSpecFeatures(tc)
            ds = syntheticDataset();
            [X, Y, info] = j2BuildFlat(ds, 'train');
            tc.verifyEqual(size(X,2), 3);  tc.verifyEqual(size(Y,2), 2);
            tc.verifyEqual(size(X,1), numel(info.k));
            % 既定ではリミット拘束の遷移も含める。excludeAtLimit=true で除外（シナリオ 1 の最後 2 遷移）
            [X2, ~] = j2BuildFlat(ds, 'train', struct('excludeAtLimit', true));
            tc.verifyEqual(size(X,1) - size(X2,1), 2);
            [Xs, Ys] = j2BuildFlat(ds, 'train', struct('features','spec'));
            tc.verifyEqual(size(Xs,2), 4);
            tc.verifyEqual(Xs(:,1).^2 + Xs(:,2).^2, ones(size(Xs,1),1), 'AbsTol', 1e-12);
            tc.verifyLessThanOrEqual(max(abs(Xs(:,4))), 1 + 1e-12);   % トルクは tauPeak で正規化
            tc.verifyEqual(size(Ys,2), 2);
            % 複数 split
            [Xa, ~] = j2BuildFlat(ds, {'train','val'});
            tc.verifyGreaterThan(size(Xa,1), size(X,1));
        end

        function validatorChecksContactScenarios(tc)
            % 接触シナリオ（contact_*）が両側に接触し、動的な接触が各側 1000 件以上であることを検査
            jp = j2Params();  ds = syntheticDataset();
            names0 = {ds.scen.name};  %#ok<NASGU>
            ds.scen(end+1) = contactScenario(jp, 'contact_drive', 'train', 3000);       % 動的な接触が十分ある
            R = validateJ2Data(ds, false);
            nm = {R.checks.name};
            tc.verifyTrue(R.checks(contains(nm, '−側リミットに接触')).pass);
            tc.verifyTrue(R.checks(contains(nm, '+側リミットに接触')).pass);
            tc.verifyTrue(R.checks(contains(nm, '動的な接触')).pass);
            % 静止して押し付けているだけの接触は「動的」に数えない
            ds2 = syntheticDataset();
            ds2.scen(end+1) = contactScenario(jp, 'contact_fall', 'train', 3000);
            ds2.scen(end).dq(:) = 0;                                                    % 全区間で静止
            R2 = validateJ2Data(ds2, false);
            tc.verifyFalse(R2.checks(contains({R2.checks.name}, '動的な接触')).pass);
        end

        function validatorFlagsBadData(tc)
            ds = syntheticDataset();
            R = validateJ2Data(ds, false);
            tc.verifyEqual(R.nErrors, 1);                           % 合成データは解析モデルと整合しない（残差比のみ NG）
            tc.verifyFalse(R.checks(end).pass);
            % トルク飽和超過・NaN を入れると検出される
            ds.scen(1).tau(3) = 1e4;
            ds.scen(2).q(4) = NaN;
            R2 = validateJ2Data(ds, false);
            names = {R2.checks(~[R2.checks.pass]).name};
            tc.verifyTrue(any(contains(names, 'tauPeak')));
            tc.verifyTrue(any(contains(names, '有限')));
        end
    end
end

function s = contactScenario(jp, pattern, split, K)
% 前半は −側、後半は +側のリミットに接触する合成シナリオ（接触区間では dq が振動する）
tol = deg2rad(0.1);
q = linspace(jp.qMin + deg2rad(5), jp.qMax - deg2rad(5), K+1)';
q(1:1200) = jp.qMin - deg2rad(0.2);  q(end-1199:end) = jp.qMax + deg2rad(0.2);   % 両側で 1200 点ずつ接触
dq = 0.5*sin(0.3*(1:K+1)');  tau = zeros(K,1);
atLimit = (q(1:end-1) <= jp.qMin + tol) | (q(2:end) <= jp.qMin + tol) | (q(1:end-1) >= jp.qMax - tol) | (q(2:end) >= jp.qMax - tol);
s = struct('name',[pattern '_001'],'pattern',pattern,'type','excitation','phase',3,'split',split,'seed',1, ...
    'q',q,'dq',dq,'tau',tau,'tauInst',tau,'qref',[],'atLimit',atLimit,'nTransition',K);
end

function ds = syntheticDataset()
% 4 本（train 2 / val 1 / test 1）の小さな合成データ（解析モデルとは無関係）
jp = j2Params();
rs = RandStream('twister','Seed',1);
splits = {'train','train','val','test'};
for k = 1:4
    K = 20;
    sc.name = sprintf('syn_%d',k);  sc.pattern = 'syn';  sc.type = 'excitation';  sc.phase = 1;
    sc.split = splits{k};  sc.seed = k;
    sc.q  = deg2rad(-100) + 0.01*cumsum(randn(rs,K+1,1));
    sc.dq = 0.1*randn(rs,K+1,1);
    sc.tau = 50*randn(rs,K,1);  sc.tauInst = sc.tau;  sc.qref = [];
    sc.atLimit = false(K,1);  sc.nTransition = K;
    if k == 1, sc.atLimit(end-1:end) = true; end
    S(k) = sc; %#ok<AGROW>
end
isTr = strcmp({S.split},'train');  X = [];  Y = [];
for k = find(isTr)
    X = [X; S(k).q(1:end-1), S(k).dq(1:end-1), S(k).tau]; Y = [Y; diff(S(k).q), diff(S(k).dq)]; %#ok<AGROW>
end
st.xMean = mean(X,1); st.xStd = std(X,0,1); st.yMean = mean(Y,1); st.yStd = std(Y,0,1);
dt = 1/jp.fsData;
sc2.DTHETA_MAX = jp.qdMax; sc2.TAU_MAX = jp.tauPeak; sc2.D_THETA_MAX = jp.qdMax*dt; sc2.D_DTHETA_MAX = (jp.tauPeak/jp.M)*dt;
ds.scen = S(:);  ds.stats = st;  ds.scale = sc2;  ds.manifest = struct();
end
