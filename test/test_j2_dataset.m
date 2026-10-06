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
            % リミット拘束の遷移は除外される（シナリオ 1 の最後 2 遷移）
            [X2, ~] = j2BuildFlat(ds, 'train', struct('excludeAtLimit', false));
            tc.verifyEqual(size(X2,1) - size(X,1), 2);
            [Xs, Ys] = j2BuildFlat(ds, 'train', struct('features','spec'));
            tc.verifyEqual(size(Xs,2), 4);
            tc.verifyEqual(Xs(:,1).^2 + Xs(:,2).^2, ones(size(Xs,1),1), 'AbsTol', 1e-12);
            tc.verifyLessThanOrEqual(max(abs(Xs(:,4))), 1 + 1e-12);   % トルクは tauPeak で正規化
            tc.verifyEqual(size(Ys,2), 2);
            % 複数 split
            [Xa, ~] = j2BuildFlat(ds, {'train','val'});
            tc.verifyGreaterThan(size(Xa,1), size(X,1));
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
    ok = ~S(k).atLimit;  Xk = [S(k).q(1:end-1), S(k).dq(1:end-1), S(k).tau];  Yk = [diff(S(k).q), diff(S(k).dq)];
    X = [X; Xk(ok,:)]; Y = [Y; Yk(ok,:)]; %#ok<AGROW>
end
st.xMean = mean(X,1); st.xStd = std(X,0,1); st.yMean = mean(Y,1); st.yStd = std(Y,0,1);
dt = 1/jp.fsData;
sc2.DTHETA_MAX = jp.qdMax; sc2.TAU_MAX = jp.tauPeak; sc2.D_THETA_MAX = jp.qdMax*dt; sc2.D_DTHETA_MAX = (jp.tauPeak/jp.M)*dt;
ds.scen = S(:);  ds.stats = st;  ds.scale = sc2;  ds.manifest = struct();
end
