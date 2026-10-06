classdef test_j2_coverage < matlab.unittest.TestCase
%TEST_J2_COVERAGE カバレッジ評価の単体テスト（合成データ、シミュレーション不要）。

    properties
        outDir
    end

    methods (TestClassSetup)
        function setup(tc)
            repo = fileparts(fileparts(mfilename('fullpath')));
            addpath(fullfile(repo,'config'), genpath(fullfile(repo,'src')));
            j2_setup_path();
            tc.outDir = tempname;
        end
    end

    methods (TestClassTeardown)
        function cleanup(tc)
            if isfolder(tc.outDir), rmdir(tc.outDir, 's'); end
        end
    end

    methods (Test)
        function metricsAndFiles(tc)
            ds = coverageDataset(false);
            M = j2Coverage(ds, tc.outDir);
            for f = {'coverage_timeseries_A.png','coverage_timeseries_B.png','coverage_timeseries_C.png', ...
                    'coverage_phase_planes.png','coverage_phase_compare.png','coverage_frequency.png','coverage_spectrogram.png'}
                tc.verifyTrue(isfile(fullfile(tc.outDir, f{1})), f{1});
            end
            tc.verifyGreaterThan(M.occupancy3D, 0);   tc.verifyLessThanOrEqual(M.occupancy3D, 1);
            tc.verifyGreaterThanOrEqual(M.occupancy2D.theta_dtheta, M.occupancy3D);   % 射影は 3 次元占有率以上
            tc.verifyEqual(M.bins, 20);
            % ベンチマークを学習系と同じ点にすると 100% 入る
            tc.verifyEqual(M.benchmarkInOccupied3D, 1, 'AbsTol', 1e-12);
            tc.verifyEqual(M.nTransition.learning, M.nTransition.benchmark * 2);
        end

        function benchmarkOutsideLearningRegionIsDetected(tc)
            ds = coverageDataset(true);                 % ベンチマークを学習系とは別の領域に置く
            M = j2Coverage(ds, tc.outDir);
            tc.verifyLessThan(M.benchmarkInOccupied3D, 0.2);
        end

        function contactMetricsAndFigure(tc)
            jp = j2Params();  ds = coverageDataset(false);
            K = 3000;  tol = deg2rad(0.1);
            q = linspace(jp.qMin + deg2rad(5), jp.qMax - deg2rad(5), K+1)';
            q(1:1000) = jp.qMin - deg2rad(0.2);  q(end-499:end) = jp.qMax + deg2rad(0.2);
            dq = 0.5*sin(0.3*(1:K+1)');  dq(1:200) = 0;                  % 先頭 200 点は静止して押し付け
            tau = zeros(K,1);
            lim = @(x) x <= jp.qMin + tol | x >= jp.qMax - tol;
            c = struct('name','contact_fall_001','pattern','contact_fall','type','excitation','phase',3, ...
                'split','train','seed',9,'q',q,'dq',dq,'tau',tau,'tauInst',tau,'qref',[], ...
                'atLimit',lim(q(1:end-1)) | lim(q(2:end)),'nTransition',K);
            ds.scen(end+1) = c;
            M = j2Coverage(ds, tc.outDir);
            tc.verifyEqual(M.contact.nScenarios, 1);
            tc.verifyTrue(isfile(fullfile(tc.outDir,'coverage_contact.png')));
            tc.verifyEqual(M.contact.flaggedTransitions, nnz(c.atLimit));
            % 動的 = 動いている接触。負側は先頭 1000 遷移のうち静止 200 を除く、正側は末尾 500
            tc.verifyEqual(M.contact.dynamicTransitions.neg, 800, 'AbsTol', 2);
            tc.verifyEqual(M.contact.dynamicTransitions.pos, 500, 'AbsTol', 2);
            tc.verifyGreaterThan(M.contact.staticPressingFraction, 0.1);
            tc.verifyEqual(M.contact.maxPenetrationDeg.neg, 0.2, 'AbsTol', 1e-6);
        end

        function frequencyMetricsFindTheBand(tc)
            % 20 Hz 帯域制限ノイズ: −40 dB 以内の被覆は 20 Hz 付近まで、ピークは 20 Hz 未満
            ds = coverageDataset(false);
            M = j2Coverage(ds, tc.outDir);
            tc.verifyLessThan(M.frequency.tau.fPeak_Hz, 20);
            tc.verifyGreaterThan(M.frequency.tau.fCover40dB_Hz, 15);
            tc.verifyLessThan(M.frequency.tau.fCover40dB_Hz, 120);
        end
    end
end

function ds = coverageDataset(shiftBenchmark)
% 学習系 2 本 + ベンチマーク 1 本（同じ軌跡を 2 本・1 本）。1 kHz、3000 遷移ずつ
jp = j2Params();  K = 3000;  fs = jp.fsData;
rs = RandStream('twister','Seed',3);
[b, a] = butter(4, 20/(fs/2));
mk = @(name, pat, split, phase, off) makeScen(name, pat, split, phase, off, rs, b, a, K, fs, jp);
S(1) = mk('bln_1', 'bln_normal', 'train', 1, 0);
S(2) = mk('bln_2', 'bln_normal', 'val',   1, 0);
S(3) = mk('ptp_1', 'ptp_v50',    'benchmark', 'benchmark', 0);
if ~shiftBenchmark
    S(3).q = S(1).q;  S(3).dq = S(1).dq;  S(3).tau = S(1).tau;       % 学習系 1 本と同一
    S(2).q = S(1).q;  S(2).dq = S(1).dq;  S(2).tau = S(1).tau;       % 学習系 2 本が同一（遷移数 2 倍）
else
    S(3).q = jp.qMax - deg2rad(2) + 0*S(3).q;  S(3).dq = 4.5*ones(K+1,1);  S(3).tau = 800*ones(K,1);
end
ds.scen = S(:);
ds.stats = struct();  ds.scale = struct();  ds.manifest = struct();
sc.DTHETA_MAX = jp.qdMax;  sc.TAU_MAX = jp.tauPeak;  sc.D_THETA_MAX = jp.qdMax/fs;  sc.D_DTHETA_MAX = (jp.tauPeak/jp.M)/fs;
ds.scale = sc;
end

function s = makeScen(name, pat, split, phase, ~, rs, b, a, K, fs, jp)
u = filtfilt(b, a, randn(rs, K, 1));  u = 200*u/max(abs(u));
w = filtfilt(b, a, randn(rs, K+1, 1));  w = 2*w/max(abs(w));
th = deg2rad(-100) + cumsum(w)/fs;  th = min(max(th, jp.qMin + 0.05), jp.qMax - 0.05);
s = struct('name',name,'pattern',pat,'type','excitation','phase',phase,'split',split,'seed',1, ...
    'q',th,'dq',w,'tau',u,'tauInst',u,'qref',[],'atLimit',false(K,1),'nTransition',K);
end
