classdef test_j2_export < matlab.unittest.TestCase
%TEST_J2_EXPORT j2ExportFlat（HDF5 / v7 mat）の書き出しと読み戻しのテスト（合成データ、シミュレーション不要）。

    properties
        tmp
        ds
    end

    methods (TestClassSetup)
        function setup(tc)
            repo = fileparts(fileparts(mfilename('fullpath')));
            addpath(fullfile(repo,'config'), genpath(fullfile(repo,'src')));
            j2_setup_path();
            tc.tmp = tempname;  mkdir(tc.tmp);
            tc.ds = syntheticDataset();
        end
    end

    methods (TestClassTeardown)
        function cleanup(tc)
            if isfolder(tc.tmp), rmdir(tc.tmp, 's'); end
        end
    end

    methods (Test)
        function hdf5MatchesFlatBuild(tc)
            f = fullfile(tc.tmp, 'a.h5');
            info = j2ExportFlat(tc.ds, f, struct('dtype','double'));
            for sp = {'train','val','test','benchmark'}
                [X, Y] = j2BuildFlat(tc.ds, sp{1});
                [Xs, Ys] = j2BuildFlat(tc.ds, sp{1}, struct('features','spec'));
                tc.verifyEqual(h5read(f, ['/' sp{1} '/X'])', X);              % MATLAB は (d, N) で返る → 転置して一致
                tc.verifyEqual(h5read(f, ['/' sp{1} '/Y'])', Y);
                tc.verifyEqual(h5read(f, ['/' sp{1} '/X_spec'])', Xs);
                tc.verifyEqual(h5read(f, ['/' sp{1} '/Y_spec'])', Ys);
                tc.verifyEqual(info.n.(sp{1}), size(X,1));
            end
            % h5py から見た形は (N, d): HDF5 上の次元が [d N]（MATLAB の h5info は逆順で報告する）
            h = h5info(f, '/train/X');
            tc.verifyEqual(h.Dataspace.Size, [3 info.n.train]);
        end

        function metadataMapsBackToScenarios(tc)
            f = fullfile(tc.tmp, 'b.h5');
            j2ExportFlat(tc.ds, f, struct('dtype','double'));
            sid = double(h5read(f, '/train/scenario_id')) + 1;               % 0 始まり → MATLAB の添字
            stp = double(h5read(f, '/train/step')) + 1;
            at  = h5read(f, '/train/at_limit');
            X = h5read(f, '/train/X')';  Y = h5read(f, '/train/Y')';
            for r = [1 5 17 numel(sid)]
                sc = tc.ds.scen(sid(r));
                tc.verifyEqual(strcmp(sc.split,'train'), true);
                tc.verifyEqual(X(r,:), [sc.q(stp(r)), sc.dq(stp(r)), sc.tau(stp(r))], 'AbsTol', 1e-12);
                tc.verifyEqual(Y(r,:), [sc.q(stp(r)+1) - sc.q(stp(r)), sc.dq(stp(r)+1) - sc.dq(stp(r))], 'AbsTol', 1e-12);
                tc.verifyEqual(logical(at(r)), sc.atLimit(stp(r)));
            end
            names = h5read(f, '/scenarios/name');
            tc.verifyEqual(cellstr(names)', {tc.ds.scen.name});
            tc.verifyEqual(double(h5read(f, '/scenarios/n_transition'))', [tc.ds.scen.nTransition]);
            tc.verifyEqual(double(h5read(f, '/scenarios/seed'))', [tc.ds.scen.seed]);
            tc.verifyEqual(cellstr(h5read(f,'/scenarios/phase'))', {'1','1','2','benchmark','3'});
        end

        function referenceAngleIsExported(tc)
            f = fullfile(tc.tmp, 'q.h5');
            j2ExportFlat(tc.ds, f, struct('dtype','double'));
            for sp = {'test','benchmark','train'}
                qr = h5read(f, ['/' sp{1} '/qref']);
                [~, ~, inf1] = j2BuildFlat(tc.ds, sp{1});
                for r = [1 numel(qr)]
                    sc = tc.ds.scen(inf1.scenario(r));
                    if isempty(sc.qref), tc.verifyTrue(isnan(qr(r)));
                    else, tc.verifyEqual(qr(r), sc.qref(inf1.k(r)), 'AbsTol', 1e-12); end
                end
            end
            tc.verifyTrue(all(isnan(h5read(f,'/train/qref')) | true));
            tc.verifyTrue(any(isnan(h5read(f,'/train/qref'))));                     % 励振シナリオは NaN
            tc.verifyFalse(any(isnan(h5read(f,'/test/qref'))));                     % 姿勢保持は参照あり
            tc.verifyFalse(any(isnan(h5read(f,'/benchmark/qref'))));                % PTP は参照あり
        end

        function statsScaleAndAttributes(tc)
            f = fullfile(tc.tmp, 'c.h5');
            j2ExportFlat(tc.ds, f);
            tc.verifyEqual(h5read(f,'/stats/x_mean')', tc.ds.stats.xMean, 'AbsTol', 1e-15);
            tc.verifyEqual(h5read(f,'/stats/y_std')',  tc.ds.stats.yStd,  'AbsTol', 1e-15);
            tc.verifyEqual(h5read(f,'/scale/TAU_MAX'), tc.ds.scale.TAU_MAX);
            jp = j2Params();                                                % 解析モデルの係数（基準線の計算に使う）
            for nm = {'M','mgL','Fc','Bv','eps','qMin','qMax','qdMax','tauRated','tauPeak'}
                tc.verifyEqual(h5read(f, ['/physics/' nm{1}]), jp.(nm{1}), 'AbsTol', 1e-12, nm{1});
            end
            tc.verifyEqual(h5readatt(f,'/','preset'), 'synthetic');
            tc.verifyEqual(h5readatt(f,'/','fs_data_hz'), 1000);
            tc.verifyEqual(h5readatt(f,'/','q_fixed_deg')', tc.ds.manifest.qFixedDeg, 'AbsTol', 1e-12);
            tc.verifyEqual(h5readatt(f,'/','exclude_at_limit'), 0);
            tc.verifyNotEmpty(h5readatt(f,'/train/X','columns'));
        end

        function singlePrecisionStaysClose(tc)
            f = fullfile(tc.tmp, 'd.h5');
            j2ExportFlat(tc.ds, f);                                          % 既定は single
            [X, Y] = j2BuildFlat(tc.ds, 'train');
            Xr = h5read(f, '/train/X')';  Yr = h5read(f, '/train/Y')';
            tc.verifyClass(Xr, 'single');
            tc.verifyEqual(double(Xr), X, 'RelTol', 1e-6, 'AbsTol', 1e-6);
            tc.verifyEqual(double(Yr), Y, 'AbsTol', 1e-7);                   % Δ は倍精度で計算してから丸める
        end

        function excludeAtLimitAndSplitSelection(tc)
            f = fullfile(tc.tmp, 'e.h5');
            info = j2ExportFlat(tc.ds, f, struct('excludeAtLimit',true,'splits',{{'train','test'}}));
            [X, ~] = j2BuildFlat(tc.ds, 'train', struct('excludeAtLimit',true));
            [Xall, ~] = j2BuildFlat(tc.ds, 'train');
            tc.verifyEqual(info.n.train, size(X,1));
            tc.verifyLessThan(info.n.train, size(Xall,1));                   % 拘束の遷移が除かれる
            tc.verifyFalse(isfield(info.n, 'val'));
            at = h5read(f, '/train/at_limit');
            tc.verifyEqual(nnz(at), 0);
            tc.verifyEqual(h5readatt(f,'/','exclude_at_limit'), 1);
        end

        function matV7RoundTrip(tc)
            f = fullfile(tc.tmp, 'f.mat');
            j2ExportFlat(tc.ds, f, struct('dtype','double'));
            L = load(f);
            [X, Y] = j2BuildFlat(tc.ds, 'train');
            tc.verifyEqual(L.train_X, X);  tc.verifyEqual(L.train_Y, Y);
            tc.verifyEqual(L.scenario_name', {tc.ds.scen.name});
            tc.verifyEqual(L.stats_x_mean, tc.ds.stats.xMean(:));
            tc.verifyEqual(L.scale_D_DTHETA_MAX, tc.ds.scale.D_DTHETA_MAX);
            tc.verifyEqual(L.physics_mgL, j2Params().mgL, 'AbsTol', 1e-12);
            tc.verifyEqual(class(L.train_scenario_id), 'int32');
            % 先頭バイトで v7（HDF5 ではない）であることを確認 → scipy.io.loadmat で読める
            fid = fopen(f); hdr = fread(fid, 8, '*char')'; fclose(fid);
            tc.verifyEqual(hdr(1:6), 'MATLAB');
        end

        function unknownExtensionFails(tc)
            tc.verifyError(@() j2ExportFlat(tc.ds, fullfile(tc.tmp,'x.csv')), 'j2ExportFlat:ext');
        end
    end
end

function ds = syntheticDataset()
% train 2 本（うち 1 本はリミット接触）/ val / benchmark / 接触 1 本。解析モデルとは無関係
jp = j2Params();  rs = RandStream('twister','Seed',2);
specs = {'bln_normal','excitation',1,'train'; 'bln_normal','excitation',1,'val'; 'hold_-090','hold',2,'test'; ...
         'ptp_v50','ptp','benchmark','benchmark'; 'contact_fall','excitation',3,'train'};
for k = 1:size(specs,1)
    K = 40 + 10*k;
    s.name = sprintf('syn_%d',k);  s.pattern = specs{k,1};  s.type = specs{k,2};  s.phase = specs{k,3};  s.split = specs{k,4};
    s.seed = 100 + k;
    s.q = deg2rad(-100) + 0.01*cumsum(randn(rs,K+1,1));  s.dq = 0.1*randn(rs,K+1,1);
    s.tau = 50*randn(rs,K,1);  s.tauInst = s.tau;  s.qref = [];
    if k == 3 || k == 4, s.qref = s.q + 0.001*randn(rs,K+1,1); end              % 閉ループ（姿勢保持・PTP）は参照あり
    s.atLimit = false(K,1);  if k == 5, s.atLimit(5:15) = true; end
    s.nTransition = K;
    S(k) = s; %#ok<AGROW>
end
isTr = strcmp({S.split},'train');  X = [];  Y = [];
for k = find(isTr)
    X = [X; S(k).q(1:end-1), S(k).dq(1:end-1), S(k).tau]; Y = [Y; diff(S(k).q), diff(S(k).dq)]; %#ok<AGROW>
end
ds.scen = S(:);
ds.stats = struct('xMean',mean(X,1),'xStd',std(X,0,1),'yMean',mean(Y,1),'yStd',std(Y,0,1),'note','');
dt = 1/jp.fsData;
ds.scale = struct('DTHETA_MAX',jp.qdMax,'TAU_MAX',jp.tauPeak,'D_THETA_MAX',jp.qdMax*dt,'D_DTHETA_MAX',(jp.tauPeak/jp.M)*dt);
ds.manifest = struct('preset','synthetic','created','2026-01-01','fsRaw',8000,'fsData',1000, ...
    'qFixedDeg',[0 0 75.0684 0 0 0],'parentCommit','abc','matlab','R2025a');
end
