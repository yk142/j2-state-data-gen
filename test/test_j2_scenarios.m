classdef test_j2_scenarios < matlab.unittest.TestCase
%TEST_J2_SCENARIOS 各シナリオ種別を実際に 8 kHz で実行して検証する（ゲート G2）。
%   励振・閉ループモデルを一時フォルダに生成し、small preset の代表シナリオを実行する。
%   実行時間は約 5 分（モデル生成 + シナリオごとの Simscape コンパイル）。

    properties
        jp
        S
        models
        tmpDir
    end

    properties (TestParameter)
        scenarioName = {'bln_normal_001','bln_high_001','bln_low_001','freefall_001', ...
            'chirp_001','step_001','micro_001','hold_-090_001','hold_+000_001', ...
            'nearlimit_001','ptp_v70_001', ...
            'contact_fall_001','contact_fall_002','contact_drive_001','contact_bln_001'};
    end

    methods (TestClassSetup)
        function build(tc)
            repo = fileparts(fileparts(mfilename('fullpath')));
            addpath(fullfile(repo,'config'), genpath(fullfile(repo,'src')));
            j2_setup_path();
            tc.jp = j2Params();
            tc.S = j2ScenarioSet('small');
            tc.tmpDir = tempname;  mkdir(tc.tmpDir);
            mE = buildJ2Excitation(struct('modelName','j2_excite_test','slxPath',fullfile(tc.tmpDir,'j2_excite_test.slx')));
            mC = buildJ2ClosedLoop(struct('modelName','j2_closed_test','slxPath',fullfile(tc.tmpDir,'j2_closed_test.slx')));
            tc.models = struct('excite',mE,'closed',mC);
        end
    end

    methods (TestClassTeardown)
        function cleanup(tc)
            for m = {tc.models.excite, tc.models.closed}
                if bdIsLoaded(m{1}), close_system(m{1},0); end
            end
            if isfolder(tc.tmpDir), rmdir(tc.tmpDir,'s'); end
        end
    end

    methods (Test)
        function runsAndRespectsLimits(tc, scenarioName)
            s = tc.S(strcmp({tc.S.name}, scenarioName));
            tc.assertNumElements(s, 1);
            raw = runJ2Scenario(s, tc.models);
            jp = tc.jp;
            tc.verifyTrue(all(isfinite(raw.q)) && all(isfinite(raw.dq)) && all(isfinite(raw.tau)));
            tc.verifyEqual(raw.t(1), 0);
            tc.verifyEqual(diff(raw.t([1 2])), 1/8000, 'AbsTol', 1e-12);          % 8 kHz
            tc.verifyLessThan(max(abs(raw.tau)), jp.tauPeak + 1e-6);               % トルク飽和
            tc.verifyLessThan(max(abs(raw.dq)), jp.qdMax);                         % 速度制限内
            tc.verifyGreaterThan(min(raw.q), jp.qMin - deg2rad(3));                % めり込み 3° 以内
            tc.verifyLessThan(max(raw.q), jp.qMax + deg2rad(3));
            tc.verifyEqual(raw.q(1), raw.q0, 'AbsTol', 1e-6);                      % 初期状態が反映
            tc.verifyEqual(raw.dq(1), raw.dq0, 'AbsTol', 1e-3);

            switch s.pattern
                case {'bln_normal','bln_high','bln_low','chirp','step','micro'}
                    % ソフトバリアが働き、リミット拘束に達しない
                    tc.verifyEqual(raw.contact, 0);
            end
            if strcmp(s.pattern,'micro')
                tc.verifyLessThan(rad2deg(max(raw.q) - min(raw.q)), 10);           % 摩擦域で大きく動かない
            end
            if strcmp(s.pattern,'step')
                tc.verifyGreaterThan(max(abs(raw.tau)), 0.99*jp.tauPeak);          % 飽和まで駆動している
            end
            if startsWith(s.pattern, 'contact_')
                tol = deg2rad(0.1);
                tc.verifyGreaterThan(raw.contact, 0);                                  % 実際にリミットへ接触した
                switch s.pattern
                    case 'contact_fall'     % 重力方向（開始側）のリミットに衝突
                        if s.side > 0, tc.verifyGreaterThanOrEqual(max(raw.q), jp.qMax - tol);
                        else, tc.verifyLessThanOrEqual(min(raw.q), jp.qMin + tol); end
                    case 'contact_drive'    % 押し付け → 反転 → 反対側へ。両側に接触
                        tc.verifyGreaterThanOrEqual(max(raw.q), jp.qMax - tol);
                        tc.verifyLessThanOrEqual(min(raw.q), jp.qMin + tol);
                end
                tc.verifyLessThan(max(abs(raw.dq)), jp.qdMax);                         % 衝突速度も上限内
            end
            if strcmp(s.type,'hold')
                tc.verifyLessThan(rad2deg(max(abs(raw.q - raw.qref))), 3);         % 外乱下でも参照の近傍
            end
            if strcmp(s.type,'ptp')
                tc.verifyLessThan(rad2deg(max(abs(raw.q - raw.qref))), 0.05);      % PTP の追従誤差
            end
        end

        function barrierPreventsLimitOvershoot(tc)
            % バリアなしでは ±400 N·m 往復で速度 25 rad/s・リミット超過だった入力（Phase 0 実測）
            jp = tc.jp;
            T = 2;  tt = (0:1/8000:T)';  u = 400*sign(sin(2*pi*tt/0.6));
            setJ2Barrier(tc.models.excite, 0.85, 0.60);
            o = runJ2Sim(tc.models.excite, tt, u, deg2rad(-90), 0);
            tc.verifyLessThan(max(abs(o.dq)), jp.qdMax);
            tc.verifyGreaterThan(min(o.q), jp.qMin);
            tc.verifyLessThan(max(o.q), jp.qMax);
        end

        function scenarioIsReproducible(tc)
            s = tc.S(strcmp({tc.S.name}, 'bln_high_001'));
            a = runJ2Scenario(s, tc.models);
            b = runJ2Scenario(s, tc.models);
            tc.verifyEqual(a.q, b.q, 'AbsTol', 1e-12);
            tc.verifyEqual(a.tau, b.tau, 'AbsTol', 1e-9);
        end
    end
end
