classdef test_j2_plant < matlab.unittest.TestCase
%TEST_J2_PLANT J2 単軸プラントの物理検証（要件 FR-2、ゲート G1）。
%   実行: results = runtests('test/test_j2_plant.m')
%   プラントは一時フォルダに 1 回だけ生成する。

    properties
        cfg
        p
        mdl
        tmpDir
    end

    properties (TestParameter)
        holdPose = {-150, -90, -45, 30, 60};     % 静止保持を確認する J2 角度 [deg]
    end

    methods (TestClassSetup)
        function buildPlant(tc)
            repo = fileparts(fileparts(mfilename('fullpath')));
            addpath(fullfile(repo,'config'), fullfile(repo,'src'), genpath(fullfile(repo,'src')));
            tc.cfg = j2_setup_path();
            tc.p = c8_params(tc.cfg.meshDir);
            tc.tmpDir = tempname;  mkdir(tc.tmpDir);
            tc.mdl = buildJ2Plant(struct('modelName','j2_plant_test', ...
                'slxPath',fullfile(tc.tmpDir,'j2_plant_test.slx')));
        end
    end

    methods (TestClassTeardown)
        function cleanup(tc)
            if bdIsLoaded(tc.mdl), close_system(tc.mdl,0); end
            if isfolder(tc.tmpDir), rmdir(tc.tmpDir,'s'); end
        end
    end

    methods (Access = private)
        function g = gravJ2(tc, q2)
            % 他軸固定姿勢での J2 重力保持トルク（c8RNEA）
            q = tc.cfg.qFixed;  q(2) = q2;  z = zeros(1,6);
            t = c8RNEA(q, z, z, tc.p);
            g = t(2);
        end
    end

    methods (Test)
        function interfaceIsScalar(tc)
            tc.verifyEqual(numel(find_system(tc.mdl,'MaskType','Revolute Joint')), 1);
            tc.verifyEqual(numel(find_system(tc.mdl,'MaskType','Weld Joint')), 6);   % base_link→base + 固定 5 軸
            tc.verifyEqual(get_param(tc.mdl,'Solver'), 'ode14x');
            tc.verifyEqual(str2double(get_param(tc.mdl,'FixedStep')), 1/8000, 'AbsTol', 1e-15);
            o = simJ2(tc.mdl, [0;0.01], 0, 0);
            tc.verifyTrue(isvector(o.q) && isvector(o.dq));
        end

        function gravityIsPendulum(tc)
            % J3 固定姿勢で J2 の重力トルクが -mgL*sin(θ) に一致し、θ=0 が平衡点
            q2 = deg2rad(-158:1:65);
            g  = arrayfun(@(a) tc.gravJ2(a), q2);
            A = -sin(q2(:));  mgL = A \ g(:);
            % qFixed(3) は 4 桁（75.0684°）に丸めているため、残差は 1e-4 N*m 程度を許容
            tc.verifyLessThan(max(abs(g(:) - A*mgL)), 1e-4);            % N*m
            tc.verifyEqual(mgL, 73.25, 'AbsTol', 0.05);
            tc.verifyEqual(tc.gravJ2(0), 0, 'AbsTol', 1e-4);
        end

        function holdsWithGravityCompensation(tc, holdPose)
            % RNEA の重力補償トルクを与えれば静止する（Simscape の重力 = c8RNEA）
            q0 = deg2rad(holdPose);
            o = simJ2(tc.mdl, [0;0.5], tc.gravJ2(q0), q0, 0);
            tc.verifyLessThan(abs(rad2deg(o.q(end) - q0)), 1e-3);        % deg
            tc.verifyLessThan(max(abs(o.dq)), 1e-3);                     % rad/s
        end

        function uprightIsUnstable(tc)
            % θ=0 から ±5° ずれた点は、tau=0 で原点から遠ざかる（倒立＝不安定平衡）
            for q0d = [5 -5]
                o = simJ2(tc.mdl, [0;3], 0, deg2rad(q0d), 0);
                tc.verifyGreaterThan(sign(q0d)*rad2deg(o.q(end)), abs(q0d));
            end
            % 保持トルクの θ 微分は負（重力トルクが変位を拡大）
            h = 1e-4;
            tc.verifyLessThan((tc.gravJ2(h) - tc.gravJ2(-h))/(2*h), 0);
        end

        function frictionMasksNearUpright(tc)
            % クーロン摩擦 > 重力トルク の領域（|θ| 小）では、tau=0 でもほとんど動かない
            o = simJ2(tc.mdl, [0;3], 0, deg2rad(0.5), 0);
            tc.verifyLessThan(abs(rad2deg(o.q(end)) - 0.5), 0.5);
        end

        function freeFallReachesNegativeLimit(tc)
            % −90° から自由落下 → 負リミット付近（−158°）に到達して静止。数値発散なし
            o = simJ2(tc.mdl, [0;6], 0, deg2rad(-90), 0);
            qmin = tc.p.qMin(2);
            tc.verifyTrue(all(isfinite(o.q)) && all(isfinite(o.dq)));
            tc.verifyLessThan(min(o.q), qmin + deg2rad(2));              % リミット到達
            tc.verifyGreaterThan(min(o.q), qmin - deg2rad(3));           % めり込み 3° 以内
            tc.verifyLessThan(abs(o.dq(end)), 0.05);                     % 静止
            tc.verifyLessThan(max(abs(o.dq)), tc.p.qdMax(2));            % 速度上限内
        end

        function freeFallReachesPositiveLimit(tc)
            o = simJ2(tc.mdl, [0;4], 0, deg2rad(60), 0);
            qmax = tc.p.qMax(2);
            tc.verifyGreaterThan(max(o.q), qmax - deg2rad(1));
            tc.verifyLessThan(max(o.q), qmax + deg2rad(3));
            tc.verifyLessThan(abs(o.dq(end)), 0.05);
        end

        function initialVelocityIsApplied(tc)
            % 初期角速度 1 rad/s を与えると、直後の dq≈1 で q が増加する
            q0 = deg2rad(-90);
            o = simJ2(tc.mdl, [0;0.05], tc.gravJ2(q0), q0, 1.0);
            tc.verifyEqual(o.dq(1), 1.0, 'AbsTol', 1e-3);
            tc.verifyGreaterThan(o.q(end), q0 + 0.02);
        end

        function randomTorqueStaysFinite(tc)
            % ±150 N·m のランダム（0.05 s 保持）トルクを 1 s。NaN・発散なし
            rng(1);
            tVec = (0:0.05:1)';  tau = 150*(2*rand(size(tVec)) - 1);
            o = simJ2(tc.mdl, tVec, tau, deg2rad(-60), 0);
            tc.verifyTrue(all(isfinite(o.q)) && all(isfinite(o.dq)));
            tc.verifyLessThan(max(abs(o.dq)), 30);                       % バリア無しでも物理的範囲
        end
    end
end
