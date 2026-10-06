classdef test_j2_scenarioset < matlab.unittest.TestCase
%TEST_J2_SCENARIOSET シナリオ表と参照軌道の単体テスト（シミュレーション不要）。

    methods (TestClassSetup)
        function setup(tc) %#ok<MANU>
            repo = fileparts(fileparts(mfilename('fullpath')));
            addpath(fullfile(repo,'config'), genpath(fullfile(repo,'src')));
            j2_setup_path();
        end
    end

    methods (Test)
        function fullPlanCounts(tc)
            S = j2ScenarioSet('full');
            cnt = @(pat) nnz(strcmp({S.pattern}, pat));
            tc.verifyEqual([cnt('bln_normal') cnt('bln_high') cnt('bln_low') cnt('freefall')], [60 30 30 6]);
            tc.verifyEqual([cnt('chirp') cnt('step') cnt('micro') cnt('nearlimit')], [10 20 20 10]);
            tc.verifyEqual(nnz(strcmp({S.type},'hold')), 7*20 + 10);     % 重要 7 姿勢 × 20 + リミット近傍 10
            tc.verifyEqual(nnz(strcmp({S.type},'ptp')), 12);             % 4 速度 × 3
            tc.verifyEqual(numel(S), 338 + 56);                              % + リミット接触 56 本（#11）
            tc.verifyEqual([cnt('contact_fall') cnt('contact_drive') cnt('contact_bln')], [12 24 20]);
        end

        function seedsAndNamesAreUnique(tc)
            for pr = {'small','full'}
                S = j2ScenarioSet(pr{1});
                tc.verifyEqual(numel(unique([S.seed])), numel(S));
                tc.verifyEqual(numel(unique({S.name})), numel(S));
            end
        end

        function reproducible(tc)
            tc.verifyEqual(j2ScenarioSet('full'), j2ScenarioSet('full'));
        end

        function splitsAreConsistent(tc)
            S = j2ScenarioSet('full');
            isBench = strcmp({S.phase},'benchmark');
            tc.verifyTrue(all(strcmp({S(isBench).split},'benchmark')));
            tc.verifyTrue(all(ismember({S(~isBench).split}, {'train','val','test'})));
            frac = mean(strcmp({S(~isBench).split},'train'));
            tc.verifyEqual(frac, 0.7, 'AbsTol', 0.05);
        end

        function initialStatesAreInsideBounds(tc)
            jp = j2Params();  S = j2ScenarioSet('full');
            isExc = strcmp({S.type},'excitation');
            q0 = [S(isExc).q0];  dq0 = [S(isExc).dq0];
            tc.verifyGreaterThan(min(q0), jp.qMin);  tc.verifyLessThan(max(q0), jp.qMax);
            tc.verifyLessThan(max(abs(dq0)), 0.6*jp.qdMax);
            % リミット近傍シナリオ（閉ループ保持）の中心姿勢はリミットの 4° 内側、両側にある
            nl = S(strcmp({S.pattern},'nearlimit'));
            gap = min(abs([nl.posture] - jp.qMax), abs([nl.posture] - jp.qMin));
            tc.verifyEqual(gap, repmat(deg2rad(4), size(gap)), 'AbsTol', 1e-12);
            tc.verifyTrue(any([nl.posture] > 0) && any([nl.posture] < 0));
            tc.verifyTrue(all(strcmp({nl.type}, 'hold')));
        end

        function existingScenariosAreUnchangedByContactAdditions(tc)
            % 接触シナリオは既存の後ろに追加され、既存のシード・初期状態・split を変えない（キャッシュを再利用できる）
            S = j2ScenarioSet('full');
            isC = startsWith({S.pattern}, 'contact_');
            tc.verifyEqual(find(isC), (numel(S)-nnz(isC)+1 : numel(S)));      % 末尾に連続して並ぶ
            tc.verifyEqual(S(1).name, 'bln_normal_001');
            jp = j2Params();
            tc.verifyEqual(S(1).dqSoftFrac, 0.60);                            % 既存の励振は速度バリア 0.60 のまま
            tc.verifyEqual(S(1).qSoftFrac, 0.85);
            tc.verifyTrue(all(strcmp({S(isC).type}, 'excitation')) && all([S(isC).phase] == 3));
        end

        function contactScenariosAreWellFormed(tc)
            jp = j2Params();  S = j2ScenarioSet('full');
            C = S(startsWith({S.pattern}, 'contact_'));
            tc.verifyTrue(all([C.qSoftFrac] > 1));                            % 位置バリアはリミットの外側
            tc.verifyEqual(unique([C.dqSoftFrac]), 0.80);
            tc.verifyGreaterThan(min([C.q0]), jp.qMin);  tc.verifyLessThan(max([C.q0]), jp.qMax);
            tc.verifyEqual(sort(unique([C.side])), [-1 1]);
            % 落下は重力方向のリミットへ向かう側から始まる（θ>0 なら +側、θ<0 なら −側）
            F = C(strcmp({C.pattern}, 'contact_fall'));
            tc.verifyEqual(sign([F.q0]), [F.side]);
            tc.verifyGreaterThan(min(abs([F.q0])), deg2rad(24));              % 摩擦で動かない領域（約 ±19°）の外
            D = C(strcmp({C.pattern}, 'contact_drive'));
            tc.verifyTrue(all([D.ampRel] >= 0.2 & [D.ampRel] <= 0.5));
            tc.verifyTrue(all([D.tFlip] >= 1.5 & [D.tFlip] <= 2.2));
            tc.verifyLessThan(max([D.tFlip]), min([D.duration]));             % 反転は継続時間内
            % 両側で同数（奇数・偶数）
            tc.verifyEqual(nnz([C.side] > 0), nnz([C.side] < 0));
        end

        function ptpTrajectoryIsC2AndStopsAtWaypoints(tc)
            WP = deg2rad([-100; -60; -120]);  T = [1.0; 1.5];
            t = (0:1/8000:sum(T))';
            [q, dq, ddq] = j2PtpTrajectory(t, WP, T);
            tc.verifyEqual(q(1), WP(1), 'AbsTol', 1e-12);
            tc.verifyEqual(q(end), WP(end), 'AbsTol', 1e-9);
            [~, iw] = min(abs(t - T(1)));
            tc.verifyEqual(dq(iw), 0, 'AbsTol', 5e-3);                    % 経由点で停止
            tc.verifyEqual(max(abs(dq)), 1.875*abs(WP(2)-WP(1))/T(1), 'RelTol', 0.02);
            tc.verifyLessThan(max(abs(diff(q))), 0.01);                   % 連続
            tc.verifyLessThan(max(abs(diff(dq))), 0.01);
        end

        function ptpSegTimeRespectsLimits(tc)
            jp = j2Params();
            WP = deg2rad([-150; 30; -150]);
            for sp = [0.2 0.8]
                T = j2PtpSegTime(WP, sp, 0.5);
                pk = 1.875*abs(diff(WP))./T;
                tc.verifyLessThanOrEqual(max(pk), sp*jp.qdMax*(1+1e-9));
            end
        end

        function feedforwardMatchesAnalyticModel(tc)
            jp = j2Params();
            % 静止: 保持トルク = −mgL·sin θ
            tc.verifyEqual(j2Feedforward(deg2rad(-90), 0, 0), jp.mgL, 'AbsTol', 1e-9);
            tc.verifyEqual(j2Feedforward(0, 0, 0), 0, 'AbsTol', 1e-12);
            % 加速: M·ddq を加算
            tc.verifyEqual(j2Feedforward(0, 0, 2) - j2Feedforward(0, 0, 0), 2*jp.M, 'AbsTol', 1e-9);
        end

        function smoothRefIsBounded(tc)
            t = (0:1/8000:3)';
            [q, dq, ddq] = j2SmoothRef(t, deg2rad(-45), deg2rad(10), 1.0, 5);
            tc.verifyLessThanOrEqual(max(abs(q - deg2rad(-45))), deg2rad(10)*1.01);
            tc.verifyLessThan(max(abs(dq)), 3);                           % 1 Hz・10° 振幅なら十分小さい
            tc.verifyLessThan(max(abs(ddq)), 30);
        end
    end
end
