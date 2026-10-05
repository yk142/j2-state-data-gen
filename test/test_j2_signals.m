classdef test_j2_signals < matlab.unittest.TestCase
%TEST_J2_SIGNALS 励振信号・Sobol 列・導出パラメータの単体テスト（シミュレーション不要）。

    methods (TestClassSetup)
        function setup(tc) %#ok<MANU>
            repo = fileparts(fileparts(mfilename('fullpath')));
            addpath(fullfile(repo,'config'), genpath(fullfile(repo,'src')));
            j2_setup_path();
        end
    end

    methods (Test)
        function sobolKnownPoints(tc)
            X = j2Sobol(8, 2);
            ref = [0 0; .5 .5; .75 .25; .25 .75; .375 .375; .875 .875; .625 .125; .125 .625];
            tc.verifyEqual(X, ref, 'AbsTol', 1e-12);
        end

        function sobolIsLowDiscrepancyAndInRange(tc)
            n = 256;
            X = j2Sobol(n, 2, 7);
            tc.verifyGreaterThanOrEqual(min(X(:)), 0);
            tc.verifyLessThan(max(X(:)), 1);
            % 8x8 の格子の占有: Sobol は 256 点でほぼ全セルを埋める（一様乱数より均一）
            occ = histcounts2(X(:,1), X(:,2), 0:1/8:1, 0:1/8:1);
            tc.verifyEqual(min(occ(:)), 4);            % 256/64=4 点ずつ（Sobol の 2D 層化）
            % seed が同じなら同一、違えば異なる
            tc.verifyEqual(j2Sobol(16,3,7), j2Sobol(16,3,7));
            tc.verifyNotEqual(j2Sobol(16,3,7), j2Sobol(16,3,8));
            for d = 1:4, tc.verifyEqual(size(j2Sobol(32,d)), [32 d]); end
        end

        function blnBandAndAmplitude(tc)
            fs = 8000;  t = (0:1/fs:10)';  fMax = 20;  amp = 100;
            u = j2Signal('bln', t, amp, struct('fMax',fMax,'seed',3));
            tc.verifyEqual(max(abs(u)), amp, 'RelTol', 1e-3);
            [pxx, f] = pwelch(u - mean(u), 4096, 2048, 4096, fs);
            inband  = trapz(f(f<=fMax*1.2), pxx(f<=fMax*1.2));
            outband = trapz(f(f>fMax*2),    pxx(f>fMax*2));
            tc.verifyLessThan(outband/inband, 1e-3);          % 帯域外のパワーは 0.1% 未満
            % 再現性
            u2 = j2Signal('bln', t, amp, struct('fMax',fMax,'seed',3));
            tc.verifyEqual(u, u2);
            u3 = j2Signal('bln', t, amp, struct('fMax',fMax,'seed',4));
            tc.verifyNotEqual(u, u3);
        end

        function chirpSweepsBand(tc)
            fs = 8000;  t = (0:1/fs:20)';
            u = j2Signal('chirp', t, 50, struct('f0',0.1,'f1',20));
            tc.verifyEqual(max(abs(u)), 50, 'RelTol', 1e-3);
            [pxx, f] = pwelch(u, 8192, 4096, 8192, fs);
            % 0.1〜20 Hz に全パワーの 99% 以上
            tc.verifyGreaterThan(trapz(f(f<=25), pxx(f<=25)) / trapz(f, pxx), 0.99);
        end

        function stepSignal(tc)
            t = (0:1/8000:1)';
            u = j2Signal('step', t, 200, struct('tOn',0.25,'sign',-1));
            tc.verifyEqual(u(t < 0.25), zeros(nnz(t < 0.25),1));
            tc.verifyEqual(unique(u(t >= 0.25)), -200);
        end

        function derivedParams(tc)
            q = j2Params();
            tc.verifyEqual(q.M,   6.086, 'AbsTol', 0.01);       % 実効慣性（armature 込み、シミュレーション実測 6.084）
            tc.verifyEqual(q.mgL, 73.25, 'AbsTol', 0.05);
            tc.verifyEqual(q.tauPeak, 3*q.tauRated, 'RelTol', 1e-12);
            tc.verifyLessThan(q.qMin, q.qTrainMin);
            tc.verifyGreaterThan(q.qMax, q.qTrainMax);
        end
    end
end
