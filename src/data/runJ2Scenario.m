function raw = runJ2Scenario(s, models)
%RUNJ2SCENARIO 1 シナリオを 8 kHz で実行し、生データを返す。
%   raw = RUNJ2SCENARIO(s, models)
%     s       j2ScenarioSet の要素
%     models  .excite（励振モデル名）、.closed（閉ループモデル名）。ロード済みであること
%   raw: .t .q .dq .tau（8 kHz、実印加トルク）、.qref .dqref（閉ループのみ、無ければ空）、
%        .q0 .dq0（実際の初期状態）、.contact（リミット拘束 = 可動範囲の端から 0.1° 以内または外側の点の割合）
%
%   type ごとの内容:
%     excitation  開ループ励振＋ソフトリミットバリア（pattern で信号が変わる）
%     hold        PD＋FF で重要姿勢まわりの滑らかな参照に追従し、BLN 外乱トルクを加える
%                 （倒立点 θ=0 近傍など、自由落下では取れない領域を安定化制御で確保）
%     ptp         PD＋FF による PTP 連続 GO 動作（ベンチマーク）
jp = j2Params();
dt = 1/jp.fsSim;
raw = struct('t',[],'q',[],'dq',[],'tau',[],'qref',[],'dqref',[],'q0',0,'dq0',0,'contact',0);

switch s.type
    case 'excitation'
        tVec = (0:dt:s.duration)';
        amp  = s.ampRel * jp.tauRef;
        switch s.pattern
            case 'step'
                rs = RandStream('twister','Seed',s.seed);
                u = j2Signal('step', tVec, jp.tauPeak, struct('tOn',0.2,'sign',1-2*(rand(rs)<0.5)));
            case 'chirp'
                u = j2Signal('chirp', tVec, amp, struct('f0',0.1,'f1',s.fMax));
            case 'micro'
                % 摩擦・粘性域: 初期姿勢の重力保持トルクのまわりに微小トルク
                u = -jp.mgL*sin(s.q0) + j2Signal('bln', tVec, amp, struct('fMax',s.fMax,'seed',s.seed));
            otherwise    % bln_normal / bln_high / bln_low / freefall
                u = j2Signal('bln', tVec, amp, struct('fMax',s.fMax,'seed',s.seed));
        end
        setJ2Barrier(models.excite, s.qSoftFrac, 0.60);
        o = runJ2Sim(models.excite, tVec, u, s.q0, s.dq0);
        raw.q0 = s.q0;  raw.dq0 = s.dq0;

    case 'hold'
        tVec = (0:dt:s.duration)';
        [qr, dqr, ddqr] = j2SmoothRef(tVec, s.posture, s.refAmp, 1.0, s.seed);
        qr = min(max(qr, jp.qMin + s.refMargin), jp.qMax - s.refMargin);  % 参照を可動範囲の内側に収める
        tauD = j2Signal('bln', tVec, s.ampRel*jp.tauRef, struct('fMax',s.fMax,'seed',s.seed+1));
        tauFf = j2Feedforward(qr, dqr, ddqr);
        o = runJ2Sim(models.closed, tVec, [qr dqr tauFf tauD], qr(1), 0);
        raw.qref = qr;  raw.dqref = dqr;  raw.q0 = qr(1);

    case 'ptp'
        rs = RandStream('twister','Seed',s.seed);
        WP = jp.qTrainMin + (jp.qTrainMax - jp.qTrainMin)*rand(rs, s.nWp, 1);
        T  = j2PtpSegTime(WP, s.speedScale, 0.5);
        % 各到達点で 0.2 s 停留（整定・定常偏差を評価できるように）
        WP2 = WP(1);  T2 = [];
        for k = 1:numel(T)
            WP2(end+1,1) = WP(k+1);  T2(end+1,1) = T(k);        %#ok<AGROW>
            WP2(end+1,1) = WP(k+1);  T2(end+1,1) = 0.2;         %#ok<AGROW>
        end
        tVec = (0:dt:min(sum(T2), s.duration))';
        [qr, dqr, ddqr] = j2PtpTrajectory(tVec, WP2, T2);
        tauFf = j2Feedforward(qr, dqr, ddqr);
        o = runJ2Sim(models.closed, tVec, [qr dqr tauFf zeros(size(tVec))], qr(1), 0);
        raw.qref = qr;  raw.dqref = dqr;  raw.q0 = qr(1);

    otherwise
        error('runJ2Scenario:type', '未知の type: %s', s.type);
end
raw.t = o.t;  raw.q = o.q;  raw.dq = o.dq;  raw.tau = o.tau;
if strcmp(s.type,'excitation'), raw.dq0 = s.dq0; end
raw.contact = mean(raw.q <= jp.qMin + deg2rad(0.1) | raw.q >= jp.qMax - deg2rad(0.1));
end
