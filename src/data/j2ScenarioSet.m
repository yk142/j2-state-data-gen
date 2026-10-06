function S = j2ScenarioSet(preset)
%J2SCENARIOSET 学習データ生成のシナリオ表を返す（引継ぎ資料の収集計画を J2 実値で再計算）。
%   S = J2SCENARIOSET(preset)   preset = 'full'（既定）| 'small'（テスト・動作確認用）
%
%   S は構造体配列。主なフィールド:
%     name, pattern, type ('excitation'|'hold'|'ptp'), phase (1|2|3 = リミット接触|'benchmark'), split
%     seed        シナリオごとの乱数シード（再現性。NFR-1）
%     duration    継続時間 [s]（ptp は軌道長で決まるため上限）
%     ampRel      励振振幅 / tauRef（tauRef = 定格トルク）、fMax 帯域 [Hz]
%     q0, dq0     初期状態 [rad], [rad/s]（hold/ptp は参照軌道から決まる）
%     qSoftFrac   ソフトリミットバリアの安全帯比率（excitation のみ）
%     posture     hold の中心姿勢 [rad]、refAmp 参照の振れ幅 [rad]、refMargin 可動範囲端からの余裕 [rad]
%     speedScale/nWp  ptp のパラメータ
%
%   振幅の基準 tauRef は定格トルク（286.8 N·m）。引継ぎ資料の TAU_MAX=300 N·m 相当。
%   リミット近傍（nearlimit）は閉ループ保持で、外乱により時々リミットに接触する。
%   リミット拘束サンプル（リミットから 0.1° 以内）は runJ2Scenario が割合を返し、
%   データセット側でフラグを付ける。
if nargin < 1, preset = 'full'; end
jp = j2Params();
c  = j2_config();
d2r = pi/180;

% 初期状態の範囲: ソフトバリア帯（0.85）の内側、リミットから 2° 以上内側
qLo = jp.qMin + (1-0.85)*(jp.qMax-jp.qMin)/2 + 2*d2r;
qHi = jp.qMax - (1-0.85)*(jp.qMax-jp.qMin)/2 - 2*d2r;
dqI = 0.6 * jp.qdMax * 0.8;                          % 初期角速度の範囲（速度バリア帯の 80%）

% ---- パターン定義: {pattern, type, phase, 件数(full), 件数(small), duration(full), duration(small), ampRel, fMax} ----
postures = deg2rad([0 30 60 -45 -90 -120 -150]);
switch preset
    case 'full'
        P = { 'bln_normal','excitation',1, 60, 10, 0.5, 20
              'bln_high',  'excitation',1, 30,  5, 1.0, 20
              'bln_low',   'excitation',1, 30, 15, 0.2,  5
              'freefall',  'excitation',1,  6,  3, 0.05, 5
              'chirp',     'excitation',2, 10, 20, 0.5, 20
              'step',      'excitation',2, 20,  2, 1.0, 0
              'micro',     'excitation',2, 20, 10, 0.05, 5 };
        nHold = 20;  durHold = 5;  nNear = 10;  durNear = 10;  nPtp = 3;  speeds = [0.2 0.4 0.6 0.8];  nWp = 6;
        nContact = [12 24 20];
    case 'small'
        P = { 'bln_normal','excitation',1, 2, 2, 0.5, 20
              'bln_high',  'excitation',1, 1, 1, 1.0, 20
              'bln_low',   'excitation',1, 1, 2, 0.2,  5
              'freefall',  'excitation',1, 1, 1.5, 0.05, 5
              'chirp',     'excitation',2, 1, 3, 0.5, 20
              'step',      'excitation',2, 1, 1, 1.0, 0
              'micro',     'excitation',2, 1, 2, 0.05, 5 };
        nHold = 1;  durHold = 2;  nNear = 2;  durNear = 2;  nPtp = 1;  speeds = [0.3 0.7];  nWp = 3;
        nContact = [2 2 2];
    otherwise
        error('j2ScenarioSet:preset', '未知の preset: %s', preset);
end

% 励振シナリオの総数 → Sobol 点を一括で割り当てる（実行順に依存しない）
nExc = sum([P{:,4}]);
U = j2Sobol(nExc, 2, c.seedBase);
S = struct('name',{},'pattern',{},'type',{},'phase',{},'split',{},'seed',{},'duration',{}, ...
    'ampRel',{},'fMax',{},'q0',{},'dq0',{},'qSoftFrac',{},'posture',{},'speedScale',{},'nWp',{}, ...
    'refAmp',{},'refMargin',{},'dqSoftFrac',{},'ampRel2',{},'tFlip',{},'side',{});
k = 0;  iExc = 0;
for r = 1:size(P,1)
    for i = 1:P{r,4}
        k = k + 1;  iExc = iExc + 1;
        s = emptyScenario();
        s.pattern = P{r,1};  s.type = P{r,2};  s.phase = P{r,3};
        s.name = sprintf('%s_%03d', s.pattern, i);
        s.seed = c.seedBase + 1000*r + i;
        s.duration = P{r,5};  s.ampRel = P{r,6};  s.fMax = P{r,7};
        s.q0  = qLo + U(iExc,1)*(qHi - qLo);
        s.dq0 = (2*U(iExc,2) - 1) * dqI;
        s.qSoftFrac = 0.85;
        if strcmp(s.pattern,'freefall'), s.dq0 = s.dq0 * 0.5; end
        S(end+1) = s; %#ok<AGROW>
    end
end

% ---- 姿勢保持 + 外乱（安定化制御、重要 7 姿勢）----
for ip = 1:numel(postures)
    for i = 1:nHold
        s = emptyScenario();
        s.pattern = sprintf('hold_%+04d', round(rad2deg(postures(ip))));
        s.type = 'hold';  s.phase = 2;
        s.name = sprintf('%s_%03d', s.pattern, i);
        s.seed = c.seedBase + 100000 + 1000*ip + i;
        s.duration = durHold;  s.ampRel = 0.3;  s.fMax = 20;
        s.posture = postures(ip);
        s.refAmp = deg2rad(10);  s.refMargin = deg2rad(3);
        S(end+1) = s; %#ok<AGROW>
    end
end

% ---- リミット近傍（閉ループ保持。参照はリミットの 4° 内側、振れ幅 ±3°、強めの外乱で時々接触）----
% 開ループの押し付けは不適だった: +側・−側とも重力がリミット方向に働くため、
% 押さなくてもリミットに張り付いたまま動かず、情報量が乏しい（接触 84〜89%）。
for i = 1:nNear
    s = emptyScenario();
    s.pattern = 'nearlimit';  s.type = 'hold';  s.phase = 2;
    s.name = sprintf('nearlimit_%03d', i);
    s.seed = c.seedBase + 150000 + i;
    side = 1 - 2*(mod(i,2)==0);                          % 奇数: 正側、偶数: 負側
    if side > 0, s.posture = jp.qMax - 4*d2r; else, s.posture = jp.qMin + 4*d2r; end
    s.duration = durNear;  s.ampRel = 0.5;  s.fMax = 20;
    s.refAmp = 3*d2r;  s.refMargin = 0.5*d2r;
    S(end+1) = s; %#ok<AGROW>
end

% ---- PTP 連続 GO ベンチマーク ----
for is = 1:numel(speeds)
    for i = 1:nPtp
        s = emptyScenario();
        s.pattern = sprintf('ptp_v%02d', round(100*speeds(is)));
        s.type = 'ptp';  s.phase = 'benchmark';
        s.name = sprintf('%s_%03d', s.pattern, i);
        s.seed = c.seedBase + 200000 + 1000*is + i;
        s.duration = 60;  s.speedScale = speeds(is);  s.nWp = nWp;
        S(end+1) = s; %#ok<AGROW>
    end
end

% ---- リミット接触（機械ストッパへの衝突・押し付け・離脱。励振モデルのバリアをリミットの外側へ移して接触を許す）----
% 追加は既存シナリオの後ろ（既存のシード・Sobol 割り当て・split を変えない）。パラメータはシナリオのシードから決める。
%   contact_fall : 自由落下で重力方向のリミットに衝突 → 反発・押し付けて静止（微小 BLN）
%   contact_drive: リミット方向へ定トルク（0.2〜0.5 τ_ref）で衝突・押し付け → tFlip で反転し反対側へ（両側に接触）
%   contact_bln  : 大振幅 BLN（1.0 τ_ref）にバイアス（0.15 τ_ref）を足してリミットへ繰り返し接触
% 速度バリアは最大角速度の 80%（衝突速度を速度上限内に収める）。静止サンプルが増えすぎないよう継続時間は短め。
kinds = {'contact_fall', 3, 5, 0.05; 'contact_drive', 4, 5, 0; 'contact_bln', 10, 5, 1.0};   % {pattern, 継続時間, fMax, ampRel}
for kd = 1:3
    for i = 1:nContact(kd)
        s = emptyScenario();
        s.pattern = kinds{kd,1};  s.type = 'excitation';  s.phase = 3;
        s.name = sprintf('%s_%03d', s.pattern, i);
        s.seed = c.seedBase + 300000 + 1000*kd + i;
        rs = RandStream('twister','Seed',s.seed);
        s.side = 1 - 2*(mod(i,2)==0);                       % 奇数: +側、偶数: −側（最初に向かう側）
        s.duration = kinds{kd,2};  s.fMax = kinds{kd,3};  s.ampRel = kinds{kd,4};
        s.qSoftFrac = 1.10;  s.dqSoftFrac = 0.80;           % 位置バリアはリミットの外側（実質無効）
        switch s.pattern
            case 'contact_fall'
                if s.side > 0, s.q0 = deg2rad(25 + 30*rand(rs)); else, s.q0 = -deg2rad(25 + 115*rand(rs)); end
                s.dq0 = (2*rand(rs) - 1) * 0.5;
            case 'contact_drive'
                s.q0 = qLo + rand(rs)*(qHi - qLo);  s.dq0 = 0;
                s.ampRel = 0.2 + 0.3*rand(rs);  s.ampRel2 = 0.2 + 0.3*rand(rs);  s.tFlip = 1.5 + 0.7*rand(rs);
            case 'contact_bln'
                s.q0 = qLo + rand(rs)*(qHi - qLo);  s.dq0 = (2*rand(rs) - 1) * dqI;
                s.ampRel2 = 0.15;                           % バイアス（側 side の向き）
        end
        S(end+1) = s; %#ok<AGROW>
    end
end

% ---- train / val / test をパターンごとに 70/15/15 で分割（シードで決定的）----
pats = unique({S.pattern}, 'stable');
for ip = 1:numel(pats)
    idx = find(strcmp({S.pattern}, pats{ip}));
    n = numel(idx);
    rs = RandStream('twister','Seed',c.seedBase + ip);
    ord = idx(randperm(rs, n));
    nTr = ceil(0.7*n);  nVa = floor(0.15*n);
    if n < 3, nTr = n; nVa = 0; end                   % 件数が少ないパターンは train に寄せる
    for j = 1:n
        if strcmp(S(ord(j)).phase,'benchmark'), S(ord(j)).split = 'benchmark';
        elseif j <= nTr, S(ord(j)).split = 'train';
        elseif j <= nTr + nVa, S(ord(j)).split = 'val';
        else, S(ord(j)).split = 'test';
        end
    end
end
end

function s = emptyScenario()
s = struct('name','','pattern','','type','','phase',1,'split','','seed',0,'duration',0, ...
    'ampRel',0,'fMax',0,'q0',0,'dq0',0,'qSoftFrac',0.85,'posture',0,'speedScale',0,'nWp',0, ...
    'refAmp',0,'refMargin',0,'dqSoftFrac',0.60,'ampRel2',0,'tFlip',0,'side',0);
end
