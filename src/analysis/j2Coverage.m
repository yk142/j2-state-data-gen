function M = j2Coverage(ds, outDir, opts)
%J2COVERAGE データセットのカバレッジを評価し、PNG と指標を出力する（要件 FR-8）。
%   M = J2COVERAGE(ds, outDir, opts)
%     ds      genJ2Dataset の出力（または data/j2_dataset_<preset>.mat を load したもの）
%     outDir  PNG の保存先（例: reports/issue-9）
%     opts.excludeAtLimit  リミット拘束の遷移を除いて評価する（既定 true。学習用フラットデータと同じ扱い）
%     opts.nBins           占有率の 1 次元あたりのビン数（既定 20。引継ぎ資料の n_bins）
%   出力 PNG:
%     coverage_timeseries_A/B/C.png   時系列（θ, θ̇, τ）。A: 広域励振、B: 構造的（ステップ・微小・保持・リミット近傍）、C: PTP ベンチマーク
%     coverage_phase_planes.png       位相平面 3 種（θ–θ̇, θ–τ, θ̇–τ）の密度。学習系（励振＋保持）
%     coverage_phase_compare.png      位相平面 3 種の学習系とベンチマークの比較
%     coverage_frequency.png          周波数カバレッジ（τ と θ̇ の PSD をパターン群別、全体、20 Hz の帯域線）
%     coverage_spectrogram.png        スペクトログラム（チャープと BLN）
%   M: 数値指標（3 次元占有率、2 次元占有率、ベンチマークが学習系の占有セルに入る割合、周波数の被覆帯域など）
if nargin < 3, opts = struct(); end
if ~isfield(opts,'excludeAtLimit'), opts.excludeAtLimit = true; end
if ~isfield(opts,'nBins'), opts.nBins = 20; end
if ~isfolder(outDir), mkdir(outDir); end
jp = j2Params();
S = ds.scen;
isBench = strcmp({S.split},'benchmark');

% ---- フラット化（学習系 = ベンチマーク以外の全 split、ベンチマーク）----
fo = struct('excludeAtLimit', opts.excludeAtLimit);
[Xl, Yl] = j2BuildFlat(ds, {'train','val','test'}, fo);     %#ok<ASGLU>
[Xb, Yb] = j2BuildFlat(ds, 'benchmark', fo);                %#ok<ASGLU>
M.nTransition = struct('learning', size(Xl,1), 'benchmark', size(Xb,1));

% ---- 数値指標 ----
nb = opts.nBins;
edges = {linspace(jp.qMin, jp.qMax, nb+1), linspace(-jp.qdMax, jp.qdMax, nb+1), linspace(-jp.tauPeak, jp.tauPeak, nb+1)};
H3 = histcnd(Xl, edges);
M.occupancy3D = nnz(H3) / numel(H3);
M.occupancy2D.theta_dtheta = occ2(Xl(:,1), Xl(:,2), edges{1}, edges{2});
M.occupancy2D.theta_tau    = occ2(Xl(:,1), Xl(:,3), edges{1}, edges{3});
M.occupancy2D.dtheta_tau   = occ2(Xl(:,2), Xl(:,3), edges{2}, edges{3});
M.occupancy1D.theta  = nnz(histcounts(Xl(:,1), edges{1}))/nb;
M.occupancy1D.dtheta = nnz(histcounts(Xl(:,2), edges{2}))/nb;
M.occupancy1D.tau    = nnz(histcounts(Xl(:,3), edges{3}))/nb;
M.bins = nb;
% ベンチマークの遷移が、学習系データの占有セルに入る割合（学習系がベンチマークの動作域を覆えているか）
if ~isempty(Xb)
    M.benchmarkInOccupied3D = inOccupied(Xb, H3, edges);
    M.benchmarkInOccupied2D.theta_dtheta = inOccupied2(Xb(:,[1 2]), Xl(:,[1 2]), edges([1 2]));
    M.benchmarkInOccupied2D.theta_tau    = inOccupied2(Xb(:,[1 3]), Xl(:,[1 3]), edges([1 3]));
    M.benchmarkInOccupied2D.dtheta_tau   = inOccupied2(Xb(:,[2 3]), Xl(:,[2 3]), edges([2 3]));
end
% θ の被覆（θ=0 近傍と両リミット近傍の点数）
M.thetaRegions.upright_pm10deg = nnz(abs(Xl(:,1)) < deg2rad(10));
M.thetaRegions.nearPosLimit_10deg = nnz(Xl(:,1) > jp.qMax - deg2rad(10));
M.thetaRegions.nearNegLimit_10deg = nnz(Xl(:,1) < jp.qMin + deg2rad(10));

% ---- 図 ----
plotTimeSeries(S, outDir, jp);
plotPhasePlanes(Xl, Xb, outDir, jp, M);
M.frequency = plotFrequency(S, outDir, jp);
plotSpectrogram(S, outDir, jp);
end

% ======================= 指標のヘルパ =======================
function H = histcnd(X, edges)
nb = numel(edges{1}) - 1;
ix = zeros(size(X,1), 3);
for d = 1:3
    ix(:,d) = discretize(X(:,d), edges{d});
end
ok = all(~isnan(ix), 2);
H = accumarray(ix(ok,:), 1, [nb nb nb]);
end

function o = occ2(a, b, ea, eb)
H = histcounts2(a, b, ea, eb);
o = nnz(H) / numel(H);
end

function f = inOccupied(Xb, H3, edges)
ix = zeros(size(Xb,1), 3);
for d = 1:3, ix(:,d) = discretize(Xb(:,d), edges{d}); end
ok = all(~isnan(ix), 2);
lin = sub2ind(size(H3), ix(ok,1), ix(ok,2), ix(ok,3));
f = mean(H3(lin) > 0);
end

function f = inOccupied2(Pb, Pl, edges)
H = histcounts2(Pl(:,1), Pl(:,2), edges{1}, edges{2});
ia = discretize(Pb(:,1), edges{1});  ib = discretize(Pb(:,2), edges{2});
ok = ~isnan(ia) & ~isnan(ib);
f = mean(H(sub2ind(size(H), ia(ok), ib(ok))) > 0);
end

% ======================= 時系列 =======================
function plotTimeSeries(S, outDir, jp)
groups = {
  'A', '広域励振',                {'bln_normal','bln_high','bln_low','freefall','chirp'}
  'B', '構造的（ステップ・微小・姿勢保持・リミット近傍）', {'step','micro','hold_+000','hold_-090','nearlimit'}
  'C', 'PTP 連続 GO ベンチマーク', {'ptp_v20','ptp_v40','ptp_v60','ptp_v80'} };
ptpPats = unique({S(startsWith({S.pattern},'ptp_')).pattern}, 'stable');   % 存在する PTP パターン（速度スケール順）
groups{3,3} = ptpPats(1:min(4,numel(ptpPats)));
for g = 1:size(groups,1)
    pats = groups{g,3};
    if isempty(pats), continue; end
    f = figure('Visible','off','Position',[0 0 1900 850],'Color','w');
    tl = tiledlayout(f, 3, numel(pats), 'TileSpacing','compact','Padding','compact');
    for j = 1:numel(pats)
        k = find(strcmp({S.pattern}, pats{j}), 1);
        if isempty(k), for i = 1:3, nexttile(tl,(i-1)*numel(pats)+j); axis off; end, continue; end
        s = S(k);  fs = jp.fsData;
        t = (0:numel(s.q)-1)'/fs;  tt = t(1:end-1);
        rows = {rad2deg(s.q), s.dq, s.tau};  tms = {t, t, tt};
        ylabs = {'\theta [deg]','\omega [rad/s]','\tau [N·m]'};
        for i = 1:3
            ax = nexttile(tl, (i-1)*numel(pats)+j);  hold(ax,'on');  grid(ax,'on');
            plot(ax, tms{i}, rows{i}, 'LineWidth', 1);
            if i == 1
                yline(ax, rad2deg(jp.qMin), 'r--');  yline(ax, rad2deg(jp.qMax), 'r--');  yline(ax, 0, 'g:');
                if ~isempty(s.qref), plot(ax, t, rad2deg(s.qref), 'k:'); end
                title(ax, sprintf('%s\n(%s)', s.pattern, s.name), 'Interpreter','none');  ylim(ax, [-175 80]);
            elseif i == 2, yline(ax, [-jp.qdMax jp.qdMax], 'r--');  ylim(ax, [-6 6]);
            else, yline(ax, [-jp.tauPeak jp.tauPeak], 'r--');  ylim(ax, [-950 950]);  xlabel(ax, 't [s]');
            end
            if j == 1, ylabel(ax, ylabs{i}); end
        end
    end
    title(tl, sprintf('時系列（θ, θ̇, τ）— %s  [1 kHz、赤破線: 可動範囲・速度上限・瞬時最大トルク、緑点線: θ=0、黒点線: 参照]', groups{g,2}));
    exportgraphics(f, fullfile(outDir, sprintf('coverage_timeseries_%s.png', groups{g,1})), 'Resolution', 100);
    close(f);
end
end

% ======================= 位相平面 =======================
function plotPhasePlanes(Xl, Xb, outDir, jp, M)
th = rad2deg(1);
f = figure('Visible','off','Position',[0 0 1800 560],'Color','w');
tl = tiledlayout(f, 1, 3, 'TileSpacing','compact','Padding','compact');
drawPlanes(tl, Xl, jp, 'パターン全体');
title(tl, sprintf('位相平面カバレッジ（学習系: 励振＋姿勢保持、%d 遷移、1 kHz）｜3 次元占有率 %.1f%%（%d^3 ビン）', ...
    size(Xl,1), 100*M.occupancy3D, M.bins));
exportgraphics(f, fullfile(outDir,'coverage_phase_planes.png'), 'Resolution', 100);  close(f);

if ~isempty(Xb)
    f = figure('Visible','off','Position',[0 0 1800 1100],'Color','w');
    tl = tiledlayout(f, 2, 3, 'TileSpacing','compact','Padding','compact');
    drawPlanes(tl, Xl, jp, '学習系');
    drawPlanes(tl, Xb, jp, 'PTP ベンチマーク');
    title(tl, sprintf('位相平面の比較: 上=学習系、下=PTP ベンチマーク｜ベンチマークの遷移のうち学習系の占有セルに入る割合: 3 次元 %.1f%%', ...
        100*M.benchmarkInOccupied3D));
    exportgraphics(f, fullfile(outDir,'coverage_phase_compare.png'), 'Resolution', 100);  close(f);
end
end

function drawPlanes(tl, X, jp, label)
d2r = 180/pi;
pairs = {[1 2], '\theta [deg]', '\omega [rad/s]'; [1 3], '\theta [deg]', '\tau [N·m]'; [2 3], '\omega [rad/s]', '\tau [N·m]'};
lim = {[jp.qMin jp.qMax]*d2r, [-jp.qdMax jp.qdMax], [-jp.tauPeak jp.tauPeak]};
sc  = [d2r 1 1];
nb = 60;
for p = 1:3
    ax = nexttile(tl);  a = pairs{p,1}(1);  b = pairs{p,1}(2);
    ea = linspace(lim{a}(1), lim{a}(2), nb+1);  eb = linspace(lim{b}(1), lim{b}(2), nb+1);
    H = histcounts2(X(:,a)*sc(a), X(:,b)*sc(b), ea, eb);
    imagesc(ax, ea([1 end]), eb([1 end]), log10(1 + H'));  set(ax,'YDir','normal');
    colormap(ax, flipud(gray(256)));  cb = colorbar(ax);  cb.Label.String = 'log_{10}(1+件数)';
    hold(ax,'on');  grid(ax,'on');
    if a == 1, xline(ax, 0, 'g:', 'LineWidth', 1.2);  xline(ax, [jp.qMin jp.qMax]*d2r, 'r--'); end
    if b == 2, yline(ax, [-jp.qdMax jp.qdMax], 'r--'); end
    if b == 3, yline(ax, [-jp.tauPeak jp.tauPeak], 'r--'); end
    if a == 2, xline(ax, [-jp.qdMax jp.qdMax], 'r--'); end
    xlabel(ax, pairs{p,2});  ylabel(ax, pairs{p,3});
    occ = nnz(H)/numel(H);
    title(ax, sprintf('%s: %s–%s（占有 %.0f%%, %d ビン格子）', label, strtok(pairs{p,2},' '), strtok(pairs{p,3},' '), 100*occ, nb));
end
end

% ======================= 周波数 =======================
function F = plotFrequency(S, outDir, jp)
fs = jp.fsData;
fam = @(p) regexprep(p, '^hold_.*', 'hold');                      % 姿勢保持 7 姿勢は 1 群にまとめる
pats = unique(cellfun(fam, {S.pattern}, 'UniformOutput', false), 'stable');
nfft = 8192;  fgrid = (0:nfft/2)' * fs/nfft;       % 0.12 Hz 刻み（短い軌跡は窓長が短く、実分解能は 1/窓長）
Ptau = zeros(numel(fgrid), numel(pats));  Pdq = Ptau;  cnt = zeros(1,numel(pats));
Pall = zeros(numel(fgrid),1);  Pdall = Pall;  nAll = 0;
for k = 1:numel(S)
    s = S(k);  if numel(s.tau) < 1024, continue; end
    gi = find(strcmp(pats, fam(s.pattern)));
    u = s.tau - mean(s.tau);  w = s.dq(1:end-1) - mean(s.dq(1:end-1));
    wl = min(nfft, 2^floor(log2(numel(u))));          % 窓は 2 のべき乗で軌跡長以下
    [pu, fu] = pwelch(u, hann(wl), floor(wl/2), nfft, fs);
    [pw, ~]  = pwelch(w, hann(wl), floor(wl/2), nfft, fs);
    Ptau(:,gi) = Ptau(:,gi) + pu;  Pdq(:,gi) = Pdq(:,gi) + pw;  cnt(gi) = cnt(gi) + 1;
    if ~strcmp(s.split,'benchmark'), Pall = Pall + pu;  Pdall = Pdall + pw;  nAll = nAll + 1; end
end
Ptau = Ptau ./ max(cnt,1);  Pdq = Pdq ./ max(cnt,1);  Pall = Pall/max(nAll,1);  Pdall = Pdall/max(nAll,1);

f = figure('Visible','off','Position',[0 0 1700 650],'Color','w');
tl = tiledlayout(f, 1, 2, 'TileSpacing','compact','Padding','compact');
cols = lines(numel(pats));
for q = 1:2
    ax = nexttile(tl);  hold(ax,'on');  grid(ax,'on');
    P = {Ptau, Pdq};  Pa = {Pall, Pdall};  nm = {'\tau [N·m]', '\omega [rad/s]'};
    for g = 1:numel(pats)
        if cnt(g) == 0, continue; end
        loglog(ax, fgrid(2:end), P{q}(2:end,g), 'Color', cols(g,:), 'LineWidth', 1, 'DisplayName', strrep(pats{g},'_','\_'));
    end
    loglog(ax, fgrid(2:end), Pa{q}(2:end), 'k', 'LineWidth', 2.5, 'DisplayName', '全体（学習系の平均）');
    xline(ax, 20, 'r--', '20 Hz', 'LineWidth', 1.2, 'HandleVisibility','off');
    set(ax,'XScale','log','YScale','log');  xlim(ax,[0.1 fs/2]);
    xlabel(ax,'周波数 [Hz]');  ylabel(ax, sprintf('PSD  [(%s)^2/Hz]', nm{q}));  title(ax, [nm{q} ' のパワースペクトル密度（パターン群別）']);
    if q == 2, legend(ax, 'Location','southwest', 'NumColumns', 2, 'FontSize', 8); end
end
title(tl, '周波数カバレッジ（1 kHz データ、Welch 法、群内平均）');
exportgraphics(f, fullfile(outDir,'coverage_frequency.png'), 'Resolution', 100);  close(f);

% 指標: 全体の PSD について、ピークから −40 dB 以内に収まる周波数範囲と、0.5〜20 Hz 帯の最小パワー
band = fgrid >= 0.5 & fgrid <= 20;               % 0.5 Hz 未満は短い軌跡の窓長では分解できない
for q = 1:2
    Pa = {Pall, Pdall};  p = Pa{q};  pk = max(p);
    above = fgrid(p >= pk*1e-4);
    nm = {'tau','dq'};
    F.(nm{q}).fCover40dB_Hz = max(above);
    F.(nm{q}).minBandPower_dB = 10*log10(min(p(band))/pk);       % 0.5〜20 Hz 内で最も弱い周波数（ピーク比）
    F.(nm{q}).fPeak_Hz = fgrid(find(p == pk, 1));
end
F.nTrajectories = nAll;
end

% ======================= スペクトログラム =======================
function plotSpectrogram(S, outDir, jp)
fs = jp.fsData;
ptpPats = unique({S(startsWith({S.pattern},'ptp_')).pattern}, 'stable');
pick = {'chirp', 'bln_normal', 'hold_-090'};
if ~isempty(ptpPats), pick{end+1} = ptpPats{ceil(numel(ptpPats)/2)}; end   % 中間の速度スケールの PTP
f = figure('Visible','off','Position',[0 0 1700 800],'Color','w');
tl = tiledlayout(f, 2, numel(pick), 'TileSpacing','compact','Padding','compact');
for j = 1:numel(pick)
    k = find(strcmp({S.pattern}, pick{j}), 1);
    if isempty(k), nexttile(tl); axis off; nexttile(tl,numel(pick)+j); axis off; continue; end
    s = S(k);
    for r = 1:2
        ax = nexttile(tl, (r-1)*numel(pick)+j);
        if r == 1, x = s.tau - mean(s.tau); nm = 'τ'; else, x = s.dq(1:end-1) - mean(s.dq(1:end-1)); nm = 'ω'; end
        wl = min(256, 2^floor(log2(numel(x)/4)));  wl = max(wl, 32);
        [sp, fq, tm] = spectrogram(x, hann(wl), floor(wl*0.75), 512, fs);
        imagesc(ax, tm, fq, 10*log10(abs(sp).^2 + eps));  set(ax,'YDir','normal');
        ylim(ax, [0 60]);  colorbar(ax);  clim(ax, [max(10*log10(abs(sp(:)).^2 + eps)) - 80, max(10*log10(abs(sp(:)).^2 + eps))]);
        xlabel(ax,'t [s]');  ylabel(ax,'Hz');  title(ax, sprintf('%s: %s（%s）', nm, s.pattern, s.name), 'Interpreter','none');
    end
end
title(tl, 'スペクトログラム（上: トルク τ、下: 角速度 ω、0〜60 Hz、dB）');
exportgraphics(f, fullfile(outDir,'coverage_spectrogram.png'), 'Resolution', 100);  close(f);
end
