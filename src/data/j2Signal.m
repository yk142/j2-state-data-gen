function u = j2Signal(type, tVec, amp, opts)
%J2SIGNAL J2 励振用のトルク信号を生成する（8 kHz の時刻列 tVec 上）。
%   u = J2SIGNAL(type, tVec, amp, opts)
%     type  'bln'    バンドリミテッドノイズ（opts.fMax [Hz]、opts.seed）
%           'chirp'  周波数掃引（opts.f0, opts.f1 [Hz]、対数掃引）
%           'step'   ステップ（opts.tOn [s]、opts.sign=±1）
%     tVec  時刻 [N x 1] [s]、amp 振幅 [N*m]（信号の最大絶対値が amp になる）
%     u     [N x 1] [N*m]
%
%   BLN は 1 kHz で生成して 4 次 Butterworth（filtfilt、零位相）で帯域制限し、
%   pchip で tVec に補間する（8 kHz で直接フィルタすると fMax/fs が小さく数値的に不利なため）。
%   最大絶対値で正規化するので、実効値は amp より小さい（ガウス性のピーク係数は約 3〜4）。
tVec = tVec(:);
switch type
    case 'bln'
        assert(isfield(opts,'fMax') && isfield(opts,'seed'), 'j2Signal:opts', 'bln には fMax と seed が必要です');
        fsG = 1000;
        T = tVec(end) - tVec(1);
        nG = max(ceil(T*fsG) + 1, 64);
        tG = tVec(1) + (0:nG-1)'/fsG;
        rs = RandStream('twister','Seed',opts.seed);
        w = randn(rs, nG, 1);
        [b, a] = butter(4, opts.fMax/(fsG/2));
        y = filtfilt(b, a, w);
        y = y / max(abs(y));
        u = amp * interp1(tG, y, tVec, 'pchip', 'extrap');
        u = max(min(u, amp), -amp);          % 補間のオーバーシュートで amp を超えないように
    case 'chirp'
        t = tVec - tVec(1);
        T = t(end);
        u = amp * chirp(t, opts.f0, T, opts.f1, 'logarithmic');
    case 'step'
        sg = 1;  if isfield(opts,'sign'), sg = opts.sign; end
        u = amp * sg * double(tVec >= opts.tOn);
    otherwise
        error('j2Signal:type', '未知の type: %s', type);
end
end
