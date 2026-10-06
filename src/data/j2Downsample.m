function D = j2Downsample(raw, fsIn, fsOut)
%J2DOWNSAMPLE 8 kHz の生データを 1 kHz の遷移データへ変換する。
%   D = J2DOWNSAMPLE(raw, fsIn, fsOut)   raw: .q .dq .tau（8 kHz）、.qref（任意）
%
%   目的は 1 ステップ遷移 x_{k+1} = f(x_k, u_k) の教師データを作ること。そのため:
%     - **状態 (q, dq) は瞬時値をそのまま間引く**（r = fsIn/fsOut 点ごと）。状態はある時刻の
%       値であり、ローパスをかけると状態間の物理的な関係（遷移）を歪めるため、
%       アンチエイリアス・フィルタは使わない（親資産の c8Downsample は加速度の差分用で目的が違う）。
%     - **トルク u_k は区間 [t_k, t_{k+1}) の平均**（8 kHz の r サンプルの平均）。
%       1 ms の間に印加されたトルクの実効値であり、ステップ入力・飽和・バリアの切替でも
%       遷移と整合する。区間先頭の瞬時値は tauInst に残す。
%   出力（K = 遷移数）:
%     D.q, D.dq   [K+1 x 1] 瞬時値（q(k+1) が x_k の次状態）
%     D.tau       [K x 1] 区間平均トルク、D.tauInst [K x 1] 区間先頭の瞬時値
%     D.qref      [K+1 x 1]（参照があれば）、D.fs 出力周波数
r = fsIn / fsOut;
assert(abs(r - round(r)) < 1e-9, 'j2Downsample:ratio', '間引き率 %g が整数ではありません', r);
r = round(r);
N = numel(raw.q);                        % 8 kHz サンプル数（t = 0 ... (N-1)/fsIn）
K = floor((N-1) / r);
idx = 1 + r*(0:K);                       % 瞬時値を取る 8 kHz の添字
D.q  = double(raw.q(idx));
D.dq = double(raw.dq(idx));
tau = double(raw.tau(:));
blk = reshape(tau(1:r*K), r, K);          % 区間 k は 8 kHz のサンプル r(k-1)+1 ... rk
D.tau     = mean(blk, 1)';
D.tauInst = tau(1 + r*(0:K-1));
D.qref = [];
if isfield(raw,'qref') && ~isempty(raw.qref), D.qref = double(raw.qref(idx)); end
D.fs = fsOut;
end
