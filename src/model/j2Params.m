function q = j2Params()
%J2PARAMS J2 単軸（固定姿勢 J3=75.0684°）の導出パラメータを返す。
%   q = J2PARAMS()
%   親資産 c8_params と c8RNEA から、J2 単軸の解析モデル
%       M·ddq = τ − τ_hold(θ) − Bv·dq − Fc·tanh(dq/ε),   τ_hold = −mgL·sin θ
%   の係数を求める（フィードフォワード・ゲイン設計・励振振幅の基準に使う）。
%
%   出力フィールド（単位は SI）:
%     M        実効慣性 [kg*m^2]（armature 込み。c8RNEA は armature を含む）
%     mgL      重力トルク振幅 [N*m]（保持トルク = −mgL·sin θ）
%     Fc, Bv, eps   クーロン摩擦 [N*m]、粘性 [N*m/(rad/s)]、tanh 幅 [rad/s]
%     qMin, qMax, qdMax        可動範囲 [rad]、最大角速度 [rad/s]
%     tauRated, tauPeak        関節側の定格・瞬時最大トルク [N*m]
%     tauRef   励振振幅の基準トルク（= tauRated。引継ぎ資料の TAU_MAX 相当）
%     qTrainMin, qTrainMax     訓練範囲（可動域の ±70%）[rad]
%     qFixed   固定軸の角度 [1x6 rad]
c = j2_setup_path();
p = c8_params(c.meshDir);
J = c.jointIdx;
z = zeros(1,6);  qf = c.qFixed;
rn = @(qv,ddq) subsref(c8RNEA(qv, z, ddq, p), struct('type','()','subs',{{J}}));
e = z;  e(J) = 1;
q.M   = rn(qf, e) - rn(qf, z);
qh = qf;  qh(J) = -pi/2;
q.mgL = rn(qh, z);                          % 保持トルクは −mgL·sin θ なので θ=−90° で +mgL
q.Fc  = p.frictionFc(J);   q.Bv = p.frictionBv(J);   q.eps = p.frictionEps(J);
q.qMin = p.qMin(J);  q.qMax = p.qMax(J);  q.qdMax = p.qdMax(J);
q.tauRated = p.tauRated(J);  q.tauPeak = p.tauPeak(J);
q.tauRef = q.tauRated;
q.qTrainMin = p.qTrainMin(J);  q.qTrainMax = p.qTrainMax(J);
q.qFixed = qf;
q.fsSim = c.fsSim;  q.fsData = c.fsData;
end
