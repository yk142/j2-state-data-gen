function tauFf = j2Feedforward(qr, dqr, ddqr)
%J2FEEDFORWARD J2 の解析モデルによるフィードフォワードトルク。
%   tauFf = J2FEEDFORWARD(qr, dqr, ddqr)
%   M·ddq = τ + mgL·sin θ − Bv·dq − Fc·tanh(dq/ε) の逆: 
%       τ_ff = M·ddq_ref − mgL·sin θ_ref + Bv·dq_ref + Fc·tanh(dq_ref/ε)
%   （重力トルクは θ>0 で θ を増やす向き、保持トルクは −mgL·sin θ。Phase 1 で検証済み）
jp = j2Params();
tauFf = jp.M*ddqr - jp.mgL*sin(qr) + jp.Bv*dqr + jp.Fc*tanh(dqr/jp.eps);
end
