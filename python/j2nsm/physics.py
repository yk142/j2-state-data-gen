"""J2 の解析モデル（基準線）。

    M·ddq = τ + mgL·sin θ − Bv·dq − Fc·tanh(dq/ε)

を、トルクを 1 ステップの間一定（ゼロ次保持）として RK4 で積分する。ストッパ（リミット）はモデルに無いため、
接触を含む区間の比較用に、可動範囲の端で位置を止め外向きの速度を 0 にする単純な非弾性ストッパを任意で付ける。
"""
import torch


class Physics:
    def __init__(self, p, dt, substeps=8):
        self.M, self.mgL, self.Fc, self.Bv, self.eps = p["M"], p["mgL"], p["Fc"], p["Bv"], p["eps"]
        self.qMin, self.qMax = p["qMin"], p["qMax"]
        self.dt, self.substeps = dt, substeps

    def accel(self, q, dq, tau):
        return (tau + self.mgL * torch.sin(q) - self.Bv * dq - self.Fc * torch.tanh(dq / self.eps)) / self.M

    def step(self, q, dq, tau, stop=True):
        """1 ステップ（dt）進める。q, dq, tau は同じ形のテンソル。"""
        h = self.dt / self.substeps
        for _ in range(self.substeps):
            k1q, k1v = dq, self.accel(q, dq, tau)
            k2q, k2v = dq + 0.5 * h * k1v, self.accel(q + 0.5 * h * k1q, dq + 0.5 * h * k1v, tau)
            k3q, k3v = dq + 0.5 * h * k2v, self.accel(q + 0.5 * h * k2q, dq + 0.5 * h * k2v, tau)
            k4q, k4v = dq + h * k3v, self.accel(q + h * k3q, dq + h * k3v, tau)
            q = q + h / 6 * (k1q + 2 * k2q + 2 * k3q + k4q)
            dq = dq + h / 6 * (k1v + 2 * k2v + 2 * k3v + k4v)
            if stop:
                hi, lo = q > self.qMax, q < self.qMin
                q = torch.where(hi, torch.full_like(q, self.qMax), q)
                q = torch.where(lo, torch.full_like(q, self.qMin), q)
                dq = torch.where(hi & (dq > 0), torch.zeros_like(dq), dq)
                dq = torch.where(lo & (dq < 0), torch.zeros_like(dq), dq)
        return q, dq
