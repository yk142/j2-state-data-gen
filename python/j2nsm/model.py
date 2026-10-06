"""Neural State Model（1 ステップ遷移の MLP）。

入力 X_spec = [sin q, cos q, dq/DTHETA_MAX, tau/TAU_MAX]（引継ぎ資料の形式）
出力（差分残差）:
  通常     : [Δq/D_THETA_MAX, Δdq/D_DTHETA_MAX]
  構造付き : [Δdq/D_DTHETA_MAX] のみ。Δq は運動学 dt·(dq + Δdq/2)（台形則）で求める
入力は train の平均・標準偏差で、出力は train の標準偏差で標準化する（バッファとして保存）。
"""
import torch
import torch.nn as nn


class NSM(nn.Module):
    def __init__(self, scale, dt, width=256, depth=3, structured=False,
                 x_mean=None, x_std=None, y_mean=None, y_std=None):
        super().__init__()
        self.scale, self.dt, self.structured = dict(scale), float(dt), bool(structured)
        self.width, self.depth = int(width), int(depth)
        out_dim = 1 if structured else 2
        layers, d = [], 4
        for _ in range(depth):
            layers += [nn.Linear(d, width), nn.SiLU()]
            d = width
        layers.append(nn.Linear(d, out_dim))
        self.net = nn.Sequential(*layers)
        z4, o4, z2, o2 = torch.zeros(4), torch.ones(4), torch.zeros(out_dim), torch.ones(out_dim)
        self.register_buffer("x_mean", z4 if x_mean is None else torch.as_tensor(x_mean, dtype=torch.float32))
        self.register_buffer("x_std", o4 if x_std is None else torch.as_tensor(x_std, dtype=torch.float32))
        self.register_buffer("y_mean", z2 if y_mean is None else torch.as_tensor(y_mean, dtype=torch.float32))
        self.register_buffer("y_std", o2 if y_std is None else torch.as_tensor(y_std, dtype=torch.float32))

    def config(self):
        return {"width": self.width, "depth": self.depth, "structured": self.structured,
                "scale": self.scale, "dt": self.dt}

    def forward(self, xs):
        """X_spec → 標準化した出力。"""
        return self.net((xs - self.x_mean) / self.x_std)

    def target(self, ys):
        """Y_spec（2 列）→ このモデルの学習目標（構造付きなら Δdq の 1 列）。"""
        return ys[:, 1:2] if self.structured else ys

    def predict_spec(self, xs):
        return self.forward(xs) * self.y_std + self.y_mean

    def step_physical(self, q, dq, tau):
        """物理量（rad, rad/s, N·m）の状態を 1 ステップ進める。q, dq, tau は同じ形のテンソル。"""
        s = self.scale
        xs = torch.stack([torch.sin(q), torch.cos(q), dq / s["DTHETA_MAX"], tau / s["TAU_MAX"]], dim=-1).float()
        y = self.predict_spec(xs).to(q.dtype)
        d_dq = y[..., -1] * s["D_DTHETA_MAX"]
        if self.structured:
            d_q = self.dt * (dq + 0.5 * d_dq)
        else:
            d_q = y[..., 0] * s["D_THETA_MAX"]
        return q + d_q, dq + d_dq

    def delta_physical(self, xs):
        """X_spec から物理量の (Δq, Δdq) を返す（1 ステップ誤差の評価用）。dq は X_spec の 3 列目から復元する。"""
        s = self.scale
        y = self.predict_spec(xs)
        d_dq = y[:, -1] * s["D_DTHETA_MAX"]
        if self.structured:
            dq = xs[:, 2] * s["DTHETA_MAX"]
            d_q = self.dt * (dq + 0.5 * d_dq)
        else:
            d_q = y[:, 0] * s["D_THETA_MAX"]
        return d_q, d_dq
