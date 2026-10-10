"""Neural State Model（1 ステップ遷移の MLP）。NSM（純粋な NN）と ResidualNSM（解析モデル＋NN の残差）。

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

    kind = "nsm"
    residual = False

    def config(self):
        return {"kind": self.kind, "width": self.width, "depth": self.depth, "structured": self.structured,
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

    def delta_physical(self, xs, x_raw=None):
        """X_spec から物理量の (Δq, Δdq) を返す（1 ステップ誤差の評価用）。dq は X_spec の 3 列目から復元する。
        x_raw（[q dq tau]）は残差型との共通インタフェースのための引数で、ここでは使わない。"""
        s = self.scale
        y = self.predict_spec(xs)
        d_dq = y[:, -1] * s["D_DTHETA_MAX"]
        if self.structured:
            dq = xs[:, 2] * s["DTHETA_MAX"]
            d_q = self.dt * (dq + 0.5 * d_dq)
        else:
            d_q = y[:, 0] * s["D_THETA_MAX"]
        return d_q, d_dq


class ResidualNSM(nn.Module):
    """解析モデル＋NN の残差型。

        次状態 = 解析モデルの 1 ステップ予測 + NN の補正

    土台はストッパ無しの解析モデル（RK4）。非弾性ストッパ版はストッパを無視するより悪いため使わない。
    NN の入力は X_spec（4 列）＋解析モデルの予測（spec 単位の Δq, Δdq の 2 列）で、出力は残差
    （真の Y_spec − 解析モデルの Y_spec）を train の標準偏差で標準化したもの。
    最後の層を 0 で初期化するので、学習の開始時点は解析モデルそのものになる。
    """

    kind = "residual"
    residual = True
    structured = False

    def __init__(self, scale, dt, phys_params, width=256, depth=3, x_mean=None, x_std=None, y_std=None):
        super().__init__()
        self.scale, self.dt = dict(scale), float(dt)
        self.phys_params = dict(phys_params)
        self.width, self.depth = int(width), int(depth)
        layers, d = [], 6
        for _ in range(depth):
            layers += [nn.Linear(d, width), nn.SiLU()]
            d = width
        last = nn.Linear(d, 2)
        nn.init.zeros_(last.weight)
        nn.init.zeros_(last.bias)
        layers.append(last)
        self.net = nn.Sequential(*layers)
        self.register_buffer("x_mean", torch.zeros(6) if x_mean is None else torch.as_tensor(x_mean, dtype=torch.float32))
        self.register_buffer("x_std", torch.ones(6) if x_std is None else torch.as_tensor(x_std, dtype=torch.float32))
        self.register_buffer("y_std", torch.ones(2) if y_std is None else torch.as_tensor(y_std, dtype=torch.float32))
        from .physics import Physics  # 循環 import を避けるため遅延 import
        self._phys = Physics(self.phys_params, self.dt)

    def config(self):
        return {"kind": self.kind, "width": self.width, "depth": self.depth, "structured": False,
                "scale": self.scale, "dt": self.dt, "phys_params": self.phys_params}

    # ---- 解析モデルの予測と特徴量 ----
    @torch.no_grad()
    def physics_spec(self, q, dq, tau):
        """解析モデル（ストッパ無し）の 1 ステップ予測を spec 単位 [Δq/D_THETA_MAX, Δdq/D_DTHETA_MAX] で返す。"""
        q, dq, tau = q.double(), dq.double(), tau.double()
        q1, v1 = self._phys.step(q, dq, tau, stop=False)
        s = self.scale
        return torch.stack([(q1 - q) / s["D_THETA_MAX"], (v1 - dq) / s["D_DTHETA_MAX"]], dim=-1).float()

    def features(self, q, dq, tau):
        s = self.scale
        xs = torch.stack([torch.sin(q), torch.cos(q), dq / s["DTHETA_MAX"], tau / s["TAU_MAX"]], dim=-1).float()
        return torch.cat([xs, self.physics_spec(q, dq, tau)], dim=-1)

    def features_from_raw(self, X, bs=262144):
        """X = [q dq tau]（N, 3）から特徴量（N, 6）を作る（学習・評価の前処理。バッチごとに計算してメモリを抑える）。"""
        X = torch.as_tensor(X, dtype=torch.float64)
        return torch.cat([self.features(X[i:i + bs, 0], X[i:i + bs, 1], X[i:i + bs, 2]) for i in range(0, len(X), bs)])

    # ---- NN ----
    def forward(self, feat):
        """特徴量（6 列）→ 標準化した残差。"""
        return self.net((feat - self.x_mean) / self.x_std)

    def predict_spec(self, feat):
        """特徴量から、真の次状態の予測（spec 単位）= 解析モデルの予測 + NN の補正。"""
        return feat[..., 4:6] + self.forward(feat) * self.y_std

    def step_physical(self, q, dq, tau):
        s = self.scale
        y = self.predict_spec(self.features(q, dq, tau)).to(q.dtype)
        return q + y[..., 0] * s["D_THETA_MAX"], dq + y[..., 1] * s["D_DTHETA_MAX"]

    def delta_physical(self, xs, x_raw):
        """x_raw = [q dq tau]（物理量、(N, 3)）から、物理量の (Δq, Δdq) を返す。xs は使わない（共通インタフェース）。"""
        s = self.scale
        x_raw = torch.as_tensor(x_raw, dtype=torch.float64)
        y = self.predict_spec(self.features(x_raw[:, 0], x_raw[:, 1], x_raw[:, 2]))
        return y[:, 0].double() * s["D_THETA_MAX"], y[:, 1].double() * s["D_DTHETA_MAX"]


def build_model(config):
    """checkpoint の config から、重みを読み込む前のモデルを作る。"""
    if config.get("kind", "nsm") == "residual":
        return ResidualNSM(config["scale"], config["dt"], config["phys_params"], config["width"], config["depth"])
    return NSM(config["scale"], config["dt"], config["width"], config["depth"], config["structured"])
