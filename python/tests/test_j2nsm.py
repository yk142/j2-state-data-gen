"""j2nsm の単体テスト（unittest）。実行: python -m unittest discover -s python/tests -v"""
import os
import sys
import tempfile
import unittest

import numpy as np
import torch

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from helpers import PHYS, SCALE, write_tiny_h5  # noqa: E402
from j2nsm import data as D  # noqa: E402
from j2nsm import evaluate_lib as ev  # noqa: E402
from j2nsm.model import NSM  # noqa: E402
from j2nsm.physics import Physics  # noqa: E402

DT = 1e-3
T = lambda x: torch.tensor(x, dtype=torch.float64)  # noqa: E731


class TestPhysics(unittest.TestCase):
    def setUp(self):
        self.ph = Physics(PHYS, DT)

    def test_upright_is_equilibrium_and_unstable(self):
        q, dq = self.ph.step(T([0.0]), T([0.0]), T([0.0]))
        self.assertAlmostEqual(float(q), 0.0, places=12)
        # θ=+0.5 rad では重力が θ を増やす向き（摩擦 Fc=23.8 < mgL·sin 0.5=35 で動き出す）
        q1, dq1 = self.ph.step(T([0.5]), T([0.0]), T([0.0]))
        self.assertGreater(float(dq1), 0.0)

    def test_gravity_hold_torque(self):
        # 保持トルク −mgL·sin θ で静止（θ = −π/2 では +mgL）
        q, dq = T([-np.pi / 2]), T([0.0])
        for _ in range(100):
            q, dq = self.ph.step(q, dq, T([PHYS["mgL"]]))
        self.assertLess(abs(float(q) + np.pi / 2), 1e-9)

    def test_friction_decelerates_and_stops(self):
        q, dq = T([0.0]), T([0.2])
        for _ in range(60):
            q, dq = self.ph.step(q, dq, T([0.0]))
        self.assertLess(abs(float(dq)), 1e-2)

    def test_stop_keeps_position_inside_limits(self):
        q, dq = T([PHYS["qMax"] - 0.01]), T([3.0])
        for _ in range(50):
            q, dq = self.ph.step(q, dq, T([500.0]))
        self.assertLessEqual(float(q), PHYS["qMax"] + 1e-12)
        self.assertLessEqual(float(dq), 1e-12)               # 外向きの速度は 0 に
        q, dq = self.ph.step(T([0.0]), T([1.0]), T([0.0]), stop=False)
        self.assertTrue(np.isfinite(float(q)))


class TestModel(unittest.TestCase):
    def test_shapes_and_buffers(self):
        for structured in (False, True):
            m = NSM(SCALE, DT, width=16, depth=2, structured=structured)
            xs = torch.randn(7, 4)
            self.assertEqual(m(xs).shape, (7, 1 if structured else 2))
            self.assertEqual(m.target(torch.randn(7, 2)).shape, (7, 1 if structured else 2))
            for k in ("x_mean", "x_std", "y_mean", "y_std"):
                self.assertIn(k, m.state_dict())             # 標準化の定数は保存される

    def test_structured_delta_q_is_kinematic(self):
        m = NSM(SCALE, DT, width=16, depth=2, structured=True)
        q, dq, tau = T([-1.0, 0.3]), T([1.5, -2.0]), T([10.0, -50.0])
        xs = torch.stack([torch.sin(q), torch.cos(q), dq / SCALE["DTHETA_MAX"], tau / SCALE["TAU_MAX"]], -1).float()
        d_q, d_dq = m.delta_physical(xs)
        np.testing.assert_allclose(d_q.detach().numpy(), (DT * (dq.numpy() + 0.5 * d_dq.detach().numpy())), rtol=1e-5, atol=1e-8)

    def test_step_physical_matches_delta_physical(self):
        for structured in (False, True):
            m = NSM(SCALE, DT, width=16, depth=2, structured=structured)
            q, dq, tau = T([-1.0, 0.3, 1.0]), T([1.5, -2.0, 0.0]), T([10.0, -50.0, 3.0])
            xs = torch.stack([torch.sin(q), torch.cos(q), dq / SCALE["DTHETA_MAX"], tau / SCALE["TAU_MAX"]], -1).float()
            d_q, d_dq = m.delta_physical(xs)
            q1, dq1 = m.step_physical(q, dq, tau)
            np.testing.assert_allclose((q1 - q).detach().numpy(), d_q.detach().numpy(), rtol=1e-4, atol=1e-7)
            np.testing.assert_allclose((dq1 - dq).detach().numpy(), d_dq.detach().numpy(), rtol=1e-4, atol=1e-7)


class TestData(unittest.TestCase):
    def test_load_and_sequences_roundtrip(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "t.h5")
            seqs = write_tiny_h5(path)
            d = D.load(path)
            self.assertAlmostEqual(d["dt"], 1e-3)
            self.assertAlmostEqual(d["scale"]["TAU_MAX"], 860.4)             # 長さ 1 の配列からスカラーへ
            self.assertAlmostEqual(d["physics"]["mgL"], 73.25)
            self.assertEqual(d["scenarios"]["pattern"][1], "contact_fall")
            tr = D.sequences(d, "train")
            self.assertEqual([s["name"] for s in tr], ["s0", "s1"])
            for s, ref in zip(tr, [seqs[0], seqs[1]]):
                np.testing.assert_allclose(s["q"], ref[1], atol=1e-5)         # float32 の丸め
                np.testing.assert_allclose(s["dq"], ref[2], atol=1e-5)
                np.testing.assert_allclose(s["tau"], ref[3], atol=1e-4)
            self.assertTrue(tr[1]["at_limit"][5:15].all() and not tr[1]["at_limit"][:5].any())
            self.assertEqual(tr[0]["group"], "bln"); self.assertEqual(tr[1]["group"], "contact")
            b = D.sequences(d, "benchmark")[0]
            self.assertFalse(np.isnan(b["qref"]).any())
            self.assertTrue(np.isnan(tr[0]["qref"]).all())

    def test_sequences_refuse_filtered_export(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "t.h5")
            write_tiny_h5(path, exclude_at_limit=True)
            with self.assertRaises(ValueError):
                D.sequences(D.load(path), "train")

    def test_group_of(self):
        self.assertEqual([D.group_of(p) for p in ("bln_high", "hold_+030", "contact_bln", "ptp_v20", "freefall", "step")],
                         ["bln", "hold", "contact", "ptp", "freefall", "step"])


class TestRollout(unittest.TestCase):
    def test_open_loop_rollout_with_exact_stepper_reproduces_trajectory(self):
        ph = Physics(PHYS, DT)
        B, H = 3, 50
        q0, dq0 = T([-1.0, -0.5, 0.2]), T([0.0, 0.5, -0.3])
        tau = T(np.tile(np.linspace(-20, 20, H), (B, 1)))
        q, dq = ev.rollout(lambda a, b, c: ph.step(a, b, c), q0, dq0, tau)
        self.assertEqual(q.shape, (B, H + 1))
        # 1 ステップずつ手で進めた結果と一致
        qq, vv = q0.clone(), dq0.clone()
        for k in range(H):
            qq, vv = ph.step(qq, vv, tau[:, k])
        np.testing.assert_allclose(q[:, -1].numpy(), qq.numpy(), atol=1e-12)

    def test_windows_cover_expected_starts(self):
        seqs = [{"q": np.arange(101.0), "dq": np.arange(101.0), "tau": np.arange(100.0),
                 "at_limit": np.r_[np.zeros(60, bool), np.ones(40, bool)], "group": "bln", "name": "a"}]
        W = ev.make_windows(seqs, horizon=40, stride=20)
        self.assertEqual(len(W["q"]), 4)                                       # 開始 0, 20, 40, 60
        self.assertEqual(W["q"].shape, (4, 41)); self.assertEqual(W["tau"].shape, (4, 40))
        self.assertEqual(W["contact"].tolist(), [False, False, True, True])    # 開始 20 は区間 20..59 で接触なし
        self.assertIsNone(ev.make_windows(seqs, horizon=200, stride=20))

    def test_closed_loop_with_exact_plant_tracks_smooth_reference(self):
        """解析モデル自身をプラントにして PD＋FF を閉ループにすると、滑らかな参照に高精度で追従する（制御器の実装の確認）。"""
        ph = Physics(PHYS, DT)
        t = np.arange(0, 3.0, DT)
        qref = -1.0 + 0.4 * np.sin(2 * np.pi * 0.5 * t)
        seq = {"q": np.r_[qref, qref[-1]], "dq": np.zeros(len(t) + 1), "qref": qref}
        seq["dq"][0] = 0.4 * 2 * np.pi * 0.5
        Q, V, U, qr = ev.closed_loop_rollout(lambda a, b, c: ph.step(a, b, c, stop=False), seq, PHYS, DT, PHYS["tauPeak"])
        self.assertTrue(np.isfinite(Q).all())
        self.assertLess(np.degrees(np.max(np.abs(Q[:-1] - qref))), 0.2)         # 追従誤差 0.2° 未満
        self.assertLess(np.max(np.abs(U)), PHYS["tauPeak"] + 1e-9)

    def test_closed_loop_flags_divergence(self):
        seq = {"q": np.zeros(201), "dq": np.zeros(201), "qref": np.zeros(200)}
        bad = lambda q, dq, tau: (q + 20.0, dq)  # noqa: E731                    # 即座に発散するプラント
        Q, V, U, _ = ev.closed_loop_rollout(bad, seq, PHYS, DT, PHYS["tauPeak"])
        self.assertTrue(np.isnan(Q[-1]))


if __name__ == "__main__":
    unittest.main()
