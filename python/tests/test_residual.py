"""ResidualNSM（解析モデル＋NN の残差型）の単体テスト。実行: python -m unittest discover -s python/tests -v"""
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
from j2nsm.model import NSM, ResidualNSM, build_model  # noqa: E402
from j2nsm.physics import Physics  # noqa: E402
from train import make_loss  # noqa: E402

DT = 1e-3
PP = {k: PHYS[k] for k in ("M", "mgL", "Fc", "Bv", "eps", "qMin", "qMax")}
T = lambda x: torch.tensor(x, dtype=torch.float64)  # noqa: E731


def make(width=16, depth=2):
    return ResidualNSM(SCALE, DT, PP, width, depth)


class TestResidualNSM(unittest.TestCase):
    def test_zero_init_equals_physics(self):
        """最後の層が 0 初期化なので、学習前は解析モデル（ストッパ無し）の 1 ステップ予測そのもの。"""
        m = make()
        q, dq, tau = T([-1.0, 0.3, 1.0, -2.0]), T([1.5, -2.0, 0.0, 0.5]), T([10.0, -50.0, 3.0, 300.0])
        q1, dq1 = m.step_physical(q, dq, tau)
        qe, de = Physics(PP, DT).step(q, dq, tau, stop=False)
        np.testing.assert_allclose(q1.detach().numpy(), qe.numpy(), atol=2e-6)     # float32 の丸め
        np.testing.assert_allclose(dq1.detach().numpy(), de.numpy(), atol=2e-6)

    def test_features_shape_and_physics_columns(self):
        m = make()
        X = np.array([[-1.0, 1.5, 10.0], [0.3, -2.0, -50.0]])
        F = m.features_from_raw(X)
        self.assertEqual(F.shape, (2, 6))
        p = m.physics_spec(T(X[:, 0]), T(X[:, 1]), T(X[:, 2]))
        np.testing.assert_allclose(F[:, 4:6].numpy(), p.numpy(), atol=1e-7)
        # X_spec の部分は MATLAB のエクスポートと同じ定義
        np.testing.assert_allclose(F[:, 2].numpy(), X[:, 1] / SCALE["DTHETA_MAX"], atol=1e-6)

    def test_delta_physical_matches_step_physical(self):
        m = make()
        torch.manual_seed(0)
        for p in m.parameters():                      # 補正が 0 でない状態でも一致すること
            torch.nn.init.normal_(p, std=0.05)
        X = np.array([[-1.0, 1.5, 10.0], [0.3, -2.0, -50.0], [1.0, 0.0, 3.0]])
        dqq, dvv = m.delta_physical(None, X)
        q1, v1 = m.step_physical(T(X[:, 0]), T(X[:, 1]), T(X[:, 2]))
        np.testing.assert_allclose(dqq.detach().numpy(), (q1 - T(X[:, 0])).detach().numpy(), atol=1e-7)
        np.testing.assert_allclose(dvv.detach().numpy(), (v1 - T(X[:, 1])).detach().numpy(), atol=1e-7)

    def test_residual_target_plus_physics_equals_truth(self):
        """学習の目標（真の Y_spec − 解析モデルの予測）に、解析モデルの予測を足すと真の Y_spec に戻る。"""
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "t.h5")
            write_tiny_h5(path)
            d = D.load(path)
            tr = d["splits"]["train"]
            m = make()
            F = m.features_from_raw(tr["X"])
            R = torch.from_numpy(tr["Y_spec"]) - F[:, 4:6]
            np.testing.assert_allclose((R + F[:, 4:6]).numpy(), tr["Y_spec"], atol=1e-7)

    def test_buffers_in_state_dict_and_build_model_roundtrip(self):
        m = make(width=8, depth=1)
        with torch.no_grad():
            m.y_std.copy_(torch.tensor([0.5, 2.0]))
            m.x_std.copy_(torch.arange(1.0, 7.0))
        cfg = m.config()
        self.assertEqual(cfg["kind"], "residual")
        m2 = build_model(cfg)
        self.assertIsInstance(m2, ResidualNSM)
        m2.load_state_dict(m.state_dict())
        np.testing.assert_allclose(m2.y_std.numpy(), [0.5, 2.0])
        np.testing.assert_allclose(m2.x_std.numpy(), np.arange(1.0, 7.0))
        self.assertIsInstance(build_model(NSM(SCALE, DT, 8, 1).config()), NSM)    # 純粋な NN も復元できる

    def test_one_step_predictions_work_for_both_kinds(self):
        """評価コードが、純粋な NN・残差型のどちらでも動き、残差型の開始時点は解析モデル（ストッパ無し）と一致する。"""
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "t.h5")
            write_tiny_h5(path)
            d = D.load(path)
            phys = Physics(PP, DT)
            for model in (NSM(SCALE, DT, 8, 1), make()):
                preds = ev.one_step_predictions(model, phys, d, "test")
                self.assertEqual(preds["model"].shape, (60, 2))
            r = ev.one_step_predictions(make(), phys, d, "test")
            np.testing.assert_allclose(r["model"], r["analytic_free"], atol=3e-6)

    def test_contact_breakdown_in_one_step_table(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "t.h5")
            write_tiny_h5(path)
            d = D.load(path)
            # train の contact_fall シナリオ（at_limit が 10 遷移）。動的／静止の内訳が全体の和になること
            tab, _ = ev.one_step_table(make(), Physics(PP, DT), d, "train")
            self.assertIn("contact:dynamic(動的な接触)", tab)
            n = tab["contact(接触あり)"]["n"]
            nd = tab["contact:dynamic(動的な接触)"]["n"]
            ns = tab.get("contact:static(静止押し付け)", {"n": 0})["n"]
            self.assertEqual(n, 10)
            self.assertEqual(nd + ns, n)


class TestLosses(unittest.TestCase):
    def test_mse_and_huber(self):
        pred, y = torch.tensor([[0.0, 3.0]]), torch.tensor([[0.5, 0.0]])
        self.assertAlmostEqual(float(make_loss("mse")(pred, y)), (0.25 + 9) / 2, places=6)
        # Huber（beta=1）: |e|<1 は 0.5e²、|e|>=1 は |e|−0.5 → (0.125 + 2.5)/2
        self.assertAlmostEqual(float(make_loss("huber")(pred, y)), (0.125 + 2.5) / 2, places=6)

    def test_huber_is_less_sensitive_to_outliers_than_mse(self):
        y = torch.zeros(100, 2)
        pred = torch.zeros(100, 2)
        pred[0] = 100.0                                                          # 外れ値 1 行
        self.assertGreater(float(make_loss("mse")(pred, y)), 100 * float(make_loss("huber")(pred, y)))

    def test_dual_combines_all_rows_and_free_rows(self):
        w = torch.tensor([0.01, 0.04])
        pred = torch.tensor([[1.0, 1.0], [2.0, 2.0], [0.0, 0.0], [0.0, 0.0]])
        y = torch.zeros(4, 2)
        free = torch.tensor([False, False, True, True])                          # 先頭 2 行は接触
        loss = make_loss("dual", w, 1.0)(pred, y, free)
        t1 = float(((pred - y) ** 2 * w).mean())
        self.assertAlmostEqual(float(loss), t1 + 0.0, places=6)                  # 接触なしの行は誤差 0 → 第 2 項は 0
        pred2 = pred.clone(); pred2[2:] = 1.0                                    # 接触なしの行に誤差を入れると第 2 項が効く
        self.assertAlmostEqual(float(make_loss("dual", w, 1.0)(pred2, y, free)),
                               float(((pred2 - y) ** 2 * w).mean()) + 1.0, places=5)
        self.assertAlmostEqual(float(make_loss("dual", w, 0.0)(pred2, y, free)), float(((pred2 - y) ** 2 * w).mean()), places=6)

    def test_dual_handles_batch_without_free_rows(self):
        pred, y = torch.ones(3, 2), torch.zeros(3, 2)
        loss = make_loss("dual", torch.ones(2), 1.0)(pred, y, torch.zeros(3, dtype=torch.bool))
        self.assertTrue(torch.isfinite(loss))


if __name__ == "__main__":
    unittest.main()
