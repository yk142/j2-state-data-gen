#!/usr/bin/env python3
"""Neural State Model を学習する。

    python python/train.py --data data/export/j2_flat_full.h5 --run baseline
    python python/train.py --data data/export/j2_flat_full.h5 --run structured --structured
    python python/train.py --data data/export/j2_flat_full.h5 --run residual --residual
    python python/train.py --run residual_robust --residual --y-scale free --loss huber
    python python/train.py --run residual_dual --residual --y-scale free --loss dual

--y-scale free : 残差の標準化に、リミット接触を除いた遷移の標準偏差を使う（既定 all は全遷移）。
                 接触の巨大な残差が標準偏差を支配し、接触なしの小さな残差が損失に効かなくなるのを避ける。
--loss huber   : Huber 損失（標準化後の残差が 1 以内は二乗、超える分は線形）。接触の外れ値に引きずられにくい。
                 ただし接触の巨大な残差（約 100σ）の勾配も抑えるため、動的な接触を学べない（#25 で判明）。
--loss dual    : 2 項の MSE。(1) 全遷移を「全遷移の標準偏差」で標準化した MSE（接触の巨大な残差を重視）
                 ＋ (2) 接触なしの遷移を「接触なしの標準偏差」で標準化した MSE（小さな残差を重視）。--y-scale free と併用。

--residual: 解析モデル（ストッパ無し）＋NN の残差型。NN の入力に解析モデルの予測を加え、
            目標は「真の次状態 − 解析モデルの予測」、最後の層を 0 初期化（開始時点は解析モデルそのもの）。

出力: data/train/<run>/{model.pt, history.json, config.json}（git 管理外）
"""
import argparse
import json
import math
import os
import sys
import time

import numpy as np
import torch

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from j2nsm.data import load  # noqa: E402
from j2nsm.model import NSM, ResidualNSM  # noqa: E402


def parse():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--data", default="data/export/j2_flat_full.h5")
    p.add_argument("--run", default="baseline", help="実行名（data/train/<run>/ に保存）")
    p.add_argument("--out-root", default="data/train")
    p.add_argument("--epochs", type=int, default=40)
    p.add_argument("--batch", type=int, default=4096)
    p.add_argument("--lr", type=float, default=2e-3)
    p.add_argument("--lr-min", type=float, default=1e-5)
    p.add_argument("--width", type=int, default=256)
    p.add_argument("--depth", type=int, default=3)
    p.add_argument("--structured", action="store_true", help="Δdq のみ予測し、Δq は運動学で求める")
    p.add_argument("--residual", action="store_true", help="解析モデル（ストッパ無し）＋NN の残差型")
    p.add_argument("--y-scale", choices=("all", "free"), default="all", help="残差型の標準化に使う標準偏差（all: 全遷移、free: 接触を除く）")
    p.add_argument("--loss", choices=("mse", "huber", "dual"), default="mse")
    p.add_argument("--dual-lambda", type=float, default=1.0, help="dual 損失の第 2 項（接触なし）の重み")
    p.add_argument("--patience", type=int, default=8, help="val 損失が改善しないエポック数の上限")
    p.add_argument("--max-train", type=int, default=0, help="学習行数の上限（動作確認用。0 は全件）")
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--threads", type=int, default=4)
    return p.parse_args()


def make_loss(kind, w=None, lam=1.0):
    """損失関数 (pred, y, free) → スカラー。free は「接触でない遷移」のマスク（dual のみ使う）。

    huber は標準化後の残差が 1 以内で二乗（0.5·e²）、超える分は線形。
    dual は mean_all(w·e²) + lam·mean_free(e²)。y は接触なしの標準偏差で標準化した残差で、
    w = (接触なしの標準偏差 / 全遷移の標準偏差)² が、第 1 項を「全遷移の標準偏差で標準化した MSE」にする。
    """
    def loss(pred, y, free=None):
        if kind == "huber":
            return torch.nn.functional.smooth_l1_loss(pred, y, beta=1.0)
        if kind == "dual":
            e2 = (pred - y) ** 2
            t1 = (e2 * w).mean()
            t2 = e2[free].mean() if bool(free.any()) else e2.new_zeros(())
            return t1 + lam * t2
        return torch.mean((pred - y) ** 2)
    return loss


@torch.no_grad()
def eval_loss(model, F, Yn, loss_fn, free=None, bs=32768):
    """特徴量 F と標準化済みの目標 Yn の損失（学習と同じ種類）。"""
    model.eval()
    tot, n = 0.0, 0
    for i in range(0, len(F), bs):
        k = Yn[i:i + bs].numel()
        tot += float(loss_fn(model(F[i:i + bs]), Yn[i:i + bs], None if free is None else free[i:i + bs])) * k
        n += k
    return tot / n


def main():
    a = parse()
    torch.manual_seed(a.seed)
    np.random.seed(a.seed)
    torch.set_num_threads(a.threads)
    out = os.path.join(a.out_root, a.run)
    os.makedirs(out, exist_ok=True)

    d = load(a.data, splits=("train", "val"))
    tr, va = d["splits"]["train"], d["splits"]["val"]
    if a.max_train and a.max_train < len(tr["X"]):
        idx = np.sort(np.random.permutation(len(tr["X"]))[:a.max_train])
        tr = {k: v[idx] for k, v in tr.items()}
    print(f"学習 {len(tr['X']):,} 行 / 検証 {len(va['X']):,} 行、残差型={a.residual}、構造付き={a.structured}、幅 {a.width} × 深さ {a.depth}")

    if a.residual:
        phys = {k: d["physics"][k] for k in ("M", "mgL", "Fc", "Bv", "eps", "qMin", "qMax")}
        model = ResidualNSM(d["scale"], d["dt"], phys, a.width, a.depth)
        t1 = time.time()
        Ftr, Fva = model.features_from_raw(tr["X"]), model.features_from_raw(va["X"])      # [X_spec, 解析モデルの予測]
        Rtr = torch.from_numpy(tr["Y_spec"]) - Ftr[:, 4:6]                                  # 残差 = 真値 − 解析モデル
        Rva = torch.from_numpy(va["Y_spec"]) - Fva[:, 4:6]
        print(f"解析モデルの予測を計算 {time.time() - t1:.0f}s。残差の標準偏差（spec 単位）{Rtr.std(0).tolist()}、"
              f"真の Y_spec の標準偏差 {torch.from_numpy(tr['Y_spec']).std(0).tolist()}")
        free = ~torch.from_numpy(tr["at_limit"])
        free_va = ~torch.from_numpy(va["at_limit"])
        if a.y_scale == "free":      # 接触を除いた遷移の標準偏差で標準化（接触の巨大な残差が支配しないように）
            ystd = Rtr[free].std(0)
            print(f"標準化に使う標準偏差（接触を除く）{ystd.tolist()}（全遷移 {Rtr.std(0).tolist()}）")
        else:
            ystd = Rtr.std(0)
        model.x_mean.copy_(Ftr.mean(0)); model.x_std.copy_(Ftr.std(0)); model.y_std.copy_(ystd)
        ytr_n, yva_n = Rtr / model.y_std, Rva / model.y_std
        w = (ystd / Rtr.std(0)) ** 2                     # dual 損失の第 1 項の重み（全遷移の標準偏差で標準化した MSE にする）
    else:
        free = free_va = w = None
        Ftr, Fva = torch.from_numpy(tr["X_spec"]), torch.from_numpy(va["X_spec"])
        model = NSM(d["scale"], d["dt"], a.width, a.depth, a.structured)
        tgt = model.target(torch.from_numpy(tr["Y_spec"]))
        model.x_mean.copy_(Ftr.mean(0)); model.x_std.copy_(Ftr.std(0))
        model.y_mean.copy_(tgt.mean(0)); model.y_std.copy_(tgt.std(0))
        ytr_n = (tgt - model.y_mean) / model.y_std
        yva_n = (model.target(torch.from_numpy(va["Y_spec"])) - model.y_mean) / model.y_std
    if a.loss == "dual" and not (a.residual and a.y_scale == "free"):
        raise SystemExit("--loss dual は --residual --y-scale free と併用してください")
    loss_fn = make_loss(a.loss, w, a.dual_lambda)
    n_param = sum(p.numel() for p in model.parameters())
    print(f"パラメータ数 {n_param:,}")
    if a.residual:   # 開始時点（NN の補正が 0 = 解析モデルそのもの）の損失
        print(f"開始時点（解析モデルそのもの）の val 損失 {eval_loss(model, Fva, yva_n, loss_fn, free_va):.6f}（損失={a.loss}）")

    opt = torch.optim.Adam(model.parameters(), lr=a.lr)
    steps_per_epoch = math.ceil(len(Ftr) / a.batch)
    total = a.epochs * steps_per_epoch
    warm = min(200, total // 10)

    def lr_at(s):
        if s < warm:
            return a.lr * (s + 1) / warm
        t = (s - warm) / max(1, total - warm)
        return a.lr_min + 0.5 * (a.lr - a.lr_min) * (1 + math.cos(math.pi * t))

    hist, best, bad, step = [], float("inf"), 0, 0
    t0 = time.time()
    for ep in range(1, a.epochs + 1):
        model.train()
        perm = torch.randperm(len(Ftr))
        run, cnt = 0.0, 0
        for i in range(0, len(Ftr), a.batch):
            idx = perm[i:i + a.batch]
            for g in opt.param_groups:
                g["lr"] = lr_at(step)
            loss = loss_fn(model(Ftr[idx]), ytr_n[idx], None if free is None else free[idx])
            opt.zero_grad(set_to_none=True)
            loss.backward()
            opt.step()
            run += float(loss.detach()) * len(idx); cnt += len(idx); step += 1
        vl = eval_loss(model, Fva, yva_n, loss_fn, free_va)
        rec = {"epoch": ep, "train_loss": run / cnt, "val_loss": vl, "lr": lr_at(step - 1), "time_s": time.time() - t0}
        hist.append(rec)
        flag = ""
        if vl < best:
            best, bad = vl, 0
            torch.save({"state_dict": model.state_dict(), "config": model.config(), "epoch": ep, "val_loss": vl},
                       os.path.join(out, "model.pt"))
            flag = " *"
        else:
            bad += 1
        print(f"epoch {ep:3d}  train {rec['train_loss']:.6f}  val {vl:.6f}  lr {rec['lr']:.2e}  {rec['time_s']:.0f}s{flag}", flush=True)
        with open(os.path.join(out, "history.json"), "w") as f:
            json.dump(hist, f, indent=1)
        if bad >= a.patience:
            print(f"早期終了（{a.patience} エポック改善なし）")
            break
    with open(os.path.join(out, "config.json"), "w") as f:
        json.dump({"args": vars(a), "n_param": n_param, "best_val_loss": best,
                   "train_time_s": time.time() - t0, "n_train": int(len(Ftr)), "n_val": int(len(Fva))}, f, indent=1)
    print(f"完了: best val {best:.6f}、学習時間 {time.time() - t0:.0f}s、保存 {out}")


if __name__ == "__main__":
    main()
