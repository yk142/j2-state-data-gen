#!/usr/bin/env python3
"""Neural State Model を学習する。

    python python/train.py --data data/export/j2_flat_full.h5 --run baseline
    python python/train.py --data data/export/j2_flat_full.h5 --run structured --structured

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
from j2nsm.model import NSM  # noqa: E402


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
    p.add_argument("--patience", type=int, default=8, help="val 損失が改善しないエポック数の上限")
    p.add_argument("--max-train", type=int, default=0, help="学習行数の上限（動作確認用。0 は全件）")
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--threads", type=int, default=4)
    return p.parse_args()


@torch.no_grad()
def eval_loss(model, X, Y, bs=32768):
    model.eval()
    tot, n = 0.0, 0
    for i in range(0, len(X), bs):
        pred = model(X[i:i + bs])
        y = (model.target(Y[i:i + bs]) - model.y_mean) / model.y_std
        tot += float(((pred - y) ** 2).sum())
        n += y.numel()
    return tot / n


def main():
    a = parse()
    torch.manual_seed(a.seed)
    np.random.seed(a.seed)
    torch.set_num_threads(a.threads)
    out = os.path.join(a.out_root, a.run)
    os.makedirs(out, exist_ok=True)

    d = load(a.data, splits=("train", "val"))
    Xtr = torch.from_numpy(d["splits"]["train"]["X_spec"])
    Ytr = torch.from_numpy(d["splits"]["train"]["Y_spec"])
    Xva = torch.from_numpy(d["splits"]["val"]["X_spec"])
    Yva = torch.from_numpy(d["splits"]["val"]["Y_spec"])
    if a.max_train and a.max_train < len(Xtr):
        idx = torch.randperm(len(Xtr))[:a.max_train]
        Xtr, Ytr = Xtr[idx], Ytr[idx]
    print(f"学習 {len(Xtr):,} 行 / 検証 {len(Xva):,} 行、構造付き={a.structured}、幅 {a.width} × 深さ {a.depth}")

    model = NSM(d["scale"], d["dt"], a.width, a.depth, a.structured)
    tgt = model.target(Ytr)
    model.x_mean.copy_(Xtr.mean(0)); model.x_std.copy_(Xtr.std(0))
    model.y_mean.copy_(tgt.mean(0)); model.y_std.copy_(tgt.std(0))
    n_param = sum(p.numel() for p in model.parameters())
    print(f"パラメータ数 {n_param:,}")

    opt = torch.optim.Adam(model.parameters(), lr=a.lr)
    steps_per_epoch = math.ceil(len(Xtr) / a.batch)
    total = a.epochs * steps_per_epoch
    warm = min(200, total // 10)

    def lr_at(s):
        if s < warm:
            return a.lr * (s + 1) / warm
        t = (s - warm) / max(1, total - warm)
        return a.lr_min + 0.5 * (a.lr - a.lr_min) * (1 + math.cos(math.pi * t))

    ytr_n = (tgt - model.y_mean) / model.y_std
    hist, best, bad, step = [], float("inf"), 0, 0
    t0 = time.time()
    for ep in range(1, a.epochs + 1):
        model.train()
        perm = torch.randperm(len(Xtr))
        run, cnt = 0.0, 0
        for i in range(0, len(Xtr), a.batch):
            idx = perm[i:i + a.batch]
            for g in opt.param_groups:
                g["lr"] = lr_at(step)
            loss = torch.mean((model(Xtr[idx]) - ytr_n[idx]) ** 2)
            opt.zero_grad(set_to_none=True)
            loss.backward()
            opt.step()
            run += float(loss.detach()) * len(idx); cnt += len(idx); step += 1
        vl = eval_loss(model, Xva, Yva)
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
                   "train_time_s": time.time() - t0, "n_train": int(len(Xtr)), "n_val": int(len(Xva))}, f, indent=1)
    print(f"完了: best val {best:.6f}、学習時間 {time.time() - t0:.0f}s、保存 {out}")


if __name__ == "__main__":
    main()
