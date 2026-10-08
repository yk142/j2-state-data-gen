#!/usr/bin/env python3
"""複数の run の評価結果（metrics.json）を 1 つの表に並べる（Markdown）。

    python python/compare.py baseline=reports/issue-25/baseline residual=reports/issue-25/residual --out docs/compare.md

列は run 名（学習モデル）と、基準線（解析モデル、状態を変えない）。最初の run の基準線を使う。
"""
import argparse
import json
import os
import sys

import numpy as np


def load(path):
    p = path if path.endswith(".json") else os.path.join(path, "metrics.json")
    with open(p) as f:
        return json.load(f)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("runs", nargs="+", help="ラベル=metrics.json（または評価ディレクトリ）")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    runs = [(r.split("=", 1)[0], load(r.split("=", 1)[1])) for r in a.runs]
    base = runs[0][1]
    head = [lab for lab, _ in runs] + ["解析（ストッパ無し）", "解析（非弾性ストッパ）", "状態を変えない"]
    out = []

    def table(title, rows, fmt="{:.2f}"):
        out.append(f"\n### {title}\n")
        out.append("| 区分 | " + " | ".join(head) + " |")
        out.append("|---|" + "---|" * len(head))
        for name, vals in rows:
            cells = [("—" if v is None or (isinstance(v, float) and np.isnan(v)) else fmt.format(v)) for v in vals]
            out.append(f"| {name} | " + " | ".join(cells) + " |")

    # ---- 1 ステップ（加速度換算 RMSE）----
    rows = []
    for split, key, label in (("test", "all", "test 全体"), ("test", "free(接触なし)", "test 接触なし"),
                              ("test", "contact(接触あり)", "test 接触あり"),
                              ("test", "contact:dynamic(動的な接触)", "　└ 動的な接触"),
                              ("test", "contact:static(静止押し付け)", "　└ 静止押し付け"),
                              ("benchmark", "all", "PTP ベンチマーク")):
        def get(m, meth):
            r = m["one_step"][split].get(key)
            return None if r is None else r[meth]["rmse_ddq_rad_s2"]
        n = base["one_step"][split].get(key, {}).get("n")
        rows.append((f"{label}（n={n:,}）" if n else label,
                     [get(m, "model") for _, m in runs] + [get(base, "analytic_free"), get(base, "analytic"), get(base, "zero")]))
    table("1 ステップ誤差（加速度換算 RMSE [rad/s²]、小さいほど良い）", rows)

    # ---- 1 秒窓の開ループ ----
    rows = []
    for sel, label in (("all", "全体"), ("free", "接触なしの窓"), ("contact", "接触を含む窓")):
        for h in ("100", "1000"):
            def getw(m, meth):
                c = m["rollout_windows"]["test"]["curves"].get(sel)
                return None if c is None else c[meth]["rmse_q_deg"].get(h)
            rows.append((f"{label}・{h} ステップ",
                         [getw(m, "model") for _, m in runs] + [None, getw(base, "analytic"), getw(base, "persistence")]))
    table("1 秒窓の開ループ（test、θ の RMSE [deg]）", rows, "{:.3f}")

    # ---- PTP の閉ループ ----
    def agg(m, meth, key, f):
        v = [r[meth][key] for r in m["closed_loop_benchmark"]]
        return f(v)
    rows = [("真値との θ の RMSE（12 本の平均）[deg]",
             [agg(m, "model", "vs_truth_rmse_q_deg", np.mean) for _, m in runs] + [None, agg(base, "analytic", "vs_truth_rmse_q_deg", np.mean), None]),
            ("最大の θ 誤差（12 本の最大）[deg]",
             [agg(m, "model", "vs_truth_max_q_deg", np.max) for _, m in runs] + [None, agg(base, "analytic", "vs_truth_max_q_deg", np.max), None]),
            ("参照への追従誤差 RMSE（平均）[deg]",
             [agg(m, "model", "track_rmse_deg", np.mean) for _, m in runs] + [None, agg(base, "analytic", "track_rmse_deg", np.mean), None]),
            ("トルクの RMSE vs 真値（平均）[N·m]",
             [agg(m, "model", "tau_rmse_N_m", np.mean) for _, m in runs] + [None, agg(base, "analytic", "tau_rmse_N_m", np.mean), None]),
            ("発散した本数（/ 12）",
             [sum(r["model"]["diverged"] for r in m["closed_loop_benchmark"]) for _, m in runs] + [None, sum(r["analytic"]["diverged"] for r in base["closed_loop_benchmark"]), None])]
    table("PTP ベンチマークの閉ループ", rows, "{:.4g}")
    truth = np.mean([r["truth_track_rmse_deg"] for r in base["closed_loop_benchmark"]])
    out.append(f"\n（参照への追従誤差の真値（記録）は平均 {truth:.4f}°）")
    text = "\n".join(out) + "\n"
    print(text)
    if a.out:
        with open(a.out, "w", encoding="utf-8") as f:
            f.write(text)


if __name__ == "__main__":
    sys.exit(main())
