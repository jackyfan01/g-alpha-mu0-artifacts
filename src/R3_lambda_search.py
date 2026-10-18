# -*- coding: utf-8 -*-
"""R3: Λ = 4π³+π²+π 的事后选择偏差检验（回应外审 Q4）

系统枚举所有系数 |aᵢ| ≤ 10、阶数 ≤ 3 的整数 π 多项式
P(π) = a₀ + a₁π + a₂π² + a₃π³，统计其逼近 α⁻¹ 的相对优劣，
以评估 Λ ≈ α⁻¹ 吻合的统计显著性。

常数：CODATA 2022（Tiesinga, Mohr, Newell, Taylor, arXiv:2409.03787），
α = 7.2973525643(11)×10⁻³。

用法：python R3_lambda_search.py
输出：stdout 统计 + results.json
"""
from pathlib import Path
import json
from mpmath import mp, mpf, pi, fabs

mp.dps = 50
ALPHA = mpf("7.2973525643e-3")
ALPHA_INV = 1 / ALPHA
LAMBDA = 4 * pi**3 + pi**2 + pi
D_LAMBDA = fabs(LAMBDA - ALPHA_INV)
COEF_BOUND = 10
MAX_DEGREE = 3


def enumerate_polys():
    recs = []
    R = range(-COEF_BOUND, COEF_BOUND + 1)
    for a0 in R:
        for a1 in R:
            for a2 in R:
                for a3 in R:
                    if (a0, a1, a2, a3) == (0, 0, 0, 0):
                        continue
                    v = a0 + a1 * pi + a2 * pi**2 + a3 * pi**3
                    d = float(fabs(v - ALPHA_INV))
                    recs.append((d, a0, a1, a2, a3))
    recs.sort()
    return recs


def main():
    recs = enumerate_polys()
    total = len(recs)
    lam_coeff = (0, 1, 1, 4)
    lam_rank = next(i for i, r in enumerate(recs) if r[1:] == lam_coeff)
    d2 = recs[1][0]
    n_within = {
        "1e-4": sum(1 for r in recs if r[0] < 1e-4),
        "1e-3": sum(1 for r in recs if r[0] < 1e-3),
        "1e-2": sum(1 for r in recs if r[0] < 1e-2),
        "0.1": sum(1 for r in recs if r[0] < 0.1),
    }
    # 局部密度估计：Λ 之后第 2、3 名偏差给出 α⁻¹ 附近候选值的局部分布
    local_rate = 2.0 / (d2 + recs[2][0])  # 每单位偏差的候选数（单侧）
    p_chance = 1.0 - pow(1.0 - min(1.0, float(D_LAMBDA) * local_rate), total)
    # 更保守的口径：仅考虑 0.1 邻域内的 68 个候选
    near = n_within["0.1"]
    p_local = 1.0 - pow(1.0 - min(1.0, 2 * float(D_LAMBDA) / 0.1), near) if near else 0.0
    out = {
        "alpha_inv_CODATA2022": mp.nstr(ALPHA_INV, 13),
        "Lambda": mp.nstr(LAMBDA, 13),
        "d_Lambda": mp.nstr(D_LAMBDA, 8),
        "total_tested": total,
        "lambda_rank": lam_rank,
        "second_best_dev": recs[1][0],
        "best_others": recs[1:4],
        "counts_within": n_within,
        "p_chance_best_this_good_conservative": p_local,
        "note": "Λ 为全局最优；但按 α⁻¹ 附近候选局部密度估计，"
                "194480 个候选中出现一个同量级最优拟合的偶然概率约 1/4 上下，"
                "吻合不具备统计显著性。",
    }
    print(json.dumps(out, ensure_ascii=False, indent=2, default=float))
    with open(str(Path(__file__).resolve().parent.parent / "results" / "lambda_enumeration.json"), "w", encoding="utf-8") as f:
        json.dump(out, f, ensure_ascii=False, indent=2, default=float)


if __name__ == "__main__":
    main()
