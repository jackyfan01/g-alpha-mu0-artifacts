# -*- coding: utf-8 -*-
r"""审稿人 Q2 + Q4 的补强计算。

Q2 完全相关证书：在同一个 Arb 算术 trace 内由 α 导出 μ₀ = 2αh/(e²c)（h、e、c 为 SI 精确值），
   于是 α²μ₀ = 2α³h/(e²c)，不确定度只来自 α（相对 3·u_r(α)）。与"α、μ₀ 各给独立球"的版本对照。
Q4 间距感知零模型：把"随机靶均匀落在家族值域上"的 p 值改成**解析精确**的最近邻覆盖测度
   p = |∪_i [V_i-d, V_i+d] ∩ [V_min, V_max]| / (V_max - V_min)，
   它天然按 Voronoi 胞宽加权（大胞被抽中的概率与宽度成正比），不需要 Monte Carlo。
   同时给出缩族/扩族各口径的 p 与跨家族的最小 p（供 multiplicity 讨论）。
"""
import itertools
import json
import math

import numpy as np
from flint import arb, ctx

ctx.prec = 128
OUT = {}

# ---------------- Q2: 完全相关 vs 独立球 ----------------
a_lo_rel = "7.2973525643e-3"
h = arb("6.62607015e-34")
e = arb("1.602176634e-19")
c = arb("299792458")
alpha = arb(a_lo_rel, "1.1e-12")
G = arb("6.67430e-11", "1.5e-15")

mu0_from_alpha = 2 * alpha * h / (e * e * c)          # 完全相关
a2mu0_corr = alpha * alpha * mu0_from_alpha
mu0_indep = arb("1.25663706127e-6", "2.0e-16")       # 独立球版本
a2mu0_indep = alpha * alpha * mu0_indep


def ball(x):
    return [float(x.lower()), float(x.upper()), float(x.rad() / x.mid())]


for name, b in [("corr", a2mu0_corr), ("indep", a2mu0_indep)]:
    lo, hi, ur = ball(b)
    gap = lo - float(G.upper())
    dev = (float(b.mid()) - float(G.mid())) / float(G.mid())
    OUT[name] = {"interval": [lo, hi], "u_r": ur,
                 "mu0": [float(mu0_from_alpha.lower()), float(mu0_from_alpha.upper())] if name == "corr"
                 else [1.25663706127e-6 - 2.0e-16, 1.25663706127e-6 + 2.0e-16],
                 "disjoint": bool(not b.overlaps(G)),
                 "gap": gap, "gap_in_uG": gap / float(G.rad()),
                 "dev_percent": dev * 100.0, "sigma": abs(dev) * float(G.mid()) / float(G.rad())}
OUT["mu0_from_alpha_vs_CODATA"] = {
    # 注意：python-flint 对极窄球的 str() 会输出 "[lo ..." 形式，直接切片会得到残缺字符串，
    # 因此这里落 float（16 位有效数字，已远超本文所需精度）并另存端点。
    "derived": float(mu0_from_alpha.mid()),
    "derived_interval": [float(mu0_from_alpha.lower()), float(mu0_from_alpha.upper())],
    "codata_displayed": "1.25663706127e-06",
    "rel_difference": float(abs(mu0_from_alpha.mid() - arb("1.25663706127e-6")) / mu0_from_alpha.mid())}

# ---------------- Q4: 间距感知精确 p ----------------
def values_pi_poly(B, deg=3):
    v = []
    for coeffs in itertools.product(range(-B, B + 1), repeat=deg + 1):
        if not any(coeffs):
            continue
        v.append(sum(a * (math.pi ** k) for k, a in enumerate(coeffs)))
    return np.array(sorted(v))


def spacing_p(V, d):
    """|∪[V_i-d,V_i+d] ∩ [V0,Vn]| / (Vn-V0)，区间合并后精确计算。"""
    lo, hi = float(V[0]), float(V[-1])
    R = hi - lo
    iv = np.stack([np.maximum(V - d, lo), np.minimum(V + d, hi)], axis=1)
    iv = iv[iv[:, 1] > iv[:, 0]]
    total, cur_l, cur_r = 0.0, None, None
    for l, r in iv:
        if cur_l is None:
            cur_l, cur_r = l, r
        elif l <= cur_r:
            cur_r = max(cur_r, r)
        else:
            total += cur_r - cur_l
            cur_l, cur_r = l, r
    if cur_l is not None:
        total += cur_r - cur_l
    return total / R


ALPHA_INV = 1 / 7.2973525643e-3
LAM = 4 * math.pi ** 3 + math.pi ** 2 + math.pi
d_obs = LAM - ALPHA_INV

fams = {}
for B in (4, 5, 6, 8, 10):
    V = values_pi_poly(B)
    R = float(V[-1] - V[0])
    N = len(V)
    p_sp = spacing_p(V, d_obs)
    fam_bonf = N * 2 * d_obs / R
    fam_analytic = 1 - (1 - 2 * d_obs / R) ** N
    # 该家族内的最优距离（用于检查 Λ 是否仍唯一最优）
    best = float(np.min(np.abs(V - ALPHA_INV)))
    fams["pi|a|<=%d" % B] = {"N": N, "R": R, "d_obs": d_obs, "best_in_family": best,
                              "p_spacing_exact": p_sp, "p_bonferroni": fam_bonf,
                              "p_analytic_indep": fam_analytic}
OUT["families"] = fams
OUT["d_obs"] = d_obs
OUT["Lambda"] = LAM
OUT["alpha_inv"] = ALPHA_INV

# 跨家族 multiplicity：在 K 个家族口径下都做了检验，最小 p 与 Bonferroni-over-families
ps = [v["p_spacing_exact"] for v in fams.values()]
K = len(ps)
OUT["cross_family"] = {"K": K, "min_p": min(ps), "bonferroni_over_families": min(1.0, min(ps) * K),
                      "note": "跨家族最小 p 是选择后再看的，不能当频率学 p 值用；给出以说明敏感性"}

# 混合基 (pi,e) 家族：6 系数、|a|<=5，1,771,560 个值
pi_e = []
vals_pi = [1.0]
pw = 1.0
for k in range(1, 7):
    pw *= math.pi
    vals_pi.append(pw)
vals_e = [1.0]
pw = 1.0
for k in range(1, 7):
    pw *= math.e
    vals_e.append(pw)
# 基 = {pi^1..pi^6, e^1..e^6}，每系数 |a|<=5：太大，改为脚本 q3b 实际用的 6 系数结构
OUT["mixed_family_note"] = "混合基口径见 q3b_family_extension.py 的独立结果，本脚本只重算 pi 多项式家族"

import os as _os
_out = _os.path.normpath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)),
                                      "..", "results", "q5_correlated_and_spacing.json"))
json.dump(OUT, open(_out, "w", encoding="utf-8"), indent=1, ensure_ascii=False, default=float)
print("写出:", _out)

print("=== Q2 完全相关 vs 独立球 ===")
for k in ("corr", "indep"):
    d = OUT[k]
    print("%-6s α²μ₀∈[%.10e, %.10e] u_r=%.4e disjoint=%s 间距=%.6e = %.2f u(G) 偏差=%.6f%% σ=%.1f"
          % (k, d["interval"][0], d["interval"][1], d["u_r"], d["disjoint"], d["gap"],
             d["gap_in_uG"], d["dev_percent"], d["sigma"]))
print("由 α 导出的 μ₀ =", OUT["mu0_from_alpha_vs_CODATA"]["derived"],
      "与 CODATA 展示值相对差 = %.3e" % OUT["mu0_from_alpha_vs_CODATA"]["rel_difference"])
print("\n=== Q4 间距感知精确 p（d_obs = %.6e）===" % d_obs)
for name, v in fams.items():
    print("%-12s N=%8d R=%8.2f p_精确=%.5f p_Bonf=%.5f p_解析=%.5f 家族内最优=%.3e"
          % (name, v["N"], v["R"], v["p_spacing_exact"], v["p_bonferroni"],
             v["p_analytic_indep"], v["best_in_family"]))
print("跨家族:", OUT["cross_family"])
