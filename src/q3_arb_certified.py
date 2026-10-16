# -*- coding: utf-8 -*-
"""Q3: 用 Arb（python-flint）定向舍入球算术重跑区间验证，替代 mpmath.mpi 的传播包络。

与 mpmath.mpi 的区别（外审 Q3 关注点）：
  · mpi = 区间端点做普通任意精度算术，误差靠"输入区间宽度"传播，不保证每步舍入被包住；
  · arb = Arb 库的 ball arithmetic (mid +/- rad)，每步运算按定向舍入扩展半径，
          结果是数学上严格的包含区间（correct rounding + validated bounds）。

常数：CODATA 2022（Tiesinga, Mohr, Newell, Taylor, arXiv:2409.03787, Web 9.0）
  alpha = 7.2973525643(11)e-3, mu0 = 1.25663706127(20)e-6, G = 6.67430(15)e-11,
  m_e = 9.1093837139(28)e-31 kg; c, h, e 精确。

用法: python q3_arb_certified.py
输出: stdout + q3_arb_results.json
"""
import json
from flint import arb, ctx

ctx.prec = 128          # 中点 128 bit ~ 38 位十进制；半径内部 30 bit 上取整（保守）


def lo(x):
    """区间的严格下界（转 float 用于报告；判定一律用 arb 原生 overlaps/contains）。"""
    return float(x.lower())


def hi(x):
    return float(x.upper())


OUT = {}

# ---------- 1. 输入球（中心值 +/- 标准不确定度） ----------
alpha = arb("7.2973525643e-3", "1.1e-12")
mu0 = arb("1.25663706127e-6", "2.0e-16")
G = arb("6.67430e-11", "1.5e-15")
c = arb("299792458")                 # 精确
h = arb("6.62607015e-34")            # 精确
hbar = h / (2 * arb.pi())
me = arb("9.1093837139e-31", "2.8e-40")

OUT["inputs"] = {"alpha": str(alpha), "mu0": str(mu0), "G": str(G),
                 "hbar": str(hbar), "me": str(me), "prec_bits": ctx.prec}

# ---------- 2. alpha^2 mu_0 与 G：Arb 原生不相交判定 ----------
a2mu0 = alpha * alpha * mu0
diff = a2mu0 - G
OUT["alpha2mu0"] = str(a2mu0)
OUT["alpha2mu0_interval"] = [lo(a2mu0), hi(a2mu0)]
OUT["G_interval"] = [lo(G), hi(G)]
OUT["disjoint_arb_native"] = bool(not a2mu0.overlaps(G))
OUT["diff_minus_G"] = str(diff)
OUT["diff_lower_bound"] = lo(diff)
OUT["diff_strictly_positive"] = bool(lo(diff) > 0)
OUT["gap_between_intervals"] = [lo(a2mu0) - hi(G)]
OUT["gap_in_units_of_G_uncertainty"] = (lo(a2mu0) - hi(G)) / float(G.rad())

# ---------- 3. 相对偏差（认证区间） ----------
rel = (a2mu0 - G) / G * 100
OUT["rel_dev_percent"] = str(rel)
OUT["rel_dev_percent_interval"] = [lo(rel), hi(rel)]

# sigma 水平（审稿人 Q4 口径：相对 G 的标准不确定度）
dev_mid = float(abs((a2mu0 - G) / G))
OUT["sigma_vs_G_u"] = dev_mid / float(G.rad() / G.mid())

# ---------- 4. Lambda vs alpha^{-1}：能否"恰好相等" ----------
Lam = 4 * arb.pi() ** 3 + arb.pi() ** 2 + arb.pi()
ainv = 1 / alpha
dlt = Lam - ainv
OUT["Lambda"] = str(Lam)
OUT["Lambda_interval"] = [lo(Lam), hi(Lam)]
OUT["alpha_inv_interval"] = [lo(ainv), hi(ainv)]
OUT["Lambda_minus_alpha_inv"] = str(dlt)
OUT["Lambda_equals_alpha_inv_possible"] = bool(dlt.contains(arb(0)))

a2L2 = alpha * alpha * Lam * Lam
OUT["alpha2_Lambda2"] = str(a2L2)
OUT["alpha2_Lambda2_equals_1_possible"] = bool(a2L2.contains(arb(1)))

# ---------- 5. 无量纲不变量 Pi1 = G me^2 / (hbar c) ----------
Pi1 = G * me * me / (hbar * c)
OUT["Pi1"] = str(Pi1)
OUT["Pi1_interval"] = [lo(Pi1), hi(Pi1)]

# ---------- 6. M* c^2 (GeV) ----------
Mp = (hbar * c / G) ** arb("1/2")
Mstar_GeV = (Mp / Lam) * c * c / arb("1.602176634e-10")
OUT["Mstar_GeV"] = str(Mstar_GeV)
OUT["Mstar_GeV_interval"] = [lo(Mstar_GeV), hi(Mstar_GeV)]

print(json.dumps(OUT, indent=2, ensure_ascii=False))
with open("q3_arb_results.json", "w", encoding="utf-8") as f:
    json.dump(OUT, f, indent=2, ensure_ascii=False)
