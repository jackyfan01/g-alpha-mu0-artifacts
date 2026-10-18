# -*- coding: utf-8 -*-
"""R4: 区间算术验证（回应外审 Q6）

用 mpmath.mpi 多精度区间对 α²μ₀ 与 G 构造包络，证明二者区间不相交，
即 0.26% 偏差不是测量/舍入误差的产物，而是给定常数下的严格事实。

诚实性说明：mpmath.mpi 提供误差传播包络，但不是像 Arb/MPFR 那样
逐舍入方向认证的严格区间；本验证给出的是"传播式包络"，
结论对界内任意扰动稳健。Lean 侧的认证浮点（FloatLib 等）留作未来工作。

常数：CODATA 2022。用法：python R4_interval_arithmetic.py
"""
from pathlib import Path
import json
from mpmath import mp, mpf, mpi, pi, fabs

mp.dps = 50

# 输入球与 q3/q5 完全同口径：一律取 CODATA **展示值的末位**（u(alpha)=1.1e-12，
# u(mu0)=2.0e-16，u(G)=1.5e-15），这样 R4 与 Arb 的差异只反映"传播包络 vs 定向舍入"
# 这一算术形式区别，而不掺杂口径区别。
ALPHA, U_ALPHA = mpf("7.2973525643e-3"), mpf("1.1e-12")     # 7.297 352 5643(11)e-3
MU0, U_MU0 = mpf("1.25663706127e-6"), mpf("2.0e-16")        # 1.256 637 061 27(20)e-6
G, U_G = mpf("6.67430e-11"), mpf("1.5e-15")                 # 6.674 30(15)e-11

A = mpi(str(ALPHA - U_ALPHA), str(ALPHA + U_ALPHA))
M = mpi(str(MU0 - U_MU0), str(MU0 + U_MU0))
Gg = mpi(str(G - U_G), str(G + U_G))

a2mu0 = A**2 * M
dev = fabs(ALPHA**2 * MU0 - G) / G

out = {
    "alpha2mu0_interval": [mp.nstr(a2mu0.a, 14), mp.nstr(a2mu0.b, 14)],
    "G_interval": [mp.nstr(Gg.a, 14), mp.nstr(Gg.b, 14)],
    "disjoint": bool(a2mu0.a > Gg.b),
    "gap_lower_minus_upper": mp.nstr(a2mu0.a - Gg.b, 6),
    "alpha2mu0_strictly_greater": bool(a2mu0.a > Gg.b),
    "rel_dev_percent": mp.nstr(dev * 100, 8),
    "widths": [mp.nstr(a2mu0.b - a2mu0.a, 6), mp.nstr(Gg.b - Gg.a, 6)],
    "caveat": "mpmath.mpi 为误差传播包络，非逐舍入认证（区别于 Arb/MPFR）；"
              "结论：区间不相交 ⇒ 偏差为严格事实，非测量不确定度所致。",
}
print(json.dumps(out, ensure_ascii=False, indent=2))
with open(str(Path(__file__).resolve().parent.parent / "results" / "mpi_crosscheck.json"), "w", encoding="utf-8") as f:
    json.dump(out, f, ensure_ascii=False, indent=2)
