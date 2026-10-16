# -*- coding: utf-8 -*-
r"""q5b：把 α 端球从"展示末位 1.1e-12"抬到 CODATA 展示的 u_r=1.6e-10（u=1.1676e-12），
   作为**上界压力测试**重跑同一条完全相关 Arb trace，检验"不相交"判定是否依赖口径选择。

为什么需要它：由 mu0 = 2*alpha*h/(e^2 c) 且 h/e/c 精确，未舍入世界里 u_r(mu0) == u_r(alpha)。
CODATA 把 u(mu0) 展示为 (20) 即 2.0e-16，两位舍入要能得到它，公共未舍入 u_r 必须 >= 1.552e-10，
相应 u(alpha) >= 1.132e-12 —— 即本文采用的展示末位 1.1e-12 可能偏窄 3%~5%。
本脚本量化这一点。
"""
import io
import json
import os

from flint import arb, ctx

ctx.prec = 128

HERE = os.path.dirname(os.path.abspath(__file__))
OUT_PATH = os.path.normpath(os.path.join(HERE, "..", "results",
                                         "q5b_alpha_ball_nist_caliber.json"))

H = arb("6.62607015e-34")            # 精确
E = arb("1.602176634e-19")           # 精确
C = arb("299792458")                 # 精确
G = arb("6.67430e-11", "1.5e-15")    # 6.674 30(15)e-11

CASES = {
    "displayed_last_digit": arb("7.2973525643e-3", "1.1e-12"),      # 正文主口径
    "nist_relative": arb("7.2973525643e-3", "1.1676e-12"),          # u_r = 1.6e-10
    "lower_bound_from_mu0_paren": arb("7.2973525643e-3", "1.1324e-12"),
}

out = {"note": "u(alpha) 三种取法下的完全相关证书；G 端球不变。", "cases": {}}
for name, a in CASES.items():
    mu0 = (2 * a * H) / (E * E * C)
    a2mu0 = a * a * mu0
    lo, hi = float(a2mu0.lower()), float(a2mu0.upper())
    gap = lo - float(G.upper())
    out["cases"][name] = {
        "u_alpha": float(a.rad()),
        "u_r_alpha": float(a.rad() / a.mid()),
        "alpha2mu0_interval": [lo, hi],
        "u_r": (hi - lo) / 2.0 / float(a2mu0.mid()),
        "disjoint": bool(not a2mu0.overlaps(G)),
        "gap": gap,
        "gap_in_uG": gap / float(G.rad()),
        "sigma": abs((float(a2mu0.mid()) - float(G.mid())) / float(G.rad())),
    }
main = out["cases"]["displayed_last_digit"]["u_r"]
for k, v in out["cases"].items():
    v["width_vs_main"] = v["u_r"] / main - 1.0

from mpmath import mp
# 关键对照：把 mu0 端球换成由 alpha 传递的半径，"独立"口径还剩多少加宽？
mu_prop = (2 * a * H) / (E * E * C)                       # 传递值球（半径 = u(mu0)）
a_disp = arb("7.2973525643e-3", "1.1e-12")
# 关键：中点仍用 CODATA **展示值** 1.25663706127e-6，只把半径换成传递值 1.89428e-16，
# 这样唯一被隔离出来的变量就是"展示舍入"。
ind_prop = a_disp * a_disp * arb("1.25663706127e-6", "1.89428e-16")
a2mu0_corr_main = 2 * a_disp ** 3 * H / (E * E * C)
out["independence_vs_display"] = {
    "u_r_indep_with_propagated_radius": float(ind_prop.rad() / ind_prop.mid()),
    "u_r_correlated": float(a2mu0_corr_main.rad() / a2mu0_corr_main.mid()),
    "relative_difference": float(ind_prop.rad() / ind_prop.mid()
                                 / (a2mu0_corr_main.rad() / a2mu0_corr_main.mid()) - 1.0),
    "u_r_indep_with_displayed_20": float(
        (a_disp * a_disp * arb("1.25663706127e-6", "2.0e-16")).rad()
        / float((a_disp * a_disp * arb("1.25663706127e-6", "2.0e-16")).mid())),
    "reading": ("独立球比完全相关宽 1.86%，其中独立性本身只值 |relative_difference|（约 1e-5 相对），"
                "余下约 1.85% 全部来自 CODATA 把 u(mu0) 向上舍入展示成 (20)。"),
}

print(json.dumps(out, ensure_ascii=False, indent=1))
io.open(OUT_PATH, "w", encoding="utf-8").write(json.dumps(out, ensure_ascii=False, indent=1) + "\n")
print("写出:", OUT_PATH)
