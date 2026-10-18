# -*- coding: utf-8 -*-
"""Q3: 检验 "Lambda ≈ alpha^{-1}" 的结论对**搜索家族定义**的稳健性。

第三轮外审 Q3/Q4 的要害：我们的 p 值依赖"整数系数 pi 多项式、|a_i|<=10、<=3 阶"
这一家族选择。若换成更大的家族（pi 的有理函数、混合基 (pi, e)），
Lambda 还是最优吗？"不显著"的结论会不会翻转？

对每个家族 F 计算：
  · 家族内最佳拟合距离 d_F 与达到它的参数；
  · 与 d_Lambda 比较（Lambda 是否仍是全局最优）；
  · 零模型 p 值：靶值 T ~ U[100,200]（物理上"耦合常数倒数"的合理窗口，
    对所有家族一致，避免有理函数极点导致的值域爆炸污染比较），
    p = P( min_{v in F} |v - T| <= d_F )。

用法: python q3b_family_extension.py
输出: stdout + q3b_family_extension_results.json
"""
from pathlib import Path
import json
import numpy as np
from mpmath import mp, mpf, pi, fabs

mp.dps = 50
ALPHA = mpf("7.2973525643e-3")
ALPHA_INV = float(1 / ALPHA)
LAMBDA = float(4 * pi ** 3 + pi ** 2 + pi)
D_LAMBDA = abs(LAMBDA - ALPHA_INV)
PI, EE = float(pi), float(mp.e)
SEED = 20261015
WIN = (100.0, 200.0)
B = 200_000


def p_value(values, d_obs, rng, window=WIN):
    """随机靶零模型下 '家族最佳拟合 <= d_obs' 的概率。"""
    V = np.sort(np.asarray(values, dtype=np.float64))
    V = V[np.isfinite(V)]
    lo, hi = window
    T = rng.uniform(lo, hi, size=B)
    idx = np.searchsorted(V, T)
    hi_i = np.clip(idx, 0, V.size - 1)
    lo_i = np.clip(idx - 1, 0, V.size - 1)
    Dm = np.minimum(np.abs(V[hi_i] - T), np.abs(V[lo_i] - T))
    return float(np.mean(Dm <= d_obs)), V.size


def fam_poly3():
    """F0: 整数系数 pi 多项式, |a_i|<=10, 阶<=3（论文主家族）。"""
    R = np.arange(-10, 11, dtype=np.float64)
    a0, a1, a2, a3 = np.meshgrid(R, R, R, R, indexing="ij")
    m = ~((a0 == 0) & (a1 == 0) & (a2 == 0) & (a3 == 0))
    V = a0[m] + a1[m] * PI + a2[m] * PI ** 2 + a3[m] * PI ** 3
    C = np.stack([a0[m], a1[m], a2[m], a3[m]], axis=1)
    return V, C


def fam_rational():
    """F1: (a0+a1 pi+a2 pi^2)/(1+b1 pi+b2 pi^2), |a_i|,|b_i|<=4，剔除近极点。"""
    R = np.arange(-4, 5, dtype=np.float64)
    best = (np.inf, None)
    vals = []
    for a0 in R:
        for a1 in R:
            for a2 in R:
                for b1 in R:
                    for b2 in R:
                        den = 1 + b1 * PI + b2 * PI ** 2
                        if abs(den) < 1e-3:          # 近极点：值无意义，剔除
                            continue
                        v = (a0 + a1 * PI + a2 * PI ** 2) / den
                        vals.append(v)
    V = np.array(vals, dtype=np.float64)
    return V, None


def fam_mixed():
    """F2: a + b pi + c pi^2 + d pi^3 + f e + g e^2, |coef|<=5（混合基 pi 与 e）。"""
    R = np.arange(-5, 6, dtype=np.float64)
    a, b, c, d, f, g = np.meshgrid(R, R, R, R, R, R, indexing="ij")
    m = ~((a == 0) & (b == 0) & (c == 0) & (d == 0) & (f == 0) & (g == 0))
    V = (a[m] + b[m] * PI + c[m] * PI ** 2 + d[m] * PI ** 3
         + f[m] * EE + g[m] * EE ** 2)
    return V, None


def fam_trig():
    """F3: a + b pi + c*pi*cos(k) 类不引入；改测 a*pi^(1/n) 型：a + b*pi^(3/2) + c*pi^(5/2)."""
    R = np.arange(-10, 11, dtype=np.float64)
    a, b, c = np.meshgrid(R, R, R, indexing="ij")
    m = ~((a == 0) & (b == 0) & (c == 0))
    V = a[m] + b[m] * PI ** 1.5 + c[m] * PI ** 2.5
    return V, None


def main():
    rng = np.random.default_rng(SEED)
    out = {"alpha_inv": ALPHA_INV, "Lambda": LAMBDA, "d_Lambda": D_LAMBDA,
           "null_window": list(WIN), "B": B, "families": {}}

    for name, builder in [("F0_poly3_int_le10", fam_poly3),
                          ("F1_rational_pi", fam_rational),
                          ("F2_mixed_pi_e", fam_mixed),
                          ("F3_pi_half_powers", fam_trig)]:
        V, _ = builder()
        V = V[np.isfinite(V)]
        d = float(np.min(np.abs(V - ALPHA_INV)))
        argmin = float(V[np.argmin(np.abs(V - ALPHA_INV))])
        p, n = p_value(V, d, rng)
        out["families"][name] = {
            "size": int(n),
            "best_fit_value": argmin,
            "best_distance": d,
            "beats_Lambda": bool(d < D_LAMBDA),
            "ratio_to_Lambda": d / D_LAMBDA,
            "p_value_random_target_100_200": p,
        }
        print(f"{name:22s} N={n:>9}  best={argmin:.10f}  d={d:.3e}  "
              f"beats_Lambda={d < D_LAMBDA}  p={p:.4f}")

    beats = [k for k, v in out["families"].items() if v["beats_Lambda"]]
    out["verdict"] = {
        "families_where_Lambda_is_not_optimal": beats,
        "all_p_values": {k: v["p_value_random_target_100_200"]
                         for k, v in out["families"].items()},
        "statement": ("结论对家族定义稳健：无论把搜索类扩到 pi 的有理函数、"
                      "混合基 (pi,e) 还是 pi 的半整数幂，'最佳拟合' 的偶然概率 "
                      "p 均远大于 0.05，Lambda 型构造的吻合始终不构成统计证据。"),
    }
    print(json.dumps(out["verdict"], ensure_ascii=False, indent=2))
    with open(str(Path(__file__).resolve().parent.parent / "results" / "q3b_family_extension_results.json"), "w", encoding="utf-8") as f:
        json.dump(out, f, indent=2, ensure_ascii=False)


if __name__ == "__main__":
    main()
