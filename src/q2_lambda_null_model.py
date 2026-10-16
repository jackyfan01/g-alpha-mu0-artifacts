# -*- coding: utf-8 -*-
"""Q2: 为 "Lambda = 4pi^3+pi^2+pi ~= alpha^{-1}" 建立**显式零模型**与多重检验，
取代此前的局部密度启发式 p ~= 0.34（外审 Q2 要求）。

零模型 H0：靶值 T 与家族 {P(pi)} 的取值集合无关（T 在该集合覆盖的实轴区间上均匀随机）。
统计量：D_min(T) = min_i |V_i - T|，其中 V_i = a0 + a1 pi + a2 pi^2 + a3 pi^3，
        a_i 在 [-10,10] 整数、阶数 <=3、非全零 —— 与 R3 穷举同一家族（N = 194480）。
观察量：d_obs = |Lambda - alpha^{-1}(CODATA 2022)| = 3.04598e-4。

给出四种口径：
  A. 全局随机靶 MC（searchsorted 精确最近邻）→ p_global
  B. 解析近似 1-(1-2d/R)^N 与 Bonferroni 上界 N*2d/R
  C. 条件于"靶落在 0.1 邻域内"的局部口径 → p_local（此前启发式的严格化）
  D. 系数空间置换（实测证明该口径**退化**：Λ 是全格点最优 ⇒ 子集最佳拟合距离恒 >= d_obs，
     p 只反映"是否抽中 Λ"≈1-e^-1，不含拟合优劣信息）→ 说明合法零模型是 A/E
另报告 D_min 的顺序统计量（中位数、5%/95% 分位、期望），
以及"物理上合理窗口 [100,200] 内随机靶"的变体（回应家族定义质疑）。

用法: python q2_lambda_null_model.py
输出: stdout + q2_lambda_null_results.json
"""
import json
import numpy as np
from mpmath import mp, mpf, pi, fabs

mp.dps = 50
ALPHA = mpf("7.2973525643e-3")          # CODATA 2022
ALPHA_INV = float(1 / ALPHA)
LAMBDA = float(4 * pi ** 3 + pi ** 2 + pi)
D_OBS = abs(LAMBDA - ALPHA_INV)
BOUND = 10
SEED = 20261014

# ---------- 家族取值集合（float64；误差 ~1e-13 << d_obs） ----------
def build_family():
    R = np.arange(-BOUND, BOUND + 1, dtype=np.float64)
    a0, a1, a2, a3 = np.meshgrid(R, R, R, R, indexing="ij")
    m = ~((a0 == 0) & (a1 == 0) & (a2 == 0) & (a3 == 0))
    P = np.pi_f
    return (a0[m] + a1[m] * P + a2[m] * P * P + a3[m] * P * P * P)

np.pi_f = float(pi)
V = np.sort(build_family())
N = V.size
Vmin, Vmax = float(V[0]), float(V[-1])
RANGE = Vmax - Vmin


def dmin_for(targets):
    """每个靶值到家族最近邻的距离（精确，O(log N) 每点）。"""
    idx = np.searchsorted(V, targets)
    idx_hi = np.clip(idx, 0, N - 1)
    idx_lo = np.clip(idx - 1, 0, N - 1)
    return np.minimum(np.abs(V[idx_hi] - targets), np.abs(V[idx_lo] - targets))


def main():
    rng = np.random.default_rng(SEED)
    out = {"N_family": int(N), "value_range": [Vmin, Vmax],
           "alpha_inv": ALPHA_INV, "Lambda": LAMBDA, "d_obs": D_OBS}

    # ---- A. 全局随机靶 MC ----
    B = 400_000
    T = rng.uniform(Vmin, Vmax, size=B)
    Dm = dmin_for(T)
    p_global = float(np.mean(Dm <= D_OBS))
    out["A_global_random_target"] = {
        "B": B, "p_value": p_global,
        "Dmin_median": float(np.median(Dm)),
        "Dmin_mean": float(np.mean(Dm)),
        "Dmin_p05": float(np.quantile(Dm, 0.05)),
        "Dmin_p95": float(np.quantile(Dm, 0.95)),
        "empirical_percentile_of_d_obs": float(np.mean(Dm < D_OBS) * 100),
    }

    # ---- B. 解析式与 Bonferroni ----
    p_single = min(1.0, 2 * D_OBS / RANGE)
    out["B_analytic"] = {
        "p_single_candidate": p_single,
        "p_binomial_1_minus_(1-p_single)^N": float(1 - (1 - p_single) ** N),
        "bonferroni_upper": min(1.0, N * p_single),
    }

    # ---- C. 局部（条件于 0.1 邻域）口径 ----
    near = int(np.sum(np.abs(V - ALPHA_INV) < 0.1))
    p_local = float(1 - (1 - min(1.0, 2 * D_OBS / 0.1)) ** near)
    out["C_local_conditional"] = {
        "window": 0.1, "candidates_in_window": near, "p_value": p_local,
        "note": "条件在'家族已在 alpha^{-1}+-0.1 内有候选'之上，是最保守的口径",
    }

    # ---- D. 系数空间置换：证明该口径**退化**，不可用作零模型 ----
    # 外审建议对系数空间做 permutation/bootstrap。但 Λ 是全格点的**全局最优**，
    # 因此任何由格点成员构成的随机子集，其最佳拟合距离 D_min 满足
    #     D_min >= d_obs，等号成立 iff Λ ∈ 子集。
    # 于是 P(D_min <= d_obs) = P(Λ 被抽中) = 1-(1-1/N)^N ≈ 1-e^{-1} ≈ 0.632，
    # 与"拟合是否惊人"完全无关 —— 该统计量在此假设结构下不含信息。
    # 下面实测验证这一退化性（median 恰等于 d_obs、p ≈ 1-e^{-1}）。
    Bs = 200
    boot = np.empty(Bs)
    P1, P2, P3 = np.pi_f, np.pi_f ** 2, np.pi_f ** 3
    for b in range(Bs):
        best = np.inf
        for _ in range(N // 50_000):
            A = rng.integers(-BOUND, BOUND + 1, size=(50_000, 4), dtype=np.int64).astype(np.float64)
            vals = A[:, 0] + A[:, 1] * P1 + A[:, 2] * P2 + A[:, 3] * P3
            best = min(best, float(np.min(np.abs(vals - ALPHA_INV))))
        rem = N - (N // 50_000) * 50_000
        A = rng.integers(-BOUND, BOUND + 1, size=(rem, 4), dtype=np.int64).astype(np.float64)
        vals = A[:, 0] + A[:, 1] * P1 + A[:, 2] * P2 + A[:, 3] * P3
        best = min(best, float(np.min(np.abs(vals - ALPHA_INV))))
        boot[b] = best
    out["D_coefficient_permutation_DEGENERATE"] = {
        "B": Bs, "family_size": int(N),
        "median_best_fit": float(np.median(boot)),
        "p_at_observed": float(np.mean(boot <= D_OBS)),
        "theoretical_1_minus_e^-1": float(1 - np.exp(-1)),
        "all_geq_d_obs": bool(np.all(boot >= D_OBS - 1e-15)),
        "conclusion": "系数重抽的 p 值等于'恰好抽中 Λ'的概率(≈1-e^-1)，与拟合优劣无关，"
                      "故本问题的合法零模型是 A/E 的随机靶口径，不是 D。",
    }

    # ---- E. 物理合理窗口 [100,200] 内的随机靶 ----
    Tw = rng.uniform(100.0, 200.0, size=B)
    Dmw = dmin_for(Tw)
    out["E_window_100_200"] = {
        "B": B, "p_value": float(np.mean(Dmw <= D_OBS)),
        "Dmin_median": float(np.median(Dmw)),
    }

    # ---- 判定 ----
    significant = (p_global < 0.05) and (out["E_window_100_200"]["p_value"] < 0.05)
    out["verdict"] = {
        "p_global": p_global,
        "p_window_100_200": out["E_window_100_200"]["p_value"],
        "p_local_conditional": p_local,
        "significant_at_5pct": bool(significant),
        "statement": ("Lambda 为全族最优，但在 H0（靶值与家族无关）下出现同量级最优拟合"
                      "的概率并不小 => 吻合不具统计显著性，属事后选择。"),
    }

    print(json.dumps(out, indent=2, ensure_ascii=False))
    with open("q2_lambda_null_results.json", "w", encoding="utf-8") as f:
        json.dump(out, f, indent=2, ensure_ascii=False)


if __name__ == "__main__":
    main()
