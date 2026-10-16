# -*- coding: utf-8 -*-
"""Q7: MDL / 描述长度先验 —— 把"家族规模惩罚"综合进单一统计量。

外审第 6 轮 #4 与第 7 轮 #3 都点名要求：给整数系数家族加一个 MDL（或 Solomonoff /
Kolmogorov 代理）先验，把"家族越大、越容易撞上"这件事折进一个统计量里，并问这会不会
定性改变 Λ 的结论。q2 已经用**随机靶零模型**的 p 值回答过同一问题（p 不显著 ⇒ 事后选择）；
本脚本用**两段码描述长度**独立复算，并显式计入"穷举了多少候选"这一选择代价。

--------------------------------------------------------------------------
方法：Rissanen 两段码（two-part code / crude MDL）
--------------------------------------------------------------------------
设想要把 α⁻¹ 传给一个知道你所用**家族**的接收方，精度 ε := u(α⁻¹)（CODATA 2022）。

  H0（无结构）：直接传 α⁻¹ 这个实数到精度 ε。代价
        L0 = log2(1/ε)                                      （比特）
  H1（用品公式）：先传公式的下标（代价 L_p 比特），再传残差
        δ := α⁻¹ − V_k，其中 V_k 是家族成员、|δ| = d_obs。
      残差是有符号量、落在宽 2·d_obs 的窗口里，仍需精度 ε，故残差代价
        log2(2·d_obs/ε)                                     （比特）
      合计
        L1 = L_p + log2(2·d_obs/ε)

  长度差
        ΔL := L0 − L1 = log2(1/(2·d_obs)) − L_p
  ΔL > 0 ⇒ 公式**压缩** ⇒ 支持 H1；ΔL ≤ 0 ⇒ 公式不划算、H0 更简。
  与 Bayes 因子的关系（Barron/Rissanen）：−log2 BF10 ≈ ΔL + O(1)。
  故 |ΔL| ≲ 2 bit 落在 O(1) 含糊区内 ⇒ "无证据"；|ΔL| ≳ 10 bit 才算真证据。

--------------------------------------------------------------------------
O(1) 约定与敏感性（主动声明，不藏）
--------------------------------------------------------------------------
残差窗口取 2·d_obs（有符号）。若改取单边 d_obs，则 ΔL 整体平移 +1 bit；
L0 的 log2(1/ε) 约定本身也只定到 O(1)。两者都在 O(1) 内，**不改变任何结论**：
先验 ΔL 在 −1.5 ~ −2.5 bit（无证据区），计入选择后 ≈ −20 bit（决定性反对）。
脚本对两种约定都报数。

--------------------------------------------------------------------------
选择代价（本脚本的关键增量）
--------------------------------------------------------------------------
L_p 只算了"写出这个公式"的代价，**没算"在多大的一堆公式里挑出它"**。
q2/R3 的家族是 (a0,a1,a2,a3)、|a_i| ≤ 10、阶 ≤ 3、非全零，规模
        N = 21^4 − 1 = 194480
穷举后取最优，等于给了假设 N 次机会。L_p 须再加 log2(N)：
        L_p^eff = L_p + log2(N)
这正是外审要的"把家族规模惩罚折进单一统计量"。

用法: python q7_mdl_prior.py
输出: stdout + q7_mdl_prior.json
"""
import json
import os

from mpmath import mp, mpf, pi, log

mp.dps = 50

HERE = os.path.dirname(os.path.abspath(__file__))
OUT_PATH = os.path.normpath(os.path.join(HERE, "..", "results", "q7_mdl_prior.json"))

# ---------------- 输入：CODATA 2022 ----------------
ALPHA = mpf("7.2973525643e-3")        # 精细结构常数（CODATA 2022）
ALPHA_U = mpf("1.1e-12")              # u(alpha) 绝对不确定度（CODATA 2022）
ALPHA_INV = 1 / ALPHA
ALPHA_INV_U = ALPHA_U / ALPHA ** 2    # u(alpha^{-1}) = u(alpha)/alpha^2

LAMBDA = 4 * pi ** 3 + pi ** 2 + pi    # Λ = 4π³+π²+π
D_OBS = abs(LAMBDA - ALPHA_INV)        # = |Λ − α⁻¹|

BOUND = 10                             # 系数上界 |a_i| ≤ B
DEG = 3                                # 阶数 ≤ 3（4 个系数 a0..a3）
N_FAMILY = (2 * BOUND + 1) ** (DEG + 1) - 1   # 21^4 − 1 = 194480


def log2(x):
    return float(log(mpf(x), 2))


def main():
    eps = ALPHA_INV_U

    # ---------------- 码长 ----------------
    L0 = log2(1 / eps)                 # 直接传 α⁻¹：log2(1/ε)
    resid_signed = log2(2 * D_OBS / eps)   # 有符号残差窗口 2·d_obs（主约定）
    resid_single = log2(D_OBS / eps)       # 单边窗口 d_obs（敏感性）

    # 公式代价 L_p：写出 Λ = 4π³+π²+π
    # 读法 A（最宽松，偏袒 H1）：只数非零系数 {4,1,1}，幂次固定取 1..3
    n_coef_A = 3
    # 读法 B（保守）：4 个系数 a0..a3 各一个整数，a0=0 也要写
    n_coef_B = 4
    L_p_A = n_coef_A * log2(2 * BOUND + 1)
    L_p_B = n_coef_B * log2(2 * BOUND + 1)
    L_search = log2(N_FAMILY)

    def delta(L_p, resid):
        return L0 - (L_p + resid)

    dL_a_priori_A = delta(L_p_A, resid_signed)
    dL_a_priori_B = delta(L_p_B, resid_signed)
    dL_with_search_A = delta(L_p_A + L_search, resid_signed)
    dL_with_search_B = delta(L_p_B + L_search, resid_signed)
    # 敏感性：单边残差窗口
    dL_a_priori_A_single = delta(L_p_A, resid_single)
    dL_with_search_A_single = delta(L_p_A + L_search, resid_single)

    # ---------------- 盈亏平衡：要多准才能压缩 ----------------
    # 先验：L_p < log2(1/(2 d_obs)) ⟺ d_obs < 2^{-L_p}/2
    d_req_apriori = mpf(2) ** (-L_p_A) / 2
    # 计入选择：L_p^eff < log2(1/(2 d_obs))
    d_req_search = mpf(2) ** (-(L_p_A + L_search)) / 2

    # ---------------- 家族界 B 扫描 ----------------
    sweep = []
    for B in range(2, 13):
        lp = n_coef_A * log2(2 * B + 1)
        nf = (2 * B + 1) ** (DEG + 1) - 1
        sweep.append({
            "B": B,
            "N_family": nf,
            "L_p_bits_3coef": lp,
            "compresses_a_priori": bool(lp < log2(1 / (2 * D_OBS))),
            "dL_a_priori_bits": delta(lp, resid_signed),
            "dL_with_selection_bits": delta(lp + log2(nf), resid_signed),
        })
    # 先验压缩的临界 B：3·log2(2B+1) < log2(1/(2 d_obs))
    threshold = log2(1 / (2 * D_OBS)) / 3
    B_crit = (2 ** threshold - 1) / 2

    def bf(dL):
        return 2.0 ** dL

    out = {
        "inputs": {
            "alpha_CODATA2022": str(ALPHA),
            "u_alpha_abs": str(ALPHA_U),
            "alpha_inv": float(ALPHA_INV),
            "u_alpha_inv": float(ALPHA_INV_U),
            "Lambda": float(LAMBDA),
            "d_obs": float(D_OBS),
            "epsilon_used": float(eps),
        },
        "family": {
            "definition": "a0 + a1*pi + a2*pi^2 + a3*pi^3, |a_i|<=10, 阶<=3, 非全零",
            "N_family": N_FAMILY,
            "log2_N": L_search,
        },
        "code_lengths_bits": {
            "L0_direct": L0,
            "residual_signed_2d": resid_signed,
            "residual_single_d": resid_single,
            "L_p_readingA_3coef": L_p_A,
            "L_p_readingB_4coef": L_p_B,
            "L_search": L_search,
            "L_p_eff_A": L_p_A + L_search,
            "L_p_eff_B": L_p_B + L_search,
        },
        "delta_L_bits": {
            "a_priori_A": dL_a_priori_A,
            "a_priori_B": dL_a_priori_B,
            "with_selection_A": dL_with_search_A,
            "with_selection_B": dL_with_search_B,
            "a_priori_A_single_window": dL_a_priori_A_single,
            "with_selection_A_single_window": dL_with_search_A_single,
        },
        "bayes_factor_BF10": {
            "a_priori_A": bf(dL_a_priori_A),
            "with_selection_A": bf(dL_with_search_A),
            "with_selection_B": bf(dL_with_search_B),
            "convention": "BF10 = 2^{dL}; −log2 BF10 ≈ dL + O(1)",
        },
        "breakeven_precision": {
            "d_required_a_priori": float(d_req_apriori),
            "d_required_with_selection": float(d_req_search),
            "d_obs_over_d_required_a_priori": float(D_OBS / d_req_apriori),
            "d_obs_over_d_required_with_selection": float(D_OBS / d_req_search),
            "note": "公式要压缩，偏差须小到什么程度；实测 d_obs 比它大几倍/几个数量级",
        },
        "family_bound_sweep": sweep,
        "B_critical_for_a_priori_compression": {
            "value": B_crit,
            "note": "3·log2(2B+1) < log2(1/(2 d_obs)) 的临界 B；超过则任何家族都不可能先验压缩",
        },
        "verdict": {
            "a_priori": (
                f"ΔL = {dL_a_priori_A:+.2f} bit（读法 A，最偏袒 H1），"
                f"落在 O(1) 含糊区内 ⇒ 与 q2 的“不显著”一致："
                f"不做选择校正时，Λ 的吻合本身不构成压缩证据。"
            ),
            "with_selection": (
                f"计入穷举 {N_FAMILY} 个候选的选择代价后 ΔL = {dL_with_search_A:+.2f} bit"
                f"（BF10 ≈ {bf(dL_with_search_A):.2e}）⇒ 强证据支持“事后选择/巧合”。"
            ),
            "sweep": (
                f"对家族上界 B 扫描：B ≥ {int(B_crit) + 1} 起，"
                f"即便完全不计选择代价也无法先验压缩 ⇒ 结论对家族定义不敏感。"
            ),
            "answers_reviewer": (
                "MDL 会不会定性改变 Λ 的判断？不会 —— 它独立复现并**强化**了 q2 零模型的"
                "结论（事后选择），只是把家族规模惩罚显式写进单一统计量。"
            ),
        },
    }

    print(json.dumps(out, indent=2, ensure_ascii=False))
    with open(OUT_PATH, "w", encoding="utf-8") as f:
        json.dump(out, f, indent=2, ensure_ascii=False)
    print("写出:", OUT_PATH)


if __name__ == "__main__":
    main()
