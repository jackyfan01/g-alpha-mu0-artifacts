# -*- coding: utf-8 -*-
"""Q8: 大统一能标 M_X 的重整化群拟合 —— 把"模型能标"从文献量级估计换成算出来的数。

外审第 6 轮 #6 点名要求：挑一到两个"标准"统一场景，用**重整化群跑动 + 阈值效应**
（并带不确定度）说明 M* 到底对不对得上；第 7 轮"缺失相关工作"又说 GUT 能标的
模型依赖与不确定性要写实。本文原表 tab:gut 用的是"文献常引量级估计"（论文 §5.5
局限 (v) 自己也承认了这点）。本脚本把那一条补成**可复现的一圈（1-loop）RGE 计算**，
并对阈值与输入不确定度分别扫。

--------------------------------------------------------------------------
方法：一圈规范耦合跑动
--------------------------------------------------------------------------
GUT 归一化下取三耦合 α_1 = (5/3)α_Y、α_2 = α_W、α_3 = α_s。一圈 RGE：
        α_i^{-1}(μ) = α_i^{-1}(μ_0) − (b_i / 2π)·ln(μ/μ_0)
标准系数（综述见 nath2007；教科书级）：
        SM  : b = (41/10, −19/6, −7)
        MSSM: b = (33/5,  1,    −3)
输入取 PDG 2022 MS-bar 值（在 M_Z 处）：
        α_em^{-1}(M_Z) = 127.952 ± 0.009
        sin²θ_W(M_Z)   = 0.23122 ± 0.00004     （MS-bar）
        α_s(M_Z)       = 0.1179  ± 0.0009
由 α_1^{-1} = (3/5)cos²θ_W·α_em^{-1}、α_2^{-1} = sin²θ_W·α_em^{-1}、α_3^{-1} = 1/α_s 定出 M_Z 处三耦合。

场景 A（MSSM，SUSY SU(5)/SO(10) 型）：M_Z→M_SUSY 用 SM 系数，M_SUSY 以上用 MSSM 系数，
解 α_1(M_X) = α_2(M_X) 给统一能标；M_SUSY 取阈值扫描（不确定的主要来源）。
场景 B（最小非超对称 SU(5)）：全程 SM 系数，检验三耦合能否相交。

--------------------------------------------------------------------------
诚实边界（主动声明）
--------------------------------------------------------------------------
· 一圈计算。已知二圈与 GUT 阈值匹配会把 M_X 整体上移数十个百分点量级（相对改变，
  不改变量级）；本文**不**在此脚本内伪造二圈系数，只把该理论不确定性显式列为附加带。
· 统一能标 M_X 不是可观测量：它随 SUSY 阈值移动达约 3 倍，远超输入测量误差（约 1%）。
· 结论只依赖"量级 + 模型依赖"，不依赖 M_X 的精确值。
"""

import json
import os
from mpmath import mp, mpf, log, exp, sqrt, pi

mp.dps = 50

HERE = os.path.dirname(os.path.abspath(__file__))
OUT_PATH = os.path.normpath(os.path.join(HERE, "..", "results", "q8_gut_running.json"))

# ------------------------------------------------------------------ 输入
MZ = mpf("91.1876")
AEM_INV = mpf("127.952")
AEM_INV_U = mpf("0.009")
SIN2W = mpf("0.23122")
SIN2W_U = mpf("0.00004")
ALPHAS = mpf("0.1179")
ALPHAS_U = mpf("0.0009")

COS2W = 1 - SIN2W
A1_INV = mpf("0.6") * COS2W * AEM_INV      # (3/5) cos²θ_W α_em^{-1}
A2_INV = SIN2W * AEM_INV                    # sin²θ_W α_em^{-1}
A3_INV = 1 / ALPHAS

# 一圈 beta 系数（GUT 归一化）
B_SM = [mpf(41) / 10, mpf(-19) / 6, mpf(-7)]
B_MSSM = [mpf(33) / 5, mpf(1), mpf(-3)]
TWOPI = 2 * pi


def run(a_inv, b, t):
    """把三个 α_i^{-1} 跑动 t = ln(μ/μ0)。"""
    return [a_inv[i] - b[i] / TWOPI * t for i in range(3)]


def mssm_unification(msusy, a1_inv_mz=None, a2_inv_mz=None, a3_inv_mz=None):
    """解 α_1(M_X)=α_2(M_X)，返回 (M_X, 三耦合@M_X)。"""
    a1 = A1_INV if a1_inv_mz is None else a1_inv_mz
    a2 = A2_INV if a2_inv_mz is None else a2_inv_mz
    a3 = A3_INV if a3_inv_mz is None else a3_inv_mz
    t_s = log(msusy / MZ)
    # M_Z → M_SUSY（SM 段）
    a1_s, a2_s, a3_s = run([a1, a2, a3], B_SM, t_s)
    # M_SUSY 以上（MSSM 段）：写成 α_i^{-1}(μ) = A_i − c_i t
    c = [b / TWOPI for b in B_MSSM]
    A = [a1_s + c[0] * t_s, a2_s + c[1] * t_s, a3_s + c[2] * t_s]
    t_x = (A[0] - A[1]) / (c[0] - c[1])          # α_1 = α_2
    m_x = MZ * exp(t_x)
    a_at_x = [A[i] - c[i] * t_x for i in range(3)]
    return m_x, a_at_x


def sm_crossing():
    """场景 B：全程 SM 系数，解 α_1=α_2，返回 (μ_cross, 三耦合@μ_cross)。"""
    c = [b / TWOPI for b in B_SM]
    # α_i^{-1}(t) = α_i^{-1}(M_Z) − c_i t ；令 1=2：t = (α1 − α2)/(c1 − c2)
    t_x = (A1_INV - A2_INV) / (c[0] - c[1])
    mu = MZ * exp(t_x)
    a_at = [A1_INV - c[0] * t_x, A2_INV - c[1] * t_x, A3_INV - c[2] * t_x]
    return mu, a_at


# ---------------------------------------------------- 输入耦合与不确定度
# D = α_1^{-1} − α_2^{-1} = (0.6cos²θ_W − sin²θ_W)·α_em^{-1}；其对输入的偏导
dD_daem = mpf("0.6") * COS2W - SIN2W
dD_dsin2w = -(mpf("1.6")) * AEM_INV            # ∂(D)/∂sin²θ_W = −0.6α_em^{-1} − α_em^{-1}
u_D = sqrt((dD_daem * AEM_INV_U) ** 2 + (dD_dsin2w * SIN2W_U) ** 2)

# 结论
MSTAR_C2 = mpf("8.909e16")      # GeV，M*c² = M_P c²/Λ（本文定理 5，Arb 区间中值）

MSUSY_SCAN_TEV = [mpf("0.3"), mpf("0.5"), mpf("1"), mpf("2"), mpf("2.8"), mpf("5"), mpf("10")]
MSUSY_TEV = mpf("1")

sweep = []
for mtev in MSUSY_SCAN_TEV:
    m_x, a_at = mssm_unification(mtev * 1000)
    sweep.append({
        "M_SUSY_TeV": float(mtev),
        "M_X_GeV": float(m_x),
        "alpha1_inv_at_MX": float(a_at[0]),
        "alpha2_inv_at_MX": float(a_at[1]),
        "alpha3_inv_at_MX": float(a_at[2]),
        "unification_quality_a3_minus_a1": float(a_at[2] - a_at[0]),
        "ratio_MstarC2_over_MX": float(MSTAR_C2 / m_x),
    })

MX_NOMINAL, A_AT_MX = mssm_unification(MSUSY_TEV * 1000)

# 输入不确定度 → M_X 的相对不确定度：t_X = (A1−A2)/(c1−c2)，u(D)=u(A1−A2)
c_mssm = [b / TWOPI for b in B_MSSM]
u_tx = u_D / (c_mssm[0] - c_mssm[1])
u_MX_rel = u_tx                                  # u(M_X)/M_X = u(t_X)

# 阈值带（扫描端点）
MX_LO = min(s["M_X_GeV"] for s in sweep)
MX_HI = max(s["M_X_GeV"] for s in sweep)

# 场景 B
MU_SM, A_SM = sm_crossing()

result = {
    "method": "1-loop gauge-coupling RGE; GUT-normalized a1=(3/5)aY, a2=aW, a3=as",
    "inputs_at_MZ": {
        "M_Z_GeV": float(MZ),
        "alpha_em_inv": float(AEM_INV), "alpha_em_inv_u": float(AEM_INV_U),
        "sin2thetaW_MS": float(SIN2W), "sin2thetaW_MS_u": float(SIN2W_U),
        "alpha_s": float(ALPHAS), "alpha_s_u": float(ALPHAS_U),
        "alpha1_inv": float(A1_INV), "alpha2_inv": float(A2_INV), "alpha3_inv": float(A3_INV),
    },
    "beta_coefficients_1loop": {
        "SM": [float(b) for b in B_SM],
        "MSSM": [float(b) for b in B_MSSM],
    },
    "scenario_A_MSSM": {
        "definition": "SM below M_SUSY, MSSM above; solve alpha1(M_X)=alpha2(M_X)",
        "M_SUSY_threshold_sweep": sweep,
        "M_X_range_GeV": [MX_LO, MX_HI],
        "M_X_nominal_GeV": float(MX_NOMINAL),
        "couplings_at_MX_nominal": {
            "alpha1_inv": float(A_AT_MX[0]),
            "alpha2_inv": float(A_AT_MX[1]),
            "alpha3_inv": float(A_AT_MX[2]),
        },
        "input_uncertainty_on_MX_rel": float(u_MX_rel),
        "threshold_spread_factor_MX_HI_over_LO": MX_HI / MX_LO,
    },
    "scenario_B_minimal_SU5": {
        "definition": "SM running throughout; solve alpha1=alpha2 for crossing",
        "crossing_scale_GeV": float(MU_SM),
        "couplings_at_crossing": {
            "alpha1_inv": float(A_SM[0]),
            "alpha2_inv": float(A_SM[1]),
            "alpha3_inv": float(A_SM[2]),
        },
        "alpha3_defect_at_crossing": float(A_SM[2] - A_SM[0]),
        "unifies": False,
    },
    "comparison_to_Mstar": {
        "Mstar_c2_GeV": float(MSTAR_C2),
        "MstarC2_over_MX_nominal": float(MSTAR_C2 / MX_NOMINAL),
        "MstarC2_over_MX_range": [float(MSTAR_C2 / MX_HI), float(MSTAR_C2 / MX_LO)],
        "Mstar_inside_MX_range": bool(MX_LO <= float(MSTAR_C2) <= MX_HI),
    },
    "honest_scope": {
        "loop_order": "1-loop only; 2-loop + GUT-threshold matching raise M_X by ~tens of percent (see two_loop_crosscheck)",
        "two_loop_crosscheck": {
            "source": "Bourilkov, JHEP 11 (2015) 117 (arXiv:1508.04176)",
            "M_GUT_GeV": "1.05e16 +/- 0.07e16",
            "M_SUSY_GeV": "2820 +650/-530",
            "note": "1-loop here at M_SUSY=2.8 TeV gives ~7.5e15; 2-loop rises it by ~40%. Even so M*c^2/M_GUT ~ 8.5.",
        },
        "dominant_uncertainty": "SUSY threshold M_SUSY (factor ~%.1f) >> input measurement error (~%.1f%%)" % (float(MX_HI / MX_LO), float(u_MX_rel * 100)),
    },
    "verdict": (
        "M*c^2 = 8.909e16 GeV lies ABOVE the entire computed MSSM unification range "
        "[%.2e, %.2e] GeV (factor %.1f-%.1f); the 'GUT-scale correspondence' does not "
        "survive RG running + thresholds. Dominant spread is the unknown SUSY threshold, "
        "not measurement error: M_X is not an observable."
        % (MX_LO, MX_HI, float(MSTAR_C2 / MX_HI), float(MSTAR_C2 / MX_LO))
    ),
}

with open(OUT_PATH, "w", encoding="utf-8") as f:
    json.dump(result, f, ensure_ascii=False, indent=2)

print("wrote", OUT_PATH)
print("inputs: a1^-1=%.4f a2^-1=%.4f a3^-1=%.4f" % (A1_INV, A2_INV, A3_INV))
print("M_X nominal (M_SUSY=1TeV) = %.4e GeV" % MX_NOMINAL)
print("couplings at M_X: (%.3f, %.3f, %.3f)" % (A_AT_MX[0], A_AT_MX[1], A_AT_MX[2]))
print("M_X range over M_SUSY in [0.3,10] TeV = [%.4e, %.4e] GeV" % (MX_LO, MX_HI))
print("threshold spread factor = %.3f ; input u(M_X)/M_X = %.4f%%" % (MX_HI / MX_LO, u_MX_rel * 100))
print("M*c^2/M_X nominal = %.3f ; range = [%.2f, %.2f]" % (MSTAR_C2 / MX_NOMINAL, MSTAR_C2 / MX_HI, MSTAR_C2 / MX_LO))
print("SM scenario: crossing at %.4e GeV, a3 defect = %.3f" % (MU_SM, A_SM[2] - A_SM[0]))
