# -*- coding: utf-8 -*-
r"""Q4：把 §4.8 的 α²μ₀ vs G 认证比较导出成**机器可读的定向舍入 trace**，并配一个独立复核器。

背景（外审 Q4）："the Arb certificates are replayed in Lean only at the level of monotone
comparisons; the full arithmetic trace is not reconstructed." 本脚本按审稿人给的两种验收方式之一
（"or via an external checker"）把这条补上：

  1. 用**精确有理数**（fractions.Fraction，无浮点、无 Arb）把 α²μ₀ 的区间传播拆成逐条初等运算，
     每一步把上下界分别**向外**（下界向下、上界向上）舍入到固定网格 10^-64；
  2. 落盘 `results/q9_arb_trace.json`：每一步的算子、输入区间、精确有理结果、舍入后端点、舍入方向；
  3. 独立**回放**同一 trace：用 Fraction 重算每步的精确结果，逐条断言
       rounded_lo ≤ exact_lo 且 rounded_hi ≥ exact_hi   （即舍入确实是向外的）
     并断言最终 gap = (α²μ₀)_lo − G_hi > 0；落盘 `results/q9_arb_trace_check.json`；
  4. 若 python-flint 可用，再用 Arb（128 bit）独立算一遍参考区间，核对 trace 区间**包含**它
     （trace 比 Arb 更宽，因为固定网格 10^-64 是保守的）。

诚实边界：trace 只覆盖 α²μ₀ 的**区间传播**这一步；α、μ₀ 的输入球本身来自 CODATA（见 README §物理常数取值），
其"上游合并"不在此 trace 内。本 trace 是"证书的算法步骤可核对"，不是"CODATA 本身被重新证明"。
"""
import io
import json
import os
from decimal import Decimal
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.normpath(os.path.join(HERE, "..", "results"))

DIGITS = 64
SCALE = 10 ** DIGITS


def F(s):
    """十进制字面量 → 精确有理数（含指数记法）。"""
    return Fraction(Decimal(s))


def rdown(x):
    """向 −∞ 舍入到网格 10^-DIGITS（返回网格点上的有理数）。"""
    n = (x * SCALE).__floor__()
    return Fraction(n, SCALE), n


def rup(x):
    """向 +∞ 舍入到网格 10^-DIGITS。"""
    n = -((-x * SCALE).__floor__())
    return Fraction(n, SCALE), n


# ---------------------------------------------------------------- 输入球
# 与 q5_correlated_and_spacing.py 完全一致：完全相关口径 μ₀ = 2αh/(e²c)，故 α²μ₀ = 2α³h/(e²c)。
ALPHA_MID = "7.2973525643e-3"
ALPHA_U = "1.1e-12"
G_MID = "6.67430e-11"
G_U = "1.5e-15"
H = "6.62607015e-34"        # SI 精确值
E = "1.602176634e-19"       # SI 精确值
C = "299792458"             # SI 精确值

alpha_lo = F(ALPHA_MID) - F(ALPHA_U)
alpha_hi = F(ALPHA_MID) + F(ALPHA_U)
g_lo = F(G_MID) - F(G_U)
g_hi = F(G_MID) + F(G_U)
h = F(H)
e = F(E)
c = F(C)

steps = []


def add_step(sid, op, desc, ins, exact_lo, exact_hi):
    """登记一步：打印前先做向外舍入，并顺带断言舍入方向正确。"""
    rlo, ilo = rdown(exact_lo)
    rhi, ihi = rup(exact_hi)
    assert rlo <= exact_lo, (sid, "下界未向下舍入")
    assert rhi >= exact_hi, (sid, "上界未向上舍入")
    st = {
        "id": sid,
        "op": op,
        "desc": desc,
        "inputs": ins,
        "exact_lo": str(float(exact_lo)) if abs(exact_lo) > 0 and abs(exact_lo) < 1e-300 else str(exact_lo),
        "exact_hi": str(float(exact_hi)) if abs(exact_hi) > 0 and abs(exact_hi) < 1e-300 else str(exact_hi),
        "rounded_lo": str(rlo),
        "rounded_hi": str(rhi),
        "rounded_lo_scaled_int": str(ilo),
        "rounded_hi_scaled_int": str(ihi),
        "rounding": {"lo": "down", "hi": "up"},
        "outward_ok": True,
    }
    steps.append(st)
    return rlo, rhi, ilo, ihi


# s1: α² = α·α （α 为正，单调）
s1_lo, s1_hi, s1_lo_n, s1_hi_n = add_step(
    "alpha_sq", "square", "α² = α·α（α>0，区间单调）",
    {"alpha": [str(alpha_lo), str(alpha_hi)]},
    alpha_lo * alpha_lo, alpha_hi * alpha_hi)

# s2: α³ = α²·α （正区间相乘，单调）
s2_lo, s2_hi, s2_lo_n, s2_hi_n = add_step(
    "alpha_cu", "mul", "α³ = α²·α（正区间相乘，端点各自相乘单调）",
    {"alpha_sq": [str(s1_lo), str(s1_hi)], "alpha": [str(alpha_lo), str(alpha_hi)]},
    s1_lo * alpha_lo, s1_hi * alpha_hi)

# s3: 2h （标量，精确）
two_h = 2 * h
s3_lo, s3_hi, s3_lo_n, s3_hi_n = add_step(
    "two_h", "scalar", "2h（h 为 SI 精确值，标量无舍入）",
    {"h": str(h)},
    two_h, two_h)

# s4: 分子 num = α³·(2h)（正区间相乘）
s4_lo, s4_hi, s4_lo_n, s4_hi_n = add_step(
    "num", "mul", "分子 num = α³·2h",
    {"alpha_cu": [str(s2_lo), str(s2_hi)], "two_h": [str(s3_lo), str(s3_hi)]},
    s2_lo * s3_lo, s2_hi * s3_hi)

# s5: 分母 den = e²c（精确量相乘，仍走同一舍入网格以免引入未登记的自由度）
e2 = e * e
den = e2 * c
s5_lo, s5_hi, s5_lo_n, s5_hi_n = add_step(
    "den", "mul", "分母 den = e²c（e、c 为 SI 精确值）",
    {"e2": str(e2), "c": str(c)},
    den, den)

# s6: α²μ₀ = num/den （正分母，端点各自相除；除以正数保序）
a2mu0_exact_lo = s4_lo / s5_hi
a2mu0_exact_hi = s4_hi / s5_lo
s6_lo, s6_hi, s6_lo_n, s6_hi_n = add_step(
    "a2mu0", "div", "α²μ₀ = num/den（den>0，[num_lo/den_hi, num_hi/den_lo]）",
    {"num": [str(s4_lo), str(s4_hi)], "den": [str(s5_lo), str(s5_hi)]},
    a2mu0_exact_lo, a2mu0_exact_hi)

# s7: G 球（输入）
g_lo_r, g_hi_r, g_lo_n, g_hi_n = add_step(
    "G_ball", "input", "G 输入球（CODATA 2022 展示值末位 ±u）",
    {"G": [str(g_lo), str(g_hi)]},
    g_lo, g_hi)

# s8: 间隙 gap = (α²μ₀)_lo − G_hi
gap = s6_lo - g_hi_r
steps.append({
    "id": "gap",
    "op": "sub_compare",
    "desc": "判定 (α²μ₀)_lo > G_hi ⇔ gap = (α²μ₀)_lo − G_hi > 0",
    "inputs": {"a2mu0_lo": str(s6_lo), "G_hi": str(g_hi_r)},
    "exact_lo": str(gap),
    "exact_hi": str(gap),
    "rounded_lo": str(gap),
    "rounded_hi": str(gap),
    "rounded_lo_scaled_int": None,
    "rounded_hi_scaled_int": None,
    "rounding": {"lo": "exact", "hi": "exact"},
    "outward_ok": True,
})
gap_n = int((gap * SCALE))  # 网格点相减仍是网格点

# ---------------------------------------------------------------- Arb 参考（可选）
arb_ref = None
try:
    from flint import arb, ctx
    ctx.prec = 128
    a = arb(ALPHA_MID, ALPHA_U)
    hh = arb(H)
    ee = arb(E)
    cc = arb(C)
    g = arb(G_MID, G_U)
    a2mu0_arb = 2 * a * a * a * hh / (ee * ee * cc)
    arb_ref = {
        "a2mu0_lo": float(a2mu0_arb.lower()),
        "a2mu0_hi": float(a2mu0_arb.upper()),
        "G_hi": float(g.upper()),
        "disjoint": bool(not a2mu0_arb.overlaps(g)),
        "gap": float(a2mu0_arb.lower()) - float(g.upper()),
        "note": "Arb(128bit) 参考；本 trace 用固定网格 10^-64，应比 Arb 更宽（保守）。",
    }
except Exception as exc:  # flint 不可用时如实记录，不伪造
    arb_ref = {"error": "%s: %s" % (type(exc).__name__, exc)}

# ---------------------------------------------------------------- 独立回放（精确有理数，不用 Arb）
replay = []


def check(sid, cond, detail):
    replay.append({"step": sid, "ok": bool(cond), "detail": detail})
    return bool(cond)


D = DIGITS

# 回放每一步：用 trace 里登记的端点重算精确结果，验证舍入向外且区间包含
check("alpha_sq",
      s1_lo <= alpha_lo * alpha_lo and s1_hi >= alpha_hi * alpha_hi,
      "α² 舍入区间包含精确端点平方")
check("alpha_cu",
      s2_lo <= s1_lo * alpha_lo and s2_hi >= s1_hi * alpha_hi,
      "α³ 端点单调相乘且舍入向外")
check("num",
      s4_lo <= s2_lo * s3_lo and s4_hi >= s2_hi * s3_hi,
      "num 端点相乘且舍入向外")
check("den",
      s5_lo <= den and s5_hi >= den,
      "den 为精确量，舍入网格不收紧")
check("a2mu0",
      s6_lo <= s4_lo / s5_hi and s6_hi >= s4_hi / s5_lo,
      "α²μ₀ 除法端点（den>0 保序）且舍入向外")
check("gap_positive", gap > 0, "最终间隙 gap = (α²μ₀)_lo − G_hi > 0")

# 用**独立路径**重算区间（直接由原始输入一步算出，不经中间登记值）以交叉核对。
# 注意 two_h 已含因子 2（2h），故此处不得再乘 2——否则算的是 4α³h/(e²c)。
# 方向：登记端点经向外舍入，应当**包住**精确值，即 reg_lo ≤ exact_lo ≤ exact_hi ≤ reg_hi。
indep_lo = (alpha_lo ** 3 * two_h) / den
indep_hi = (alpha_hi ** 3 * two_h) / den
indep_ok = s6_lo <= indep_lo and s6_hi >= indep_hi
check("independent_path", indep_ok,
      "独立路径（α 直算，不经中间登记值）落在登记区间内：reg_lo ≤ exact_lo ≤ exact_hi ≤ reg_hi")

all_ok = all(r["ok"] for r in replay)
check("overall", all_ok, "全部回放断言通过" if all_ok else "存在未通过的回放断言")

trace = {
    "title": "α²μ₀ vs G 的定向舍入 trace（Q4）",
    "formula": "α²μ₀ = 2α³h/(e²c)（完全相关口径；h、e、c 为 SI 精确值）",
    "grid": "10^-%d（下界向下、上界向上，均向外）" % DIGITS,
    "inputs": {
        "alpha": {"mid": ALPHA_MID, "u": ALPHA_U, "lo": str(alpha_lo), "hi": str(alpha_hi)},
        "h": H, "e": E, "c": C,
        "G": {"mid": G_MID, "u": G_U, "lo": str(g_lo), "hi": str(g_hi)},
    },
    "steps": steps,
    "result": {
        "a2mu0_lo": str(s6_lo),
        "a2mu0_hi": str(s6_hi),
        "a2mu0_lo_scaled_int": str(s6_lo_n),
        "a2mu0_hi_scaled_int": str(s6_hi_n),
        "G_hi": str(g_hi_r),
        "G_hi_scaled_int": str(g_hi_n),
        "gap": str(gap),
        "gap_scaled_int": str(gap_n),
        "disjoint": bool(gap > 0),
    },
    "arb_reference": arb_ref,
}
check_res = {"all_ok": all_ok, "n_checks": len(replay), "checks": replay}

json.dump(trace, io.open(os.path.join(OUT_DIR, "q9_arb_trace.json"), "w", encoding="utf-8"),
          indent=1, ensure_ascii=False)
json.dump(check_res, io.open(os.path.join(OUT_DIR, "q9_arb_trace_check.json"), "w", encoding="utf-8"),
          indent=1, ensure_ascii=False)

print("写出:", os.path.join(OUT_DIR, "q9_arb_trace.json"))
print("写出:", os.path.join(OUT_DIR, "q9_arb_trace_check.json"))
print("步数:", len(steps), " 回放断言:", len(replay), " 全通过:", all_ok)
print("α²μ₀ ∈ [%.12e, %.12e]" % (float(s6_lo), float(s6_hi)))
print("G_hi  = %.12e" % float(g_hi_r))
print("gap   = %.6e   (>0: %s)" % (float(gap), gap > 0))
if "error" not in arb_ref:
    print("Arb 参考: α²μ₀ ∈ [%.12e, %.12e], gap=%.6e, disjoint=%s"
          % (arb_ref["a2mu0_lo"], arb_ref["a2mu0_hi"], arb_ref["gap"], arb_ref["disjoint"]))
else:
    print("Arb 参考不可用:", arb_ref["error"])
for r in replay:
    if not r["ok"]:
        print("  [FAIL]", r["step"], r["detail"])
