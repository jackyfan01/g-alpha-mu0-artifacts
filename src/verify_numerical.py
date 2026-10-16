#!/usr/bin/env python3
"""
Stage 1.5: Python 数值重算验证
独立使用 mpmath (dps=30) 高精度重算论文中全部数值断言。
论文: G_alpha_squared_mu_zero_paper.tex
"""

from mpmath import mp, mpf, pi, sqrt, fabs

# Helper to format mpf values
def fmt(val, prec=10):
    """Format mpf value to string with given decimal places."""
    return mp.nstr(val, prec)

def fmt_sci(val, prec=6):
    """Format mpf value in scientific notation."""
    f = float(val)
    return f"{f:.{prec}e}"

def fmt_fix(val, prec=6):
    """Format mpf value in fixed notation."""
    f = float(val)
    return f"{f:.{prec}f}"

# 设置精度
mp.dps = 30

print("=" * 80)
print("Stage 1.5: Python mpmath (dps=30) 独立数值重算验证")
print("=" * 80)

# ============================================================
# 物理常数 (CODATA 2022)
# ============================================================
c      = mpf(299792458)                          # m/s, 精确
mu0    = mpf('1.25663706127e-6')  # N/A^2, CODATA 2022, u_r ~ 1.6e-10 (2019 SI 后非精确)
alpha  = mpf('7.2973525643e-3')  # CODATA 2022
G      = mpf('6.67430e-11')                      # m^3/(kg s^2)
hbar   = mpf('1.054571817e-34')                  # J s

# 导出量
eps0   = 1 / (mu0 * c**2)  # 导出量；u_r(eps0) = u_r(mu0)
M_P    = sqrt(hbar * c / G)                       # kg, 普朗克质量（非约化）
Lambda = 4 * pi**3 + pi**2 + pi                   # 几何常数
M_star = M_P / Lambda                              # kg

print("\n--- 基本常数 ---")
print(f"c       = {c}")
print(f"mu0     = {mu0}")
print(f"eps0    = {eps0}")
print(f"alpha   = {alpha}")
print(f"G       = {G}")
print(f"hbar    = {hbar}")
print(f"M_P     = {M_P}")
print(f"Lambda  = {Lambda}")
print(f"M_star  = {M_star}")

# 用于 GeV 转换: 1 GeV = 1.602176634e-10 J  (1 eV = 1.602176634e-19 J)
GeV_to_J = mpf('1.602176634e-10')  # J/GeV

results = []

# ============================================================
# 定理1: α²/(c²ε₀) = α²μ₀  (代数恒等式检验)
# ============================================================
print("\n" + "=" * 80)
print("定理1: α²/(c²ε₀) = α²μ₀  (代数恒等式)")
print("=" * 80)

lhs_th1 = alpha**2 / (c**2 * eps0)
rhs_th1 = alpha**2 * mu0
diff_th1 = fabs(lhs_th1 - rhs_th1)

print(f"  α²/(c²ε₀) = {lhs_th1}")
print(f"  α²μ₀      = {rhs_th1}")
print(f"  差值      = {diff_th1}")
print(f"  判定: 这是纯代数恒等式 (由 ε₀ = 1/(μ₀c²) 直接推出)")
print(f"  论文也承认此点: '此式为由定义保证的代数恒等式，与物理无关'")

results.append({
    'name': '定理1',
    'claim': 'α²/(c²ε₀) = α²μ₀',
    'type': '代数恒等式',
    'verdict': 'REVIEW - 表演性诚实'
})

# ============================================================
# 定理2: α²μ₀ = 6.6918×10⁻¹¹, 与G偏差0.262%
# ============================================================
print("\n" + "=" * 80)
print("定理2: SI 中 G ≈ α²μ₀ (数值验证)")
print("=" * 80)

alpha2_mu0 = alpha**2 * mu0
G_claimed  = mpf('6.6743e-11')   # 论文声称的G值
G_actual   = mpf('6.67430e-11')  # CODATA 2022

alpha2_mu0_claimed = mpf('6.6918e-11')
dev_pct_claimed    = mpf('0.262')

# 独立计算
dev_abs = fabs(alpha2_mu0 - G_actual)
dev_rel = dev_abs / G_actual * 100  # 百分比

print(f"  α²μ₀ (重算)    = {alpha2_mu0}")
print(f"  α²μ₀ (论文)    = {alpha2_mu0_claimed}")
print(f"  α²μ₀ 偏差      = {fabs(alpha2_mu0 - alpha2_mu0_claimed)}")
print(f"  G (CODATA)     = {G_actual}")
print(f"  G (论文)       = {G_claimed}")
print(f"  |α²μ₀ - G|     = {dev_abs}")
print(f"  相对偏差 (重算) = {dev_rel} %")
print(f"  相对偏差 (论文) = {dev_pct_claimed} %")
print(f"  偏差差值       = {fabs(dev_rel - dev_pct_claimed)} %")

# 论文阈值检查
threshold_th2 = mpf('0.27')  # 论文Lean代码中使用的阈值
print(f"\n  论文判定阈值   = {threshold_th2} %")
print(f"  实际偏差       = {dev_rel} %")
print(f"  阈值 - 偏差    = {threshold_th2 - dev_rel} %")

# 判定
if fabs(alpha2_mu0 - alpha2_mu0_claimed) < mpf('1e-14'):
    th2_val_match = "PASS"
else:
    th2_val_match = "FAIL"

if dev_rel < threshold_th2:
    th2_dev_pass = "PASS"
else:
    th2_dev_pass = "FAIL"

print(f"\n  数值匹配: {th2_val_match}")
print(f"  偏差阈值: {th2_dev_pass}")

results.append({
    'name': '定理2',
    'paper_claim': f'α²μ₀ = {alpha2_mu0_claimed}, 偏差 = {dev_pct_claimed}%',
    'recomputed': f'α²μ₀ = {fmt_sci(alpha2_mu0)}, 偏差 = {fmt_fix(dev_rel,6)}%',
    'verdict': 'PASS' if (th2_val_match == 'PASS' and th2_dev_pass == 'PASS') else 'FAIL'
})

# ============================================================
# 定理4: Λ = 137.0363, α⁻¹ = 137.0360, |Λ-137|=0.0363, |Λ-α⁻¹|=0.0003
# ============================================================
print("\n" + "=" * 80)
print("定理4: Λ ≈ 137 ≈ α⁻¹")
print("=" * 80)

alpha_inv = 1 / alpha

Lambda_claimed   = mpf('137.0363')
alpha_inv_claimed = mpf('137.0360')
diff_137_claimed  = mpf('0.0363')
diff_alpha_claimed = mpf('0.0003')

Lambda_computed = Lambda
alpha_inv_computed = alpha_inv
diff_137 = fabs(Lambda - 137)
diff_alpha = fabs(Lambda - alpha_inv)

print(f"  Λ (重算)       = {Lambda_computed}")
print(f"  Λ (论文)       = {Lambda_claimed}")
print(f"  Λ 偏差         = {fabs(Lambda_computed - Lambda_claimed)}")
print(f"  α⁻¹ (重算)     = {alpha_inv_computed}")
print(f"  α⁻¹ (论文)     = {alpha_inv_claimed}")
print(f"  α⁻¹ 偏差       = {fabs(alpha_inv_computed - alpha_inv_claimed)}")
print(f"  |Λ - 137| (重算) = {diff_137}")
print(f"  |Λ - 137| (论文) = {diff_137_claimed}")
print(f"  |Λ - α⁻¹| (重算) = {diff_alpha}")
print(f"  |Λ - α⁻¹| (论文) = {diff_alpha_claimed}")

# 阈值检查
th4_threshold_1 = mpf('0.1')
th4_threshold_2 = mpf('0.01')

print(f"\n  阈值1: |Λ-137| < {th4_threshold_1}, 实际 = {diff_137}")
print(f"  阈值2: |Λ-α⁻¹| < {th4_threshold_2}, 实际 = {diff_alpha}")

th4_pass_1 = diff_137 < th4_threshold_1
th4_pass_2 = diff_alpha < th4_threshold_2
print(f"  判定1: {'PASS' if th4_pass_1 else 'FAIL'}")
print(f"  判定2: {'PASS' if th4_pass_2 else 'FAIL'}")

results.append({
    'name': '定理4',
    'paper_claim': 'Λ=137.0363, α⁻¹=137.0360, |Λ-137|=0.0363, |Λ-α⁻¹|=0.0003',
    'recomputed': f'Λ={fmt_fix(Lambda_computed,10)}, α⁻¹={fmt_fix(alpha_inv_computed,10)}, |Λ-137|={fmt_fix(diff_137,10)}, |Λ-α⁻¹|={fmt_fix(diff_alpha,10)}',
    'verdict': 'PASS' if (th4_pass_1 and th4_pass_2) else 'FAIL'
})

# ============================================================
# 定理5: M*c² = 8.909×10¹⁶ GeV
# ============================================================
print("\n" + "=" * 80)
print("定理5: M*c² 对应 GUT 能标")
print("=" * 80)

# M*c² in Joules, then convert to GeV
M_star_c2_J = M_star * c**2
M_star_c2_GeV = M_star_c2_J / GeV_to_J

GUT_claimed = mpf('8.909e16')   # GeV, 论文声称
GUT_approx  = mpf('8.9e16')     # GeV, 论文提到的近似值

print(f"  M* (kg)        = {M_star}")
print(f"  M*c² (J)       = {M_star_c2_J}")
print(f"  M*c² (GeV)     = {M_star_c2_GeV}")
print(f"  论文声称 (GeV) = {GUT_claimed}")
print(f"  偏差 (GeV)     = {fabs(M_star_c2_GeV - GUT_claimed)}")

dev_rel_GUT = fabs(M_star_c2_GeV - GUT_claimed) / GUT_claimed * 100
print(f"  相对偏差       = {dev_rel_GUT} %")

# 与 GUT 近似值的偏差
dev_GUT_approx = fabs(M_star_c2_GeV - GUT_approx) / GUT_approx * 100
print(f"  与 8.9e16 的相对偏差 = {dev_GUT_approx} %")

th5_pass = dev_rel_GUT < mpf('1')  # 论文说 < 1%
print(f"\n  论文判定: 相对偏差 < 1%")
print(f"  实际判定: {'PASS' if th5_pass else 'FAIL'}")

results.append({
    'name': '定理5',
    'paper_claim': f'M*c² = {GUT_claimed} GeV',
    'recomputed': f'M*c² = {fmt_sci(M_star_c2_GeV)} GeV',
    'verdict': 'PASS' if th5_pass else 'FAIL'
})

# ============================================================
# 定理6: α²Λ² = 1.000004, G/(α²/M*²) = 0.999996
# ============================================================
print("\n" + "=" * 80)
print("定理6: α²Λ² 与 G/(α²/M*²) 自洽性")
print("=" * 80)

alpha2_Lambda2 = alpha**2 * Lambda**2
ratio = 1 / alpha2_Lambda2  # G / (α²/M*²) = 1/(α²Λ²)

alpha2_Lambda2_claimed = mpf('1.000004')
ratio_claimed = mpf('0.999996')

print(f"  α²Λ² (重算)      = {alpha2_Lambda2}")
print(f"  α²Λ² (论文)      = {alpha2_Lambda2_claimed}")
print(f"  α²Λ² 偏差        = {fabs(alpha2_Lambda2 - alpha2_Lambda2_claimed)}")
print(f"  1/(α²Λ²) (重算)  = {ratio}")
print(f"  1/(α²Λ²) (论文)  = {ratio_claimed}")
print(f"  偏差             = {fabs(ratio - ratio_claimed)}")

# |α²Λ² - 1|
dev_from_1 = fabs(alpha2_Lambda2 - 1)
threshold_th6 = mpf('0.002')

print(f"\n  |α²Λ² - 1| (重算) = {dev_from_1}")
print(f"  论文阈值         = {threshold_th6}")
print(f"  判定: {'PASS' if dev_from_1 < threshold_th6 else 'FAIL'}")

results.append({
    'name': '定理6',
    'paper_claim': f'α²Λ² = {alpha2_Lambda2_claimed}, G/(α²/M*²) = {ratio_claimed}',
    'recomputed': f'α²Λ² = {fmt_fix(alpha2_Lambda2,10)}, G/(α²/M*²) = {fmt_fix(ratio,10)}',
    'verdict': 'PASS' if dev_from_1 < threshold_th6 else 'FAIL'
})

# ============================================================
# 额外检查: 表演性诚实 & ≈偷换为=
# ============================================================
print("\n" + "=" * 80)
print("额外检查: 表演性诚实 & ≈偷换为=")
print("=" * 80)

print("\n--- 检查1: 定理1 是否为表演性诚实 ---")
print("  定理1 (α²/(c²ε₀) = α²μ₀) 是纯代数恒等式。")
print("  由 ε₀ = 1/(μ₀c²) 直接代入即得。")
print("  论文自己也承认: '此式为由定义保证的代数恒等式，与物理无关'")
print("  但论文仍将其列为6个定理之一并'通过Lean验证'。")
print("  判定: 存在表演性诚实 —— 验证了一个无需验证的恒等式，")
print("  给读者造成'经过严格检验'的印象，但实际上是重言式。")

print("\n--- 检查2: ≈ 是否被偷换为 = ---")
print("  论文标题使用 G = α²μ₀ (等号)，但:")
print(f"  实际相对偏差 = {fmt_fix(dev_rel,4)}%")
print("  这不是精确等式，而是近似关系。")
print("  论文摘要中写 'G ≈ α²μ₀' (用approx)，但定理编号和标题中用 '='")
print("  定理2标题: 'G ≈ α²μ₀ (数值验证)' - 这里用了approx，较为诚实")
print("  但定理6标题: 'G = α²/M*² 的自洽性' - 用了等号")
print(f"  实际 α²Λ² = {fmt_fix(alpha2_Lambda2,10)} ≠ 1 (精确)")
print("  判定: 定理6存在≈偷换为=的倾向，但论文正文中有'高度自洽'等")
print("  限定词，且Lean验证使用了阈值判定，不完全算偷换。")

print("\n--- 检查3: 阈值设置是否合理 ---")
print(f"  定理2: 偏差 = {fmt_fix(dev_rel,6)}%, 阈值 = 0.27%")
print(f"  阈值余量 = {fmt_fix(mpf('0.27') - dev_rel,6)} 个百分点")
print("  问题: 0.27% 阈值仅比实际偏差大 ~0.008 个百分点。")
print("  这是一个'恰好通过'的阈值，如果常数精度稍有变化就可能不通过。")
print("  判定: 阈值设置偏紧，有'恰好让结果通过'的嫌疑。")

print(f"\n  定理4: |Λ-137| = {fmt_fix(diff_137,6)}, 阈值 = 0.1")
print(f"  阈值余量 = {fmt_fix(mpf('0.1') - diff_137,6)}")
print("  这个阈值相当宽松，0.0363 远小于 0.1。合理。")

print(f"\n  定理4: |Λ-α⁻¹| = {fmt_fix(diff_alpha,6)}, 阈值 = 0.01")
print(f"  阈值余量 = {fmt_fix(mpf('0.01') - diff_alpha,6)}")
print("  阈值也较宽松。合理。")

print(f"\n  定理6: |α²Λ²-1| = {fmt_sci(dev_from_1)}, 阈值 = 0.002")
print("  阈值极其宽松（实际偏差 ~4e-6，阈值 2e-3）。合理。")

# ============================================================
# 汇总
# ============================================================
print("\n" + "=" * 80)
print("汇总")
print("=" * 80)

for r in results:
    print(f"\n  {r['name']}:")
    print(f"    论文声称: {r.get('paper_claim', r.get('claim', 'N/A'))}")
    print(f"    独立重算: {r.get('recomputed', 'N/A')}")
    print(f"    判定: {r['verdict']}")

print("\n" + "=" * 80)
print("最终结论")
print("=" * 80)
print("""
1. 数值重算结果: 全部数值断言经独立重算验证，数值本身正确。
   - 定理2: α²μ₀ = 6.6918e-11 正确，偏差 0.2616% 正确
   - 定理4: Λ = 137.0363 正确，各差值正确
   - 定理5: M*c² = 8.909e16 GeV 正确
   - 定理6: α²Λ² = 1.000004 正确

2. 表演性诚实: 存在。
   - 定理1 (α²/(c²ε₀) = α²μ₀) 是纯代数恒等式/重言式，
     不需要任何数值验证。将其列为"定理"并通过"Lean验证"
     制造了"经过严格检验"的假象。
   - 论文自己也承认这是恒等式，但仍占一个定理编号。

3. ≈偷换为=: 部分存在。
   - 论文标题 G = α²μ₀ 使用等号，但实际是近似 (偏差0.26%)。
   - 定理6标题 "G = α²/M*² 的自洽性" 使用等号，
     但实际 α²Λ² = 1.000004 ≠ 1。
   - 不过论文正文中使用了"approx"和"高度自洽"等限定词，
     不算严重偷换，但标题有误导性。

4. 阈值合理性: 定理2的0.27%阈值可疑。
   - 实际偏差 0.2616%，阈值 0.27%，余量仅 0.008 个百分点。
   - 这是一个"恰好让结果通过"的阈值，有 cherry-picking 嫌疑。
   - 其他定理的阈值设置合理。
""")

print("\n验证完成。")
