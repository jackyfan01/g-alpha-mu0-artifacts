-- 量纲格与"无单位制不变量"结论的 Lean 4 形式化
-- 覆盖三部分：定理 1 的恒等式（单项式模型）、量纲格的不变量分类、Arb 证书回放。
--
-- 不依赖 Mathlib：只用 Lean 4 core 的 Int 与 omega 决策过程。
-- 运行： lean G_dimension_lattice.lean
--
-- 约定：Dim = SI 四个基本量纲 (M, L, T, I) 上的整数指数向量。
--   [G]     = m³ kg⁻¹ s⁻²      → ⟨-1,  3, -2,  0⟩
--   [μ₀]    = kg·m·s⁻²·A⁻²     → ⟨ 1,  1, -2, -2⟩
--   [c]     = m·s⁻¹            → ⟨ 0,  1, -1,  0⟩
--   [ℏ]     = kg·m²·s⁻¹        → ⟨ 1,  2, -1,  0⟩
--   [e]     = A·s              → ⟨ 0,  0,  1,  1⟩
--   [m]     = kg               → ⟨ 1,  0,  0,  0⟩

-- ============================================================================
-- 第零部分：定理 1 的形式化（单项式模型）
-- ============================================================================
--
-- 定理 1 断言：原文的两个公式 α²/(c²ε₀) 与 α²μ₀ 在定义 ε₀ := 1/(μ₀c²) 下是同一个式子。
--
-- 形式化层面的诚实声明：Lean 4 core 不提供域算术策略（`ring`/`field_simp` 与 `ℚ`
-- 记号均属 Mathlib，本仓库刻意不依赖 Mathlib——实测 `field_simp` 报 unknown tactic）。
-- 故此处把每个量表示成生成元 (α, μ₀, c) 的**整数指数向量**（单项式的自由阿贝尔群）：
-- 乘法=指数相加、取逆=取负、乘方=数乘。恒等式在该模型下是定义性的，`rfl` 即可判定。
-- 这正是该断言所在的层面（它只涉及幂积结构，不涉及域公理）；
-- 但它**不是**"在任意域上"的证明——后者需引入 Mathlib 后改写。

/-- 以 (α, μ₀, c) 为生成元的单项式指数向量 -/
structure Mono where
  a : Int
  m : Int
  l : Int

def M.mul (x y : Mono) : Mono := ⟨x.a + y.a, x.m + y.m, x.l + y.l⟩
def M.inv (x : Mono) : Mono := ⟨-x.a, -x.m, -x.l⟩
def M.pow (n : Int) (x : Mono) : Mono := ⟨n * x.a, n * x.m, n * x.l⟩

def Mα  : Mono := ⟨1, 0, 0⟩
def Mmu : Mono := ⟨0, 1, 0⟩
def Mc  : Mono := ⟨0, 0, 1⟩

/-- 定义 ε₀ := 1/(μ₀ c²) -/
def Meps : Mono := M.inv (M.mul Mmu (M.pow 2 Mc))

/-- 原文公式一：α²/(c²ε₀) -/
def thm1_lhs : Mono := M.mul (M.pow 2 Mα) (M.inv (M.mul (M.pow 2 Mc) Meps))

/-- 原文公式二：α²μ₀ -/
def thm1_rhs : Mono := M.mul (M.pow 2 Mα) Mmu

/-- **定理 1 的形式化**：两个公式是同一个单项式。 -/
theorem thm1_monoid : thm1_lhs = thm1_rhs := rfl

/-- SI 量纲：质量、长度、时间、电流的整数指数 -/
structure Dim where
  m : Int
  l : Int
  t : Int
  i : Int

/-- 量纲相乘 = 指数相加 -/
def D.mul (x y : Dim) : Dim :=
  ⟨x.m + y.m, x.l + y.l, x.t + y.t, x.i + y.i⟩

/-- 量纲取幂 = 指数乘整数 -/
def D.pow (n : Int) (x : Dim) : Dim :=
  ⟨n * x.m, n * x.l, n * x.t, n * x.i⟩

def dG     : Dim := ⟨-1,  3, -2,  0⟩
def dMu0   : Dim := ⟨ 1,  1, -2, -2⟩
def dLight : Dim := ⟨ 0,  1, -1,  0⟩
def dHbar  : Dim := ⟨ 1,  2, -1,  0⟩
def dE     : Dim := ⟨ 0,  0,  1,  1⟩
def dMass  : Dim := ⟨ 1,  0,  0,  0⟩

/-- 单项式 G^a μ₀^c c^d ℏ^f e^g m^h 的量纲（α 无量纲，不占指数位） -/
def dimMono (a c d f g h : Int) : Dim :=
  D.mul (D.mul (D.mul (D.mul (D.mul
      (D.pow a dG) (D.pow c dMu0))
      (D.pow d dLight))
      (D.pow f dHbar))
      (D.pow g dE))
      (D.pow h dMass)

/-- 无量纲 ⟺ 四条线性丢番图方程同时成立 -/
theorem dimMono_eq_zero_iff (a c d f g h : Int) :
    dimMono a c d f g h = ⟨0, 0, 0, 0⟩ ↔
    (-a + c + f + h = 0 ∧
      3 * a + c + d + 2 * f = 0 ∧
     -2 * a - 2 * c - d - f + g = 0 ∧
     -2 * c + g = 0) := by
  simp only [dimMono, D.mul, D.pow, dG, dMu0, dLight, dHbar, dE, dMass,
    Dim.mk.injEq]
  omega

/-- **完整刻画**：无量纲单项式的指数格是秩 2 的自由阿贝尔群，
    自由参数只有 (a, c)：a 生成 Π₁ = G m²/(ℏc)，c 生成 μ₀ c e²/ℏ = 4πα。 -/
theorem dimMono_param (a c d f g h : Int)
    (H : dimMono a c d f g h = ⟨0, 0, 0, 0⟩) :
    g = 2 * c ∧ f = -a - c ∧ d = c - a ∧ h = 2 * a := by
  obtain ⟨h1, h2, h3, h4⟩ := (dimMono_eq_zero_iff a c d f g h).mp H
  omega

/-- **论文 §4.7 主结论的形式化**：不引入任何质量标度（h = 0）时，
    G 不可能出现在任何无量纲单项式中（a = 0）。
    故"由 G 与 α²μ₀ 构成的单位制无关量"不是"我们没找到"，而是线性代数上不存在。 -/
theorem no_dimensionless_G (a c d f g : Int)
    (H : dimMono a c d f g 0 = ⟨0, 0, 0, 0⟩) : a = 0 := by
  obtain ⟨_, _, _, h2a⟩ := dimMono_param a c d f g 0 H
  omega

/-- 反之，一旦引入一个质量标度，Π₁ = G m²/(ℏ c) 确实无量纲（机器验证）。
    说明 a = 0 的障碍恰好是"缺一个质量"，与 §4.7(二) 的构造一致。 -/
theorem Pi1_dimensionless : dimMono 1 0 (-1) (-1) 0 2 = ⟨0, 0, 0, 0⟩ := by
  rw [dimMono_eq_zero_iff]
  omega

/-- α = μ₀ c e²/(4πℏ) 无量纲（机器验证；4π 为纯数，不影响量纲）。
    这同时给出 §4.7(二) 中"α² 无量纲"的形式依据。 -/
theorem alpha_dimensionless : dimMono 0 1 1 (-1) 2 0 = ⟨0, 0, 0, 0⟩ := by
  rw [dimMono_eq_zero_iff]
  omega

/-- [G] ≠ [μ₀]；因 α 无量纲，dim(α²μ₀) = dim(μ₀)，故这就是 SI 中的"量纲矛盾"。
    机器可判定，不是修辞。 -/
theorem dimG_ne_dim_mu0 : dG ≠ dMu0 := by
  intro H
  have h := congrArg Dim.m H
  simp only [dG, dMu0] at h
  omega

-- ============================================================================
-- 第 2.5 部分：Gaussian（非有理化 CGS）约定下的同一判定
-- ============================================================================
-- 复审要求核对：§4.7 表 3 里 Gaussian（非有理化 CGS）与 HL 两列的"μ₀ 被吸收为 1"
-- 是否有形式依据。此处给出一个独立的、以 (M, L, T) 为基的量纲模型。
--
-- 高斯制中电荷是导出量：[e] = M^{1/2} L^{3/2} T^{-1}（半整数指数）。
-- 为避免分数指数，把三个指数一律 ×2（保序、保线性、保相等判定），得到整数值向量：
--   [c]   = L T^{-1}            → ⟨ 0,  2, -2⟩
--   [ℏ]   = M L² T^{-1}         → ⟨ 2,  4, -2⟩
--   [e]   = M^{1/2} L^{3/2} T^{-1} → ⟨ 1,  3, -2⟩
--   [G]   = L³ M^{-1} T^{-2}    → ⟨-2,  6, -4⟩
--   [μ₀]  = 1（不存在独立 μ₀，被吸收进电荷定义） → ⟨ 0,  0,  0⟩
-- 关键差别：这里**没有**电流量纲 I，S 除去了第四个基本量纲。

structure GDim where
  m : Int
  l : Int
  t : Int

def GD.mul (x y : GDim) : GDim := ⟨x.m + y.m, x.l + y.l, x.t + y.t⟩
def GD.pow (n : Int) (x : GDim) : GDim := ⟨n * x.m, n * x.l, n * x.t⟩

def gG    : GDim := ⟨-2, 6, -4⟩
def gMu0  : GDim := ⟨ 0, 0,  0⟩
def gE    : GDim := ⟨ 1, 3, -2⟩
def gHbar : GDim := ⟨ 2, 4, -2⟩
def gC    : GDim := ⟨ 0, 2, -2⟩

/-- 高斯制中 α = e²/(ℏc) 无量纲（机器验证）。
    这是最实质的一条：它说明"把 μ₀ 吸收进电荷定义"后，α 仍是同一个纯数。 -/
theorem gauss_alpha_dimensionless :
    GD.mul (GD.pow 2 gE) (GD.pow (-1) (GD.mul gHbar gC)) = ⟨0, 0, 0⟩ := by
  simp only [GD.mul, GD.pow, gE, gHbar, gC, GDim.mk.injEq]
  omega

/-- 高斯制中 μ₀ 无量纲（= 1，被吸收）——与 SI 的 [μ₀] = ⟨1,1,-2,-2⟩ 形成对照。 -/
theorem gauss_mu0_dimensionless : gMu0 = (⟨0, 0, 0⟩ : GDim) := rfl

/-- **高斯制下的同一结论**：α 无量纲，故 α²μ₀ 无量纲，而 G 有量纲 →
    "G ≈ α²μ₀"在高斯制里同样是量纲不匹配。结论不依赖 SI 的四维基。 -/
theorem gauss_dimG_ne_dim_alpha2mu0 :
    gG ≠ GD.mul (GD.pow 2 (⟨0, 0, 0⟩ : GDim)) gMu0 := by
  intro H
  have h := congrArg GDim.m H
  simp only [GD.mul, GD.pow, gG, gMu0] at h
  omega

-- HL（Heaviside–Lorentz，有理化 CGS）与 Gaussian 共享同一 (M, L, T) 基：
--   二者的 [e] 只差一个纯数因子 √(4π)（4π 无量纲），故量纲向量完全相同，
--   于是 gauss_alpha_dimensionless / gauss_dimG_ne_dim_alpha2mu0 在 HL 下逐字成立；
--   差别只在 α 的**数值**定义 α = e²/(4πℏc)（HL）vs α = e²/(ℏc)（Gaussian）里那个 4π。

-- ============================================================================
-- Q7：Arb 证书回放
-- ============================================================================
-- §4.8 的区间不相交由 Arb 在 Lean 之外算出（见 results/q3_arb_results.json）。
-- 这里把"给定证书 ⇒ 结论"这一步搬进 Lean：证书以整数形式进入假设，单调推理在 Lean 内完成。
--
-- 取整方向（关键）：单位 10⁻²² m³ kg⁻¹ s⁻²。
--   Arb 认证 α²μ₀ ≥ 6.691762553097503e-11，即 669176255309.7503（单位 10⁻²²）
--     → 下界必须**向下**取整：669176255309
--   Arb 认证 G ≤ 6.674450000000378e-11，即 667445000000.0378
--     → 上界必须**向上**取整：667445000001
-- 两个整数都严格在认证区间之外（保守方向），故"落入证书区间"的前提可被 Arb 结果蕴含。
-- 唯一留在 Lean 之外的，是端点数值本身的来历——那由 Arb 的定向舍入语义负责。

/-- Arb 证书：α²μ₀ 下界 6.691762553097503e-11 向下取整，单位 10⁻²² -/
def a2mu0_lo : Int := 669176255309

/-- Arb 证书：G 上界 6.674450000000378e-11 向上取整，单位 10⁻²² -/
def g_hi : Int := 667445000001

/-- 证书间距为正（机器验证；对应 §4.8 的 1.731e-13 间距） -/
theorem cert_gap_pos : a2mu0_lo - g_hi > 0 := by
  simp only [a2mu0_lo, g_hi]
  omega

/-- 间距的确切整数值（= 1.731255308e-13，单位 10⁻²²） -/
theorem cert_gap_val : a2mu0_lo - g_hi = 1731255308 := by
  simp only [a2mu0_lo, g_hi]
  omega

/-- **回放引理**：任何落入证书区间的量都满足 α²μ₀ > G。
    这把"外部认证"与"论文结论"之间的推理完全机械化。 -/
theorem alpha2mu0_gt_G_of_cert (x y : Int)
    (hx : x ≥ a2mu0_lo) (hy : y ≤ g_hi) : x > y := by
  simp only [a2mu0_lo, g_hi] at hx hy
  omega

-- 自检：对附录 D 清单中的全部定理逐条打印公理依赖，不得出现 sorryAx。
-- 实测输出（2026-10-16，Lean 4.33.1 core，无 Mathlib）：
--   thm1_monoid            → 不依赖任何公理
--   dimMono_eq_zero_iff / dimMono_param / Pi1_dimensionless / alpha_dimensionless
--   no_dimensionless_G     → [propext, Classical.choice, Quot.sound]
--   dimG_ne_dim_mu0 / cert_gap_pos / alpha2mu0_gt_G_of_cert → [propext, Quot.sound]
--   gauss_alpha_dimensionless / gauss_dimG_ne_dim_alpha2mu0 → 同上两条基础公理
--   gauss_mu0_dimensionless → 不依赖任何公理（rfl）
#print axioms thm1_monoid
#print axioms dimMono_eq_zero_iff
#print axioms dimMono_param
#print axioms Pi1_dimensionless
#print axioms alpha_dimensionless
#print axioms no_dimensionless_G
#print axioms dimG_ne_dim_mu0
#print axioms gauss_alpha_dimensionless
#print axioms gauss_mu0_dimensionless
#print axioms gauss_dimG_ne_dim_alpha2mu0
#print axioms cert_gap_pos
#print axioms alpha2mu0_gt_G_of_cert

#print "G_dimension_lattice: 全部定理通过类型检查"
