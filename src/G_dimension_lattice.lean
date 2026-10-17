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

/-- 量纲取逆 = 指数取负 -/
def D.neg (x : Dim) : Dim :=
  ⟨-x.m, -x.l, -x.t, -x.i⟩

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
-- Q2：单项式模型的闭包性质与超越函数类型规则
-- ============================================================================
-- 外审 Q2 问：除单项式外，哪些表达式类（乘积、幂、有理函数、保持无量纲性的超越映射）在
-- 下一步形式化中是可处理的？表示超越常数/函数会遇到什么障碍？
--
-- 本节的答复分三层：
--   (一) **闭包**：无量纲单项式在乘法、取逆、整数幂下封闭。这三条合起来即"有理函数闭包"
--        （有理函数 = 乘 + 逆，见 dimensionless_mul_closed / dimensionless_inv_closed），
--        也就是下一层可直接接入的表达式类。
--   (二) **超越函数类型规则**：exp/log/sin 等要求**无量纲**实参且输出无量纲。Lean core 无实数
--        与超越函数，故此处置于**量纲层面**：transTy 把"参数必须无量纲、结果无量纲"实现为
--        量纲上的部分函数，transTy_spec 证明它恰好只在零量纲上有定义，
--        transTy_rejects_dimensional 证明它对任何非零量纲一律拒绝。这条规则与 (一) 合成，
--        即得"超越映射把无量纲子幺半群映回自身"。
--   (三) **障碍**（诚实边界）：① 分数指数/√ 需要 ℚ，整数指数格容纳不下——no_fractional_exponent
--        把"1/2 不是整数"机械化为"方程 2n = 1 在 ℤ 上无解"；这正是为什么 α 的 √ 形式
--        必须先在实数量纲层消解，而不能直接进入格模型。② 超越常数的**数值**（π、e 的十进制展开）
--        与量纲无关，本模型只断言它们无量纲、因而不改变任何单项式的量纲；其数值内容不在 Lean core
--        的射程内，须由 Arb/区间算术承担（见第 Q4 部分）。

/-- 以 (G, μ₀, c, ℏ, e, m) 为生成元的单项式指数向量；α 无量纲，不占指数位。 -/
structure MonoD where
  g  : Int
  mu : Int
  c  : Int
  hb : Int
  e  : Int
  m  : Int

def MonoD.mul (x y : MonoD) : MonoD :=
  ⟨x.g + y.g, x.mu + y.mu, x.c + y.c, x.hb + y.hb, x.e + y.e, x.m + y.m⟩

def MonoD.inv (x : MonoD) : MonoD :=
  ⟨-x.g, -x.mu, -x.c, -x.hb, -x.e, -x.m⟩

/-- 自然数幂（结构递归；整数幂见 zpow） -/
def MonoD.npow : Nat → MonoD → MonoD
  | 0, _ => ⟨0, 0, 0, 0, 0, 0⟩
  | k + 1, x => MonoD.mul (MonoD.npow k x) x

/-- 整数幂：非负走 npow，负走逆（逆的闭包已单列）。 -/
def MonoD.zpow (n : Int) (x : MonoD) : MonoD :=
  if 0 ≤ n then MonoD.npow n.toNat x else MonoD.inv (MonoD.npow (-n).toNat x)

/-- 把六个指数喂给 dimMono 得量纲。 -/
def dimOf (x : MonoD) : Dim := dimMono x.g x.mu x.c x.hb x.e x.m

def zeroD : Dim := ⟨0, 0, 0, 0⟩

def Dimless (d : Dim) : Prop := d = zeroD

/-- 乘法的量纲同态性（这是下面所有闭包律的种子）。 -/
theorem dimOf_mul (x y : MonoD) :
    dimOf (MonoD.mul x y) = D.mul (dimOf x) (dimOf y) := by
  simp only [dimOf, MonoD.mul, dimMono, D.mul, D.pow, dG, dMu0, dLight, dHbar, dE, dMass,
    Dim.mk.injEq]
  omega

/-- 取逆的量纲反同态性。 -/
theorem dimOf_inv (x : MonoD) :
    dimOf (MonoD.inv x) = D.neg (dimOf x) := by
  simp only [dimOf, MonoD.inv, dimMono, D.mul, D.pow, D.neg, dG, dMu0, dLight, dHbar, dE,
    dMass, Dim.mk.injEq]
  omega

/-- 闭包律 1：无量纲 × 无量纲 = 无量纲。 -/
theorem dimensionless_mul_closed {x y : MonoD}
    (hx : Dimless (dimOf x)) (hy : Dimless (dimOf y)) :
    Dimless (dimOf (MonoD.mul x y)) := by
  simp only [Dimless, dimOf, MonoD.mul, dimMono, D.mul, D.pow, zeroD, dG, dMu0, dLight, dHbar,
    dE, dMass, Dim.mk.injEq] at hx hy ⊢
  omega

/-- 闭包律 2：无量纲的逆仍无量纲（与闭包律 1 合成即"有理函数闭包"）。 -/
theorem dimensionless_inv_closed {x : MonoD} (hx : Dimless (dimOf x)) :
    Dimless (dimOf (MonoD.inv x)) := by
  simp only [Dimless, dimOf, MonoD.inv, dimMono, D.mul, D.pow, zeroD, dG, dMu0, dLight,
    dHbar, dE, dMass, Dim.mk.injEq] at hx ⊢
  omega

/-- 闭包律 3：无量纲的自然数幂无量纲（对 k 归纳，每一步用闭包律 1）。 -/
theorem dimensionless_npow_closed {x : MonoD} (hx : Dimless (dimOf x)) (k : Nat) :
    Dimless (dimOf (MonoD.npow k x)) := by
  induction k with
  | zero =>
      simp only [MonoD.npow, Dimless, dimOf, dimMono, D.mul, D.pow, zeroD, dG, dMu0, dLight,
        dHbar, dE, dMass, Dim.mk.injEq]
      omega
  | succ k ih =>
      simp only [MonoD.npow]
      exact dimensionless_mul_closed ih hx

/-- 闭包律 4：无量纲的整数幂无量纲（负指数经逆，复用闭包律 2、3）。 -/
theorem dimensionless_zpow_closed {x : MonoD} (hx : Dimless (dimOf x)) (n : Int) :
    Dimless (dimOf (MonoD.zpow n x)) := by
  unfold MonoD.zpow
  split
  · exact dimensionless_npow_closed hx _
  · exact dimensionless_inv_closed (dimensionless_npow_closed hx _)

/-- 超越函数的类型规则（量纲层面）：仅当实参无量纲时有定义，且结果无量纲。
    这里把它实现为量纲上的部分函数——对零量纲返回零量纲，否则返回 none。 -/
def transTy (d : Dim) : Option Dim :=
  if d.m = 0 ∧ d.l = 0 ∧ d.t = 0 ∧ d.i = 0 then some zeroD else none

/-- **类型规则的特化**：exp/log/sin 有定义 ⟺ 实参量纲为零，且输出量纲为零。 -/
theorem transTy_spec (d : Dim) : transTy d = some zeroD ↔ d = zeroD := by
  unfold transTy
  by_cases h : d.m = 0 ∧ d.l = 0 ∧ d.t = 0 ∧ d.i = 0
  · rw [if_pos h]
    constructor
    · intro _
      obtain ⟨hm, hl, ht, hi⟩ := h
      cases d with
      | mk m l t i => simp_all [zeroD]
    · intro _
      rfl
  · rw [if_neg h]
    constructor
    · intro hc
      exact absurd hc (by simp)
    · intro hd
      exact absurd (by rw [hd]; exact ⟨rfl, rfl, rfl, rfl⟩) h

/-- **类型规则的否定形式**：任何非零量纲的实参都不能喂给超越函数。 -/
theorem transTy_rejects_dimensional (d : Dim) (hd : d ≠ zeroD) :
    transTy d = none := by
  unfold transTy
  by_cases h : d.m = 0 ∧ d.l = 0 ∧ d.t = 0 ∧ d.i = 0
  · rw [if_pos h]
    exfalso
    apply hd
    obtain ⟨hm, hl, ht, hi⟩ := h
    cases d with
    | mk m l t i => simp_all [zeroD]
  · rw [if_neg h]

/-- **障碍的机械化**：分数指数落在 ℤ 格之外——方程 2n = 1 在 ℤ 上无解。
    这解释了为什么 √α 之类必须先搬到实数量纲层，而不能在整数指数格内表示。 -/
theorem no_fractional_exponent : ¬ (∃ n : Int, 2 * n = 1) := by
  intro h
  obtain ⟨n, hn⟩ := h
  omega

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

-- ============================================================================
-- Q4：定向舍入 trace 的 Lean 侧回放内核
-- ============================================================================
-- 外审 Q4 问：能否导出带定向舍入标注的机器可读 trace，并在 Lean 内用一个小型可信内核
-- （或外部复核器）验证它？此处采取"外部复核器 + Lean 内核"双轨：
--
--   • 外部复核器：src/q9_arb_trace.py 用精确有理数（fractions.Fraction，无浮点、无 Arb）
--     把 α²μ₀ = 2α³h/(e²c) 的区间传播拆成 8 步，每步把上下界分别向外舍入到 10⁻⁶⁴ 网格，
--     落盘 results/q9_arb_trace.json（每步算子/输入/精确值/舍入方向）与
--     results/q9_arb_trace_check.json（逐条回放断言）。实测 8/8 通过。
--     该脚本还用 python-flint 的 Arb（128 bit）独立算一份参考区间，核对 trace 区间**包含**它。
--
--   • Lean 内核（本节）：先把"定向舍入的区间运算 ⇒ 结论"这一层**规则**机械化——
--     Ival.mul_contains / Ival.sq_contains 证明非负区间的乘、平方是对精确值的**外包围**；
--     再给出**假设更弱、结论更强**的整数形式：把 SI 精确值与 CODATA 展示值代入后，
--     α、h、e、c、G 都写成整数比
--         α_lo = 72973525632 / 10¹³,   h = 662607015 / 10⁴²,
--         e    = 1602176634 / 10²⁸,    c = 299792458,   G_hi = 667445 / 10¹⁶
--     判定式 α²μ₀ > G_hi 等价于纯整数比较
--         2·72973525632³·662607015·10¹⁶ > 667445·10²⁵·1602176634²·299792458
--     两端皆为常量，q4_exact_gap 由 decide 在核内直接判定（不依赖任何公理），
--     因而**连舍入都不需要**——这是比 trace 更强的一条交叉验证。
--
-- 诚实边界：α 的展示值（7.2973525643e-3 ± 1.1e-12）本身来自 CODATA 上游合并，不在本节射程内；
-- 本节证明的是"给定该球，结论成立"，与第 Q7 部分的范围界定一致。

/-- 整数端点区间 [lo, hi]。 -/
structure Ival where
  lo : Int
  hi : Int

def Ival.Contains (x : Ival) (v : Int) : Prop := x.lo ≤ v ∧ v ≤ x.hi

/-- 非负区间（下界 ≥ 0）。 -/
def Ival.nonneg (x : Ival) : Prop := 0 ≤ x.lo

/-- 区间乘：端点各自相乘（仅在非负时保序，由 mul_contains 给出健全性）。 -/
def Ival.mul (x y : Ival) : Ival := ⟨x.lo * y.lo, x.hi * y.hi⟩

/-- 区间平方。 -/
def Ival.sq (x : Ival) : Ival := ⟨x.lo * x.lo, x.hi * x.hi⟩

/-- **乘法健全性**：若 a ∈ x、b ∈ y 且 x、y 非负，则 a·b ∈ x·y。
    这正是 trace 里 alpha_sq / alpha_cu / num 三步所依赖的规则。 -/
theorem Ival.mul_contains {x y : Ival} {a b : Int}
    (hx1 : x.lo ≤ a) (hx2 : a ≤ x.hi)
    (hy1 : y.lo ≤ b) (hy2 : b ≤ y.hi)
    (hx : x.nonneg) (hy : y.nonneg) :
    (Ival.mul x y).Contains (a * b) := by
  unfold Ival.Contains Ival.mul Ival.nonneg at *
  constructor
  · exact Int.mul_le_mul hx1 hy1 (by omega) (by omega)
  · exact Int.mul_le_mul hx2 hy2 (by omega) (by omega)

/-- **平方健全性**（乘法健全性的特例，trace 的 alpha_sq 步）。 -/
theorem Ival.sq_contains {x : Ival} {a : Int}
    (hx1 : x.lo ≤ a) (hx2 : a ≤ x.hi) (hx : x.nonneg) :
    (Ival.sq x).Contains (a * a) := by
  have h := Ival.mul_contains hx1 hx2 hx1 hx2 hx hx
  simpa only [Ival.sq, Ival.mul, Ival.Contains] using h

/-- **精确整数判据**（外审 Q4 的加强版）：代入 SI 精确值与 CODATA 展示值后，
    α²μ₀ > G_hi 变成一条常量整数不等式，由 decide 在核内判定，不依赖任何公理。
    注意 h 的尾数有 9 位，故 h = 662607015/10⁴²（不是 10⁴¹）。 -/
theorem q4_exact_gap :
    (2 * 72973525632 ^ 3 * 662607015 : Int) * 10 ^ 16
      > 667445 * (10 ^ 25 * 1602176634 ^ 2 * 299792458) := by
  decide

/-- 与 Q7 回放引理同型：任何落入该整数区间的一对量仍满足 x > y。 -/
theorem q4_gt_of_bounds (x y : Int)
    (hx : x ≥ (2 * 72973525632 ^ 3 * 662607015 : Int) * 10 ^ 16)
    (hy : y ≤ 667445 * (10 ^ 25 * 1602176634 ^ 2 * 299792458)) : x > y := by
  have h := q4_exact_gap
  omega

/-- **与导出工件直接对接**：trace JSON（results/q9_arb_trace.json）里登记的
    10⁻⁶⁴ 网格整数端点 result.a2mu0_lo_scaled_int 与 result.G_hi_scaled_int 之差为正。
    这使 Lean 内核可对**外部复核器落盘的具体数字**下结论，而非只对公式下结论。 -/
theorem q4_trace_endpoints_gap :
    (669176255311555437642179286715984120030000500084569479 : Int)
      > 667445000000000000000000000000000000000000000000000000 := by
  decide

-- 自检：对附录 D 清单中的全部定理逐条打印公理依赖，不得出现 sorryAx。
-- 实测输出（2026-10-17，Lean 4.33.1 core，无 Mathlib），共 26 条：
--   thm1_monoid / gauss_mu0_dimensionless                 → 不依赖任何公理
--   q4_exact_gap / q4_trace_endpoints_gap                 → 不依赖任何公理（decide 核内判定）
--   transTy_spec / transTy_rejects_dimensional            → [propext]
--   dimG_ne_dim_mu0 / gauss_dimG_ne_dim_alpha2mu0
--     / cert_gap_pos / alpha2mu0_gt_G_of_cert
--     / no_fractional_exponent / Ival.mul_contains
--     / Ival.sq_contains / q4_gt_of_bounds                → [propext, Quot.sound]
--   其余（dimMono_* / Pi1 / alpha / no_dimensionless_G
--     / gauss_alpha_dimensionless / dimOf_mul / dimOf_inv
--     / dimensionless_{mul,inv,npow,zpow}_closed）         → [propext, Classical.choice, Quot.sound]
-- 全部落在白名单 [propext, Classical.choice, Quot.sound] 之内，**无 sorryAx**。
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
-- Q2（闭包与超越函数类型规则）
#print axioms dimOf_mul
#print axioms dimOf_inv
#print axioms dimensionless_mul_closed
#print axioms dimensionless_inv_closed
#print axioms dimensionless_npow_closed
#print axioms dimensionless_zpow_closed
#print axioms transTy_spec
#print axioms transTy_rejects_dimensional
#print axioms no_fractional_exponent
-- Q4（定向舍入 trace 的 Lean 侧回放内核）
#print axioms Ival.mul_contains
#print axioms Ival.sq_contains
#print axioms q4_exact_gap
#print axioms q4_gt_of_bounds
#print axioms q4_trace_endpoints_gap

#print "G_dimension_lattice: 全部定理通过类型检查"
