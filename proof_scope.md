# 证明覆盖面声明（proof scope）

配套论文：`G_alpha_squared_mu_zero_paper.tex`（36 页）
声称登记表：`claims.yaml`（C1–C14）
形式化工件：`src/G_alpha_squared_mu_zero.lean`、`src/G_dimension_lattice.lean`
数值工件：`results/*.json`

本文件主动申报"机器检查覆盖到哪、没覆盖到哪、边界条件是什么"。凡未在下文登记的断言，均不得被引用为"已形式化验证"。

## 已验证内容

| 编号 | 内容 | 验证形态 | 强度 |
|---|---|---|---|
| C1 | α²/(c²ε₀) = α²μ₀（定理 1） | Lean `thm1_monoid`，`by rfl` | **形式证明**，`#print axioms` 报告不依赖任何公理 |
| C3 | [G] ≠ [α²μ₀]（定理 3 的量纲部分） | Lean `dimG_ne_dim_mu0` | 形式证明（量纲向量化后为整数向量不等） |
| C7 | 不引入质量标度则 G 不进入任何无量纲单项式 | Lean `no_dimensionless_G` + `dimMono_param` | 形式证明，**限单项式（幂积）类** |
| C9(部分) | Π₁、α 的无量纲性 | Lean `Pi1_dimensionless`、`alpha_dimensionless` | 形式证明 |
| C8(后半) | 给定证书 ⇒ α²μ₀ > G | Lean `cert_gap_pos`、`alpha2mu0_gt_G_of_cert` | 形式证明（证书端点作为整数输入，来历在 Lean 之外） |
| C2,C4,C5,C6,C8,C9 | 各数值断言 | Lean `Float` + Python mpmath(dps=30/50) + Arb 球算术 | **数值/区间验证**，不是证明 |
| C8(补充) | 完全相关口径的第二枚 Arb 证书（μ₀ 由 α 在同一 trace 内导出） | Arb（`q5_correlated_and_spacing.py`）+ `results/q5_correlated_and_spacing.json` | 数值/区间验证；与独立球口径结论一致（宽度差 1.86%，源于 CODATA 对 u(μ₀) 的向上舍入展示） |
| C8(压力) | α 端球口径的三档压力测试（u(α)=1.10/1.132/1.168e-12） | Arb（`q5b_alpha_ball_nist_caliber.py`）+ `results/q5b_alpha_ball_nist_caliber.json` | 三档均判不相交，间距均 115.42 u(G) |
| C10(历史) | 各版 CODATA（2010/2014/2018/2022）输入球下的 α²μ₀ vs G 证书 | Arb（`q6_codata_history.py`）+ `results/q6_codata_history.json` | 四版全部不相交；间距 1.71–1.74e-13 几乎不动，σ 随 G 精度升（21→56→115）；反向阈值：G 中心须移 2616 ppm（=1.74×1998 峰值）才重合 |
| C11(MDL) | Λ≈α⁻¹ 的 MDL/描述长度先验（家族规模惩罚折进单一统计量） | `q7_mdl_prior.py` + `results/q7_mdl_prior.json` | 两段码：先验 ΔL=−2.50 bit（BF10=0.177，落在 O(1) 含糊区）；含选择 ΔL=−20.07 bit（BF10=9.1e-7，强支持事后选择）；B≥6 起即使不计选择也无法先验压缩 |
| C12(GUT-RG) | 单圈 RGE + M_SUSY 阈值扫描：M_X 区间 vs M*c² | `q8_gut_running.py` + `results/q8_gut_running.json` | 数值验证（**仅 1 圈**）；M_X ∈ [4.96,14.10]×10¹⁵ GeV，阈值主导（张度 2.84 倍）远大于测量（0.95%）；M*c²=8.909×10¹⁶ GeV 高出全区间 6.3–17.9 倍，故"与 GUT 能标对应"不成立 |
| C4(补充) | 同一均匀靶零模型的解析精确 p（覆盖测度 = 最近邻被估量的闭式） | 区间合并逐段求长（无 MC，dps=60 复核） | 数值验证；与 MC 在小数点后三位一致（差 0.17 倍标准误） |
| C13 | 无量纲单项式的四条闭包律（乘/逆/自然数幂/整数幂 ⇒ 有理函数闭包）与超越函数类型规则 transTy | Lean `dimensionless_{mul,inv,npow,zpow}_closed`、`dimOf_mul`、`dimOf_inv`、`transTy_spec`、`transTy_rejects_dimensional`、`no_fractional_exponent` | 形式证明，**限"单项式经四则运算与整数幂复合、再经保无量纲超越映射"的有界表达式类** |
| C14 | α²μ₀ 的定向舍入 trace（10⁻⁶⁴ 网格、8 步、含舍入方向标注）+ 独立精确有理数复核器；Lean 侧区间内核与精确整数判据 | `src/q9_arb_trace.py` → `results/q9_arb_trace.json`、`results/q9_arb_trace_check.json`（8/8）；Lean `Ival.mul_contains`、`Ival.sq_contains`、`q4_exact_gap`、`q4_gt_of_bounds`、`q4_trace_endpoints_gap` | trace 为数值/区间验证；`q4_exact_gap`、`q4_trace_endpoints_gap` 为**形式证明且不依赖任何公理** |

Lean 全部定理经 `#print axioms` 逐条核验（2026-10-17 复测，`src/G_dimension_lattice.lean` 末尾打印二十六条），按公理集分组如下，**二十六条均无 `sorryAx`**：

- **不依赖任何公理（4 条）**：`thm1_monoid`；`gauss_mu0_dimensionless`（`rfl`）；`q4_exact_gap` 与 `q4_trace_endpoints_gap`（由 `decide` 在核内直接判定常量整数不等式）。
- **仅 `[propext]`（2 条）**：`transTy_spec`；`transTy_rejects_dimensional`。
- **`[propext, Classical.choice, Quot.sound]`（12 条）**：五条量纲判定链 `dimMono_eq_zero_iff`、`dimMono_param`、`Pi1_dimensionless`、`alpha_dimensionless`、`no_dimensionless_G`；Gaussian 的 `gauss_alpha_dimensionless`；§4.7(五) 的 `dimOf_mul`、`dimOf_inv` 与四条闭包律 `dimensionless_mul_closed`、`dimensionless_inv_closed`、`dimensionless_npow_closed`、`dimensionless_zpow_closed`。
- **`[propext, Quot.sound]`（8 条）**：`dimG_ne_dim_mu0`、`gauss_dimG_ne_dim_alpha2mu0`、`cert_gap_pos`、`alpha2mu0_gt_G_of_cert`、`no_fractional_exponent`，以及 §4.8 的区间内核 `Ival.mul_contains`、`Ival.sq_contains`、`q4_gt_of_bounds`。

## 未覆盖内容与缺口

1. **C8/C14 的端点来历**：Arb 区间由 Lean 之外的库按定向舍入算出，Lean 内只回放"单调比较"（`cert_gap_pos`/`alpha2mu0_gt_G_of_cert`）与"给定 trace 端点 ⇒ 结论"（`q4_trace_endpoints_gap`），**未**重放 Arb 的运算步骤，也**未**导出区间依赖图。C14 补上的是"运算步骤可复核"（机器可读 trace + 独立精确有理数复核器 + Lean 区间内核），不是"Arb 的依赖图已导出"。另：trace 中的除法步以"精确有理端点相除"实现（分母为正故保序），Lean 内核未把该步公理化。
2. **任意表达式类的量纲不变性**：C7 只覆盖单项式；C13 把证据范围扩到"单项式经四则运算与整数幂复合、再经保无量纲超越映射"的**有界**表达式类（四条闭包律 + 超越函数类型规则）。两类障碍已在机器上定位且**未**消除：**分数指数**（如 √α）需 ℚ，整数指数格容纳不下（`no_fractional_exponent` 把 ½∉ℤ 化为 `2n=1` 在 ℤ 上无解）；**超越常数的数值内容**（π、e 的十进制展开）超出 Lean core，交由区间算术承担。故本文**不**声称已机械化"任何形式的可能关系都不存在"。
3. **μ₀ 量纲重诠释的物理可接受性**：属物理论证，不在形式化范围内（定理 3 的"需重新诠释"部分因此只是解析论证）。
4. **CODATA 调整本身**：本文认证的是"给定 CODATA 输入下的推导"，未重新证明 CODATA 的数据合并过程。
5. **定理 4 的家族边界**：穷举与零模型的结论依赖人为指定的家族（整数系数、π 的 ≤3 次、|aᵢ|≤10）与均匀随机靶；未给出与家族无关的先验框架。
6. **原文其余断言**（波动方程等）：本文未审。
7. **C12 的 RGE 阶数与阈值模型**：只跑单圈 β 函数，未含 2 圈修正；SUSY 阈值以单一 M_SUSY 台阶近似，未含逐粒子阈值的 Yukawa/规范贡献。结论只依赖"量级 + 模型依赖性"（M*c² 高于整个 1 圈区间一个量级以上），故对 2 圈修正稳健；但 M_X 的具体数值不应被当作可观测量引用。

## 边界条件与隐式常数

- 单位制：SI（§4.2、§4.8 的一切数值）与 Heaviside–Lorentz 自然单位制（§4.3、§4.5）。二者结论不可跨口径引用：`116σ` 仅在 SI 数值口径成立。
- 常数版本：CODATA 2022。若改用 CODATA 2014（u_r(G)=4.7×10⁻⁵）或 1998（1.5×10⁻³），显著性随 u(G) 线性缩放，见 §4.8 稳健性段。
- 精度设定：Arb 中点精度 128 bit；`str(arb)` 的显示位数低于端点，一切报告值取 `lower()/upper()`。
- 证书整数化：区间端点以 10⁻²² 为单位**向外**取整（floor 下界、ceil 上界），保证包含性；向内取整会使证书失效（本文自查中实际发生过一次）。
- **原文可及性**：被审计的原文（参考文献 [1]）链接已于 2026-10-16 复测为"该内容已被发布者删除"；本文的转述依据 2026 年 9 月获取的可见文本，外部读者可能无法逐字复核——凡依赖原文措辞的判断（如"双重证据"表述）均应在这一前提下阅读。
- **显著性口径**：§4.4 的 p 值随家族边界双向漂移（缩到 |a_i|<=4 时 p=0.011，扩到 (pi,e) 混合基时 p=0.73），因此 p 不是内禀量，本文只报敏感性区间而不取单点结论。间距感知的解析版本给出同一形状（0.0111→0.1270 单调），但它仍是**给定家族内**的量：把“家族本身也是挑的”计入，对 K=5 个口径取最小 p 再 Bonferroni 得 0.055，该数字只用于说明敏感性，不作为频率学 p 值使用。
- **两枚证书的相关性口径**：独立球（α 与 μ₀ 各取各的球，u_r(α²μ₀)=4.6063×10⁻¹⁰）与完全相关（μ₀=2αh/(e²c)，h/e/c 精确，u_r=3u_r(α)=4.5222×10⁻¹⁰）文中均给出；因 α²μ₀ 侧比 G 侧窄五个数量级，两口径对“是否重合”的判定不产生分歧。
- 阈值语义：`0.27%`、`0.01`、`0.002` 等阈值只是"通过性判据"，不构成证据强度；实质结论由 §4.4 零模型与 §4.8 认证区间给出。
- 随机性可复现性：零模型与家族扩展实验固定种子 20261014 / 20261015，重跑逐位一致。
- **输入球口径**：本文一切区间计算使用 CODATA **展示值的末位**（u(α)=1.1e-12、u(μ₀)=2.0e-16、u(G)=1.5e-15），R4/mpmath 交叉核对已与 q3/q5 对齐到同一组球。由 (20) 与 (11) 的两位展示舍入反推，公共未舍入 u_r ∈ [1.552, 1.576)e-10，即本文 u(α)=1.1e-12 比真值**窄 2.9%--4.5%**（方向不保守）；影响已由 q5b 三档压力测试覆盖（判定不变）。
- **展示区间的取整方向**：正文所有区间按所列位数**向外**取整，故展示区间恒包含认证区间；本文自查中曾有三处（6.6917625531、1.73126、0.2594%）为向内取整，已改正。
