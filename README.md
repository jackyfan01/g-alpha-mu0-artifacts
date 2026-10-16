# $G \approx \alpha^2\mu_0$ — Lean 4 辅助验证与可复现工件

论文：**$G \approx \alpha^2\mu_0$ 量纲与数值关系的 Lean 4 辅助验证**
（对塔里芙 (2026)《$G=\alpha^2\mu_0$ 量纲与数值矛盾的完整解决》的多层次检验）

本仓库收录论文的全部验证代码与机器可读结果，使任一数值声称都可独立重跑核对。

## 快速复现

```bash
# 1) Python 环境（二选一）
conda env create -f environment.yml && conda activate galpha-mu0
#   或
pip install -r requirements.txt

# 2) 数值验证（在 src/ 目录下运行，结果 JSON 落在当前目录）
cd src
python verify_numerical.py        # 论文全部数值断言的独立重算（mpmath, dps=30）
python R3_lambda_search.py        # Lambda 全族穷举（194,480 个 pi 多项式）
python q2_lambda_null_model.py    # 显式零模型 + p 值 + 顺序统计量 + 置换退化性
python q3b_family_extension.py    # 家族稳健性（有理函数 / 混合基 / 半整数幂）
python q3_arb_certified.py        # Arb 定向舍入认证区间（不相交 / Lambda != alpha^-1）
python R4_interval_arithmetic.py  # mpmath 区间交叉核对
python q5_correlated_and_spacing.py  # 完全相关口径的第二证书 + 同一零模型的解析精确 p
python q5b_alpha_ball_nist_caliber.py # alpha 端球口径的三档压力测试

# 3) Lean 4 验证（无需 mathlib）
lean src/G_alpha_squared_mu_zero.lean    # 打印 6 项判定，全部 PASS
lean src/G_dimension_lattice.lean        # 量纲格定理 + 证书回放，并打印 axioms 自检
```

## 版本锁定

| 组件 | 版本 | 用途 |
|---|---|---|
| Python | 3.13.7 | 全部脚本（3.9+ 应可运行） |
| mpmath | 1.3.0 | 任意精度重算、`mpi` 区间交叉核对 |
| numpy | 2.4.3 | 零模型 Monte Carlo、最近邻 searchsorted |
| python-flint | 0.9.0（捆绑 Arb 2.x） | §4.8 的定向舍入球算术 |
| Lean | 4.33.1（`lean-toolchain`） | 定理 1 形式证明 + 数值检查 |

随机性说明：`q2_lambda_null_model.py` 固定 `SEED = 20261014`，重跑应逐位复现表中 $p$ 值。

## 物理常数取值

CODATA 2022（Mohr, Newell, Taylor, Tiesinga, *Rev. Mod. Phys.* **97**, 025002 (2025),
arXiv:2409.03787）：

| 常数 | 值 | 相对不确定度 |
|---|---|---|
| $\alpha$ | 7.297 352 5643(11)×10⁻³ | 1.6×10⁻¹⁰（NIST 展示值） |
| $\mu_0$ | 1.256 637 061 27(20)×10⁻⁶ N/A² | 1.6×10⁻¹⁰（NIST 展示值） |
| $\varepsilon_0$ | 8.854 187 8188(14)×10⁻¹² F/m | 1.6×10⁻¹⁰（NIST 展示值） |
| $G$ | 6.674 30(15)×10⁻¹¹ m³ kg⁻¹ s⁻² | 2.2×10⁻⁵ |
| $c,\ h,\ e$ | 精确值（2019 SI 定义） | — |

注意：2019 SI 之后 $\mu_0$、$\varepsilon_0$ **不再是精确值**，其不确定度由 $\alpha$ 传递，严格有 $u_r(\mu_0)=u_r(\varepsilon_0)=u_r(\alpha)$。

**展示口径（重要，且对本文不利的一面也在此登记）**：仓库脚本与论文一律取 CODATA **展示值的末位**（$u(\alpha)=1.1\times10^{-12}$、$u(\mu_0)=2.0\times10^{-16}$、$u(G)=1.5\times10^{-15}$），R4/mpmath 交叉核对已与 q3/q5 用同一组输入球。由 $(20)$ 与 $(11)$ 的两位展示舍入反推，公共未舍入 $u_r\in[1.552,\,1.576)\times10^{-10}$，即本文的 $u(\alpha)=1.1\times10^{-12}$ **比真值窄 2.9\%--4.5\%**（方向不保守）。$u_r(\mu_0)$、$u_r(\varepsilon_0)$ 按各自括号反除为 $1.6\times10^{-10}$，而解析传递值恒等于 $u_r(\alpha)=1.51\times10^{-10}$——这个矛盾只是展示舍入，论文表 1 注 a 有完整说明。压力测试见 `src/q5b_alpha_ball_nist_caliber.py`：$u(\alpha)=1.10/1.132/1.168\times10^{-12}$ 给 $u_r(\alpha^2\mu_0)=4.52/4.66/4.80\times10^{-10}$（宽度 $+0\%/+2.95\%/+6.15\%$），**三档全部判不相交**，间距均为 $1.7312553\times10^{-13}=115.42\,u(G)$。

## 目录结构

```
src/
  G_alpha_squared_mu_zero.lean   Lean 4：定理1形式证明 + 定理2/4/5/6数值检查
  G_dimension_lattice.lean       Lean 4：量纲格不变量分类 + Arb 证书回放（无 Mathlib）
  verify_numerical.py            论文全部数值断言的独立重算
  R3_lambda_search.py            Lambda 全族穷举（定理4）
  q2_lambda_null_model.py        零模型显著性检验（定理 4：事后选择偏差）
  q3b_family_extension.py        家族稳健性：有理函数 / 混合基 (pi,e) / 半整数幂
  q3_arb_certified.py            Arb 认证区间（§4.8 严格判定）
  q5_correlated_and_spacing.py   完全相关口径的第二证书 + 解析精确 p（覆盖测度闭式）
  q5b_alpha_ball_nist_caliber.py alpha 端球口径的三档压力测试
  R4_interval_arithmetic.py      mpmath 区间交叉核对
results/
  lambda_enumeration.json        穷举结果：排名、次优、各阈值计数
  q2_lambda_null_results.json    五种口径的 p 值与顺序统计量
  q3b_family_extension_results.json  四个家族的稳健性比较
  q3_arb_results.json            Arb 区间、不相交判定、sigma 水平
  q5_correlated_and_spacing.json 两版证书（独立球 / 完全相关）+ 五个家族的 p_exact
  q5b_alpha_ball_nist_caliber.json 三档 alpha 端球压力测试（4.52/4.66/4.80e-10，均不相交）
  mpi_crosscheck.json            mpmath 区间交叉核对结果
```

Lean 侧自检：`lean src/G_dimension_lattice.lean` 打印

```
'no_dimensionless_G' depends on axioms: [propext, Classical.choice, Quot.sound]
'alpha2mu0_gt_G_of_cert' depends on axioms: [propext, Quot.sound]
```

即无 `sorryAx`，只依赖 Lean 三条基础公理。

## 主要可复现结论

- 代数恒等式 $\alpha^2/(c^2\varepsilon_0)=\alpha^2\mu_0$：Lean `by rfl` 形式证明成立。
- SI 中 $|\alpha^2\mu_0-G|/G = 0.2616\%$，Arb 认证区间不相交，间距下界 $1.73\times10^{-13}$，
  相当于 $G$ 标准不确定度的 **116σ**。
- $\Lambda = 4\pi^3+\pi^2+\pi = 137.03630$ 是 194,480 个低复杂度 $\pi$ 多项式中的**全局最优**
  （次优差 12 倍），但在"靶值与格点无关"的零模型下 $p \approx 0.13$（各口径 0.123–0.340），
  **统计不显著**。
- 系数空间置换检验在此类"事后挑选"问题上**退化**：$P \approx 1-e^{-1}$ 只反映是否抽中 $\Lambda$。
- $\Lambda \neq \alpha^{-1}$、$\alpha^2\Lambda^2 \neq 1$ 均为 Arb 认证判定（区间不含 0 / 不含 1）。
- 完全相关口径的第二证书（$\mu_0$ 由 $\alpha$ 在同一 trace 内导出）给 $u_r(\alpha^2\mu_0)=4.52\times10^{-10}$；独立球口径为 $4.61\times10^{-10}$，宽 $1.86\%$。**这 1.86\% 不是独立性的代价**：把 $\mu_0$ 端球换成传递值 $1.894\times10^{-16}$ 后两口径只差 $10^{-5}$ 量级（工件 	exttt{independence\_vs\_display}），差额来自 CODATA 把 $u(\mu_0)$ 向上舍入展示成 $(20)$（`q5_correlated_and_spacing.py`）。
- 同一零模型的**解析精确** $p$：均匀随机靶落进第 $i$ 个 Voronoi 胞的概率本就是 $w_i/R$，所以解析闭式 $\lvert\bigcup_i[V_i\pm d]\cap[V_{\min},V_{\max}]\rvert/R$ 与 Monte Carlo 估的是同一个量；五个家族给 0.0111 / 0.0198 / 0.0322 / 0.0698 / 0.1270，与 MC 的 0.12708 在**小数点后三位**一致（差 $9.2\times10^{-5}$，为主口径 MC 标准误 $5.3\times10^{-4}$ 的 $0.17$ 倍）。易错点：对称的 $\min(2d,w_i)/R$ 是错的（本家族有 $17.15\%$ 的相邻间距小于 $2d_{\mathrm{obs}}$，会给出 0.1310），正确式是 $\sum_i(\min(d,L_i)+\min(d,R_i))/R$。

## 许可与引用

代码采用 MIT（`LICENSE-CODE`），论文文本、结果 JSON 与证明工件采用 CC BY 4.0（`LICENSE-DATA`）。引用本仓库时请同时引用论文 `paper/G_alpha_squared_mu_zero_paper.pdf` 与上述 CODATA 2022 出处。


## 代码与数据可用性

- 本仓库是论文全部数值与形式化断言的可复现工件：`src/` 为脚本与 Lean 源，`results/` 为机器可读结果，
  `paper/` 为论文的 `.tex`/`.md`/`.pdf` 三形态。
- 许可证：**代码 MIT**（`LICENSE-CODE`），**论文文本、结果 JSON 与证明工件 CC BY 4.0**（`LICENSE-DATA`）。
- 环境：`environment.yml`（conda）或 `requirements.txt`；Lean 侧无需 Mathlib，工具链由 `lean-toolchain` 固定。
- 复现命令见上节；`results/*.json` 均为脚本直接产出，未手工编辑。

## 被审计的原文（重要说明）

本仓库检验的对象是公众号"塔里芙"(2026) 的文章《G=α²μ₀ 量纲与数值矛盾的完整解决》
（原文链接见 `paper/` 内 `.tex`/`.pdf` 的参考文献 [1]），并顺带检验同期流传的"四大常数关系"方案
（作者"张祥前"，参考文献 [4]）。两点如实记录：

1. 截至 2026-10-16，**原文 [1] 的链接返回"该内容已被发布者删除"**，本仓库的转述依据 2026 年 9 月获取的文本；
2. 本仓库的结论是对**该公式**（G 与 α²μ₀ 的数值关系）的审计，被审计方为具名公开文章；
   结论范围与局限见 `paper/` 与 `proof_scope.md`，请勿在引用时脱离其限定条件。
