import Lake
open Lake DSL

/-- 本仓库不依赖 Mathlib：Lean 侧只做两件形式化工作——定理 1 的代数恒等式
    （`thm1_monoid`）与 §4.7 量纲格的不变量分类（`no_dimensionless_G` 等），
    其余为 `Float` 数值检查。因此也可直接逐个运行：
      lean src/G_alpha_squared_mu_zero.lean
      lean src/G_dimension_lattice.lean
    工具链由同目录的 `lean-toolchain` 固定为 leanprover/lean4:lean-4.33.1。 -/
package galpha_mu0 where
  -- 版本号交由 lake 默认值；Lean 4.33 的 Lake DSL 不接受字符串字面量于此字段。

@[default_target]
lean_lib G_alpha_squared_mu_zero where
  srcDir := "./src"
  roots := #[`G_alpha_squared_mu_zero, `G_dimension_lattice]
