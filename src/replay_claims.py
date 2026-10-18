# -*- coding: utf-8 -*-
"""Stage 4.8 / K9 复算入口：把 claims.yaml 的数值与 Lean 声称变成可当场复算的判据。

结论行约定（与 proof_deliverables_gate.py 的 K9 契约对齐）：
  stdout 末行必须打印 `VERDICT: PASS` 或 `VERDICT: FAIL <项,项>`；
  退出码 0 而没有结论行，按"缺证据不等于证据"判 FAIL。

三条判据来源，缺一不可：
  · 数值：python-flint 的 arb（定向舍入球算术）从 CODATA 2022 **独立复算**；
    被对照的数字由正则从论文正文抠出（md 与 tex 必须抠到同一个数），
    而不是写在脚本里——脚本里的常量只会证明脚本自己。
  · Lean：真跑 `lake build`，把 `#print axioms` 实测与论文附录 D 的分桶逐条对照。
  · 台账：MC 类不可当场重跑的（194480 项枚举 + 种子），做"论文 ↔ results JSON"一致性，
    并在明细里写明它只证一致、不证复算。

用法:
  python replay_claims.py                # 全量（20 项检查）
  python replay_claims.py --claim C4     # 只跑 C4 对应检查（门控按 claim 复算）
"""
import io
import json
import re
import subprocess
import sys
from pathlib import Path

from flint import arb, ctx

ctx.prec = 128
HERE = Path(__file__).resolve().parent              # .../G_alpha_mu0_repo/src
REPO = HERE.parent                                   # .../G_alpha_mu0_repo
RESULTS = REPO / "results"
TOPIC = REPO.parent                                  # .../papers/g_alpha_mu0
# 论文有两条存放路径：可复现仓库自带的 paper/ 与主题目录的交付三件套。
# 两处都在时以 paper/ 为准（那才是公开仓库里读者能取到的件）；都缺则判红。
PAPER_CANDIDATES = [(REPO / "paper", TOPIC)]


def _pick(name):
    for folder in PAPER_CANDIDATES[0]:
        p = folder / name
        if p.exists():
            return p
    return TOPIC / name


PAPER_MD = _pick("G_alpha_squared_mu_zero_paper.md")
PAPER_TEX = _pick("G_alpha_squared_mu_zero_paper.tex")

# 每条 claim 对应哪些检查前缀（一处笔误不该把全篇判红）。
# 检查 id 一律以 claim 号开头，跳过判定就只看这一条前缀，不易误伤。
CLAIM_CHECKS = {c: [c] for c in ("C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9")}
SELECT = None
CLAIM_ID = None
argv = sys.argv[1:]
if "--claim" in argv:
    CLAIM_ID = argv[argv.index("--claim") + 1]
    if CLAIM_ID not in CLAIM_CHECKS:
        print(f"VERDICT: FAIL 未知 claim {CLAIM_ID}（可选 {sorted(CLAIM_CHECKS)}）")
        raise SystemExit(2)
    SELECT = (CLAIM_ID,)      # 必须是 tuple：str.startswith 不接受 list，传 list 会直接 TypeError

out = {}
fails = []


def check(cid, ok, detail):
    if SELECT is not None and not cid.startswith(SELECT):
        print(f"CHECK {cid} SKIP")
        return True
    out[cid] = {"pass": bool(ok), "detail": detail}
    print(f"CHECK {cid} {'PASS' if ok else 'FAIL'}  {detail}")
    if not ok:
        fails.append(cid)
    return ok


def read(path):
    return io.open(path, encoding="utf-8").read() if path.exists() else None


def paper_number(cid, pattern):
    """从论文正文抠出该声称的数字；md 与 tex 必须抠到完全相同的 (值, 指数)。

    抠不到（锚点失效、论文改版没同步）一律判红——"没找到"绝不能当"通过"。
    """
    rx = re.compile(pattern)
    got = {}
    for tag, path in (("md", PAPER_MD), ("tex", PAPER_TEX)):
        text = read(path)
        if text is None:
            return None, f"{tag} 正文缺失（{path.name}）"
        m = rx.search(text)
        got[tag] = (m.group(1), m.group(2) if m.re.groups >= 2 else None) if m else None
    if got["md"] is None or got["tex"] is None:
        return None, f"锚点未命中 md={got['md']} tex={got['tex']}（pattern={pattern}）"
    if got["md"] != got["tex"]:
        return None, f"md/tex 不同源：md={got['md']} tex={got['tex']}"
    val, exp = got["md"]
    return (val, exp), f"md/tex 同源抠到 {val}" + (f"×10^{exp}" if exp else "")


def in_ball(ball, value, half):
    """论文字面值四舍五入到最后一位 ⇒ 与球中点的距离须 ≤ 半格 + 计算半径。"""
    tol = arb(half) + arb(str(ball.rad()))
    ok = abs(ball.mid() - arb(value)) <= tol
    return ok, f"实算={ball.str(16)} 判定容差=±{half}"


# ── CODATA 2022 输入球（与 q3_arb_certified.py 同一常数源，另行书写以便交叉）────
alpha = arb("7.2973525643e-3", "1.1e-12")
mu0 = arb("1.25663706127e-6", "2.0e-16")
G = arb("6.67430e-11", "1.5e-15")
c = arb("299792458")
h = arb("6.62607015e-34")
hbar = h / (2 * arb.pi())
me = arb("9.1093837139e-31", "2.8e-40")
a2mu0 = alpha * alpha * mu0
Lam = 4 * arb.pi() ** 3 + arb.pi() ** 2 + arb.pi()

# ── 数值声称：抠论文 → 对复算球 ───────────────────────────────────────────────
rel = (a2mu0 - G) / G * 100
got = paper_number("C2", r"相对偏差为\s*\$([\d.]+)\\%")
if got[0] is None:
    check("C2 相对偏差", False, got[1])
else:
    ok, d = in_ball(rel, got[0][0], "5e-5")
    check("C2 相对偏差", ok, f"{got[1]}；{d}")

got = paper_number("C4", r"4\\pi\^3\s*\+\s*\\pi\^2\s*\+\s*\\pi\s*=\s*([\d.]+)")
if got[0] is None:
    check("C4 Lambda 值", False, got[1])
else:
    ok, d = in_ball(Lam, got[0][0], "5e-6")
    check("C4 Lambda 值", ok, f"{got[1]}；{d}")

Mp = (hbar * c / G) ** arb("1/2")
Mstar_GeV = (Mp / Lam) * c * c / arb("1.602176634e-10")
got = paper_number("C5", r"M\^\*c\^2\s*=\s*([\d.]+)\\times10\^\{(\d+)\}")
if got[0] is None:
    check("C5 Mstar", False, got[1])
else:
    lit, exp = got[0]
    ball = Mstar_GeV / arb("1e" + exp)
    ok, d = in_ball(ball, lit, "5e-4")
    check("C5 Mstar", ok, f"{got[1]}；按 1e{exp} 归一后 {d}")

a2L2 = alpha * alpha * Lam * Lam
got = paper_number("C6", r"\\alpha\^2\\Lambda\^2\s*=\s*([\d.]+)\s*\\pm")
if got[0] is None:
    check("C6 alpha2Lambda2", False, got[1])
else:
    ok, d = in_ball(a2L2, got[0][0], "5e-10")
    check("C6 alpha2Lambda2", ok, f"{got[1]}；{d}")
check("C6 alpha2Lambda2 不含 1", bool(not a2L2.contains(arb(1))),
      f"α²Λ²={a2L2.str(16)}，认证区间是否含 1 = {a2L2.contains(arb(1))}")

dlt = Lam - 1 / alpha
got = paper_number("C9", r"\\Lambda\s*-\s*\\alpha\^\{-1\}\\?\s*\|?\s*&?\s*=?\s*([\d.]+)\\times10\^\{(-?\d+)\}")
if got[0] is None:
    check("C9 Lambda减alpha_inv", False, got[1])
else:
    lit, exp = got[0]
    ball = dlt / arb("1e" + exp)
    ok, d = in_ball(ball, lit, "5e-4")
    check("C9 Lambda减alpha_inv", ok, f"{got[1]}；按 1e{exp} 归一后 {d}")
check("C9 Lambda减alpha_inv 不含 0", bool(not dlt.contains(arb(0))),
      f"Λ−α⁻¹={dlt.str(14)}（含 0 = {dlt.contains(arb(0))}）")

# ── C8 §4.8：认证区间不相交（纯复算，不依赖论文文字）+ 论文 σ 口径 ─────────────
gap = float(a2mu0.lower()) - float(G.upper())
sigma_gap = gap / float(G.rad())
dev_sigma = float(abs((a2mu0 - G) / G)) / (float(G.rad()) / float(G.mid()))
check("C8 区间不相交", bool(not a2mu0.overlaps(G)),
      f"gap={gap:.6e} = {sigma_gap:.2f} u(G)；α²μ₀={a2mu0.str(12)}，G={G.str(12)}")
md = read(PAPER_MD) or ""
tex = read(PAPER_TEX) or ""
cm, ct = md.count("116"), tex.count("116")
check("C8 σ=116 与论文口径", cm == ct and cm >= 1 and abs(dev_sigma - 116) <= 0.5,
      f"复算 σ={dev_sigma:.2f}（论文口径 116，容 ±0.5）；「字面 116」md={cm} tex={ct} 次需相等且 ≥1")

# ── C4 附：零模型 p 走台账（MC 不可当场重跑，只证论文与台账一致）───────────────
q2 = RESULTS / "q2_lambda_null_results.json"
try:
    p_ledger = json.loads(q2.read_text(encoding="utf-8"))["A_global_random_target"]["p_value"]
except (OSError, ValueError, KeyError) as exc:
    check("C4 零模型 p 台账", False, f"台账读取失败：{q2.name}（{type(exc).__name__}: {exc}）")
else:
    lit = "0.12708"
    same = md.count(lit) >= 1 and md.count(lit) == tex.count(lit)
    check("C4 零模型 p 台账",
          abs(float(p_ledger) - 0.12708) <= 5e-6 and same,
          f"台账 p={p_ledger}，论文写 {lit}（md={md.count(lit)}/tex={tex.count(lit)} 次）；"
          "此项只证论文与台账一致，194480 项枚举 + 种子 MC 不在此复算")

# ── Lean：真跑 lake build，与论文附录 D 分桶逐条对照 ───────────────────────────
BUILD = {}
if SELECT is None or SELECT and SELECT[0] in ("C1", "C3", "C7"):
    proc = subprocess.run(["lake", "build"], cwd=str(REPO), capture_output=True,
                          text=True, encoding="utf-8", timeout=900)
    BUILD = {"exit": proc.returncode,
             "text": (proc.stdout or "") + (proc.stderr or "")}
    tail = [l for l in BUILD["text"].splitlines() if l.strip()]
    check("C1 lake build 退出码 0", proc.returncode == 0,
          f"exit={proc.returncode}；输出尾：{tail[-1][:80] if tail else '无输出'}")
    check("C3 lake build 退出码 0", proc.returncode == 0, "同上（同一份构建复用）")
    check("C7 lake build 退出码 0", proc.returncode == 0, "同上（同一份构建复用）")


def norm_axioms(s):
    return re.sub(r"\s+", "", s).strip("[]") or "NONE"


def build_axiom_map(text):
    got = {}
    for ln in text.splitlines():
        m = re.search(r"'([A-Za-z][\w.]*)'\s+(depends on axioms:.*|does not depend on any axioms)", ln)
        if m:
            got[m.group(1)] = "NONE" if "does not depend" in m.group(2) else \
                norm_axioms(m.group(2).split(":", 1)[1])
    return got


def paper_axiom_buckets(md_text):
    """附录 D 的分桶行 → {声明: 期望公理集}；只认 build 真报过的名字。"""
    bucket = re.compile(r"^(?:-\s+|\*\s+).*(不依赖任何公理|仅依赖\s*`([^`]+)`|依赖\s*`([^`]+)`)")
    emap = {}
    for ln in md_text.splitlines():
        m = bucket.match(ln)
        if not m:
            continue
        want = "NONE" if m.group(1) == "不依赖任何公理" else norm_axioms(m.group(2) or m.group(3) or "")
        for tok in re.findall(r"`([^`]+)`", ln):
            if "{" in tok or tok.startswith("[") or tok.startswith("\\"):
                continue
            emap[tok] = want
    return emap


LEAN_DECL = {"C1": "thm1_monoid", "C3": "dimG_ne_dim_mu0", "C7": "no_dimensionless_G"}
for cid_exp, decl in LEAN_DECL.items():
    if SELECT is not None and not cid_exp.startswith(SELECT):
        continue
    if "text" not in BUILD:
        check(f"{cid_exp} 公理依赖", False, "本轮未跑 lake build（SELECT 与 build 未同开，属脚本 bug）")
        continue
    got_all = build_axiom_map(BUILD["text"])
    want_all = paper_axiom_buckets(md)
    if decl not in got_all:
        check(f"{cid_exp} 公理依赖", False, f"lake build 输出里没有 {decl} 的 #print axioms 行")
        continue
    if decl not in want_all:
        check(f"{cid_exp} 公理依赖", False, f"附录 D 未给 {decl} 分桶（实测 {got_all[decl]}）")
        continue
    check(f"{cid_exp} 公理依赖", got_all[decl] == want_all[decl],
          f"{decl} 实测={got_all[decl]}；论文附录 D={want_all[decl]}")

if SELECT is None and BUILD:
    got_all = build_axiom_map(BUILD["text"])
    want_all = paper_axiom_buckets(md)
    comparable = {k: v for k, v in want_all.items() if k in got_all}
    bad = sorted(k for k, v in comparable.items() if v != got_all[k])
    check("C0 附录 D 全量对照", bool(comparable) and not bad,
          f"附录 D 列 {len(want_all)} 个声明名，与 build 可比对 {len(comparable)} 个，分歧 {bad or 0}；"
          f"附录 D 有而 build 无={sorted(set(want_all) - set(got_all))}")

# ── 汇总 ─────────────────────────────────────────────────────────────────────
print(f"判定所用论文：{PAPER_MD.parent.name}/{PAPER_MD.name}（tex 同源比对：{PAPER_TEX.name}）")
print(f"\n共 {len(out)} 项检查，失败 {len(fails)} 项：{fails or '无'}")
print("VERDICT: PASS" if not fails else f"VERDICT: FAIL {','.join(fails)}")
raise SystemExit(1 if fails else 0)
