# -*- coding: utf-8 -*-
r"""q6：历史 CODATA 回溯 —— 用**每一版各自公布**的 α、μ₀、G 输入球重跑同一条证书，
检验"α²μ₀ 与 G 的认证区间不相交"这一判定是否只对 CODATA 2022 成立。

两轮外审都点名要这条（R1）：§4.8 里原有的"对 G 未来调整的稳健性"只是**文字**论述，
本脚本把它变成实跑的区间轨迹。

口径要点（这是本脚本最容易被搞错的地方）：
  * CODATA 2010 / 2014 处于**旧 SI**，μ₀ ≡ 4π×10⁻⁷ 是**定义**（精确，零不确定度）。
    此时 α²μ₀ = α²·(4π×10⁻⁷)，不确定度只来自 α（相对 2·u_r(α)）。
    **不能**用 μ₀ = 2αh/(e²c) 去"传递"——那条关系在新 SI 下才是定义级恒等式；
    在旧 SI 下它只是一条用测得的 α 成立的相容性关系，直接用会把 μ₀ 的零不确定度
    人为撑大，属于口径错误。
  * CODATA 2018 / 2022 处于**新 SI**，μ₀ = 2αh/(e²c)（h、e、c 精确），
    于是 α²μ₀ = 2α³h/(e²c)，相对不确定度 3·u_r(α)（完全相关）。
    同时给出"α、μ₀ 各给独立球"的版本作对照，量化相关性的价值。

来源（已逐个核验）：
  * 2010  Mohr/Taylor/Newell, Rev. Mod. Phys. 84, 1527 (2012), arXiv:1203.5425
  * 2014  Mohr/Newell/Taylor, Rev. Mod. Phys. 88, 035009 (2016), arXiv:1507.07956
  * 2018  Tiesinga/Mohr/Newell/Taylor, Rev. Mod. Phys. 93, 025010 (2021)
  * 2022  由 NIST CUU 2022 表（本文 §2 表 tab:constants）给出
  数值原文取自 NIST 归档 ASCII 表
  https://physics.nist.gov/cuu/Constants/ArchiveASCII/allascii_YYYY.txt
  （2010/2014/2018 均可取到；2022 无归档 ASCII，直接用本文 §2 表）
"""
import json
import math
import os

from flint import arb, ctx

ctx.prec = 128

HERE = os.path.dirname(os.path.abspath(__file__))
OUT_PATH = os.path.normpath(os.path.join(HERE, "..", "results",
                                         "q6_codata_history.json"))

H = arb("6.62607015e-34")   # 精确
E = arb("1.602176634e-19")  # 精确
C = arb("299792458")        # 精确

MU0_EXACT = arb("12.56637061435917295385057353312e-7")  # 4π×10⁻⁷，旧 SI 定义值

# 各版：[α 中点, u(α), μ0 中点, u(μ0) 或 None=精确, G 中点, u(G), SI, 文献]
EDITIONS = {
    "2010": {"alpha": ("7.2973525698e-3", "2.4e-12"),
             "mu0": None,                                # 精确 4πe-7
             "G": ("6.67384e-11", "8.0e-15"),
             "si": "old", "ref": "Mohr/Taylor/Newell RMP 84, 1527 (2012)"},
    "2014": {"alpha": ("7.2973525664e-3", "1.7e-12"),
             "mu0": None,
             "G": ("6.67408e-11", "3.1e-15"),
             "si": "old", "ref": "Mohr/Newell/Taylor RMP 88, 035009 (2016)"},
    "2018": {"alpha": ("7.2973525693e-3", "1.1e-12"),
             "mu0": ("1.25663706212e-6", "1.9e-16"),
             "G": ("6.67430e-11", "1.5e-15"),
             "si": "new", "ref": "Tiesinga/Mohr/Newell/Taylor RMP 93, 025010 (2021)"},
    "2022": {"alpha": ("7.2973525643e-3", "1.1e-12"),
             "mu0": ("1.25663706127e-6", "2.0e-16"),
             "G": ("6.67430e-11", "1.5e-15"),
             "si": "new", "ref": "NIST CUU 2022（本文 §2 表）"},
}


def interval_of(x):
    return [float(x.lower()), float(x.upper())]


rows = {}
for ed, spec in EDITIONS.items():
    a = arb(*spec["alpha"])
    G = arb(*spec["G"])
    a2 = a * a

    if spec["mu0"] is None:
        mu0 = MU0_EXACT
        a2mu0 = a2 * mu0
        trace = "alpha2_times_exact_mu0"
        u_r_alpha2mu0 = float((a2mu0.rad() / a2mu0.mid()))
        # 独立球 = 同一式（μ0 精确，无"独立 vs 相关"之分）
        a2mu0_indep = a2mu0
    else:
        mu0 = arb(*spec["mu0"])
        a2mu0 = (2 * a ** 3 * H) / (E * E * C)      # 完全相关 trace
        a2mu0_indep = a2 * mu0                      # α、μ₀ 各给独立球
        trace = "two_alpha_cubed_h_over_e2c"
        u_r_alpha2mu0 = float(a2mu0.rad() / a2mu0.mid())

    lo, hi = interval_of(a2mu0)
    gap = lo - float(G.upper())
    rows[ed] = {
        "si": spec["si"], "ref": spec["ref"],
        "alpha": float(a.mid()), "alpha_inv": float(1 / a.mid()),
        "u_alpha": float(a.rad()), "u_r_alpha": float(a.rad() / a.mid()),
        "mu0": float(mu0.mid()), "u_mu0": (0.0 if spec["mu0"] is None else float(mu0.rad())),
        "mu0_exact": spec["mu0"] is None,
        "G": float(G.mid()), "u_G": float(G.rad()),
        "u_r_G": float(G.rad() / G.mid()),
        "trace": trace,
        "alpha2mu0_interval": [lo, hi],
        "alpha2mu0_mid": float(a2mu0.mid()),
        "u_r_alpha2mu0": u_r_alpha2mu0,
        "G_interval": interval_of(G),
        "disjoint": bool(not a2mu0.overlaps(G)),
        "gap": gap,
        "gap_in_uG": gap / float(G.rad()),
        "sigma": abs((float(a2mu0.mid()) - float(G.mid())) / float(G.rad())),
        # α²μ₀ 中心相对 G 中心的偏移（各自版本内）
        "dev_percent": (float(a2mu0.mid()) - float(G.mid())) / float(G.mid()) * 100.0,
    }
    if spec["mu0"] is not None:
        rows[ed]["alpha2mu0_independent"] = {
            "interval": interval_of(a2mu0_indep),
            "u_r": float(a2mu0_indep.rad() / a2mu0_indep.mid()),
            "u_r_ratio_vs_correlated": float(
                (a2mu0_indep.rad() / a2mu0_indep.mid()) / u_r_alpha2mu0 - 1.0),
        }

out = {
    "note": ("每版各自公布的 α/μ₀/G 输入球下的 α²μ₀ vs G 证书。旧 SI（2010/2014）μ₀ 为定义值，"
             "α²μ₀ 不确定度仅来自 α；新 SI（2018/2022）μ₀ 由 α 传递，用完全相关 trace "
             "α²μ₀ = 2α³h/(e²c)。"),
    "constants_exact": {"h": "6.62607015e-34", "e": "1.602176634e-19", "c": "299792458",
                        "mu0_old_SI_definition": "4*pi*1e-7 = 1.2566370614359172953850573533e-6"},
    "editions": rows,
}

# ---------------- 稳健性：G 要动多少才重合 ----------------
# 判据：要使两区间相交，需 G.upper() >= α²μ₀.lower()（中心上移），
# 或在中心不动时把 u(G) 撑到 (α²μ₀.lo - G.mid)。
robust = {}
for ed, r in rows.items():
    a2lo = r["alpha2mu0_interval"][0]
    Gmid = r["G"]
    required_shift = (a2lo - Gmid) / Gmid               # 中心值相对上移
    required_uG = (a2lo - Gmid)                            # 中心不动时所需半宽
    robust[ed] = {
        "alpha2mu0_lower": a2lo,
        "G_mid": Gmid,
        "max_G_center_shift_percent": required_shift * 100.0,
        "max_G_center_shift_in_uG_current": (a2lo - Gmid) / r["u_G"],
        "required_uG_for_overlap": required_uG,
        "uG_inflation_needed_vs_current": required_uG / r["u_G"],
        "u_r_G_needed_for_overlap": required_uG / Gmid,
    }
out["robustness"] = robust

# ---------------- G 的历史轨迹 ----------------
Ghist = [{"edition": ed, "G": rows[ed]["G"], "u_G": rows[ed]["u_G"],
          "u_r_G": rows[ed]["u_r_G"]} for ed in ("2010", "2014", "2018", "2022")]
out["G_history"] = Ghist
out["G_center_shift"] = {
    "2014_vs_2010_ppm": (rows["2014"]["G"] - rows["2010"]["G"]) / rows["2010"]["G"] * 1e6,
    "2018_vs_2014_ppm": (rows["2018"]["G"] - rows["2014"]["G"]) / rows["2014"]["G"] * 1e6,
    "2022_vs_2018_ppm": (rows["2022"]["G"] - rows["2018"]["G"]) / rows["2018"]["G"] * 1e6,
    "2022_vs_2014_ppm": (rows["2022"]["G"] - rows["2014"]["G"]) / rows["2014"]["G"] * 1e6,
    "2022_vs_2010_ppm": (rows["2022"]["G"] - rows["2010"]["G"]) / rows["2010"]["G"] * 1e6,
}
# 最大历史中心位移（相邻版）与它在"最长跨度"上的累计
out["G_center_shift"]["max_adjacent_abs_ppm"] = max(
    abs(out["G_center_shift"][k]) for k in ("2014_vs_2010_ppm", "2018_vs_2014_ppm"))

# 三档压力测试：把 G 端球独立地替换成"历史最坏"——这是对**未来**调整的外推上界
press = {}
a22 = rows["2022"]
G22 = rows["2022"]["G"]
for name, spec in {
    "G_at_2010_value_with_2010_u": ("6.67384e-11", "8.0e-15"),
    "G_at_2014_value_with_2014_u": ("6.67408e-11", "3.1e-15"),
    "G_2022_center_u_inflated_10x": ("6.67430e-11", "1.5e-14"),
    "G_2022_center_u_inflated_100x": ("6.67430e-11", "1.5e-13"),
    # 直接对接 §4.8 正文提到的"1998 年前后 u_r(G) 冲到约 1500 ppm"
    "G_2022_center_u_at_1998_1500ppm": ("6.67430e-11", repr(G22 * 1.5e-3)),
}.items():
    Gx = arb(*spec)
    a2mu0_22 = (2 * arb(*EDITIONS["2022"]["alpha"]) ** 3 * H) / (E * E * C)
    gap = float(a2mu0_22.lower()) - float(Gx.upper())
    press[name] = {
        "G": float(Gx.mid()), "u_G": float(Gx.rad()),
        "u_r_G": float(Gx.rad() / Gx.mid()),
        "disjoint": bool(not a2mu0_22.overlaps(Gx)),
        "gap": gap,
        "gap_in_uG": gap / float(Gx.rad()) if float(Gx.rad()) else None,
        "sigma": abs((float(a2mu0_22.mid()) - float(Gx.mid())) / float(Gx.rad())),
    }
out["G_future_stress_2022_alpha"] = press

# ---------------- 重合阈值：二分求"中心不动时 u(G) 需多大" ----------------
a2mu0_22 = (2 * arb(*EDITIONS["2022"]["alpha"]) ** 3 * H) / (E * E * C)
lo_edge = float(a2mu0_22.lower())
lo_u, hi_u = 0.0, 1e-10            # G 半宽下界/上界；在 lo_u 处 disjoint，在 hi_u 处重合
for _ in range(200):
    mid_u = 0.5 * (lo_u + hi_u)
    if arb(repr(G22), repr(mid_u)).overlaps(a2mu0_22):
        hi_u = mid_u
    else:
        lo_u = mid_u
thr_u = hi_u
out["overlap_threshold_2022_alpha"] = {
    "u_G_threshold": thr_u,
    "u_r_G_threshold": thr_u / G22,
    "u_r_G_threshold_ppm": thr_u / G22 * 1e6,
    "uG_inflation_vs_CODATA2022": thr_u / rows["2022"]["u_G"],
    "reading": ("中心值钉在 CODATA 2022 的 6.67430e-11 不动时，G 的半宽须膨胀到 %.2e"
                "（u_r=%.0f ppm，即当前 u_r 的 %.0f 倍）才与 α²μ₀ 区间相接。"
                % (thr_u, thr_u / G22 * 1e6, thr_u / rows["2022"]["u_G"])),
}
# 历史最大相邻位移（2014→2010 的 36 ppm）相对该阈值还差多少
out["overlap_threshold_2022_alpha"]["max_hist_adjacent_shift_ppm"] = \
    abs(out["G_center_shift"]["max_adjacent_abs_ppm"])
out["overlap_threshold_2022_alpha"]["ratio_threshold_over_max_hist_shift"] = \
    (thr_u / G22) / abs(out["G_center_shift"]["max_adjacent_abs_ppm"] * 1e-6)

json.dump(out, open(OUT_PATH, "w", encoding="utf-8"), indent=1, ensure_ascii=False)
print("写出:", OUT_PATH)

print("\n=== 各版 CODATA 证书 ===")
hdr = "%-5s %-4s %-14s %-10s %-13s %-11s %-13s %-8s %-9s %-7s"
print(hdr % ("ed", "SI", "alpha_inv", "u_r(alpha)", "u_r(a2mu0)", "G", "u_r(G)", "disjoint",
             "gap", "gap/uG"))
for ed in ("2010", "2014", "2018", "2022"):
    r = rows[ed]
    print(hdr % (ed, r["si"], "%.10f" % r["alpha_inv"], "%.2e" % r["u_r_alpha"],
                 "%.3e" % r["u_r_alpha2mu0"], "%.8e" % r["G"], "%.2e" % r["u_r_G"],
                 r["disjoint"], "%.4e" % r["gap"], "%.2f" % r["gap_in_uG"]))

print("\n=== G 要动多少才重合（以各版 α²μ₀ 下界为准）===")
for ed in ("2010", "2014", "2018", "2022"):
    b = robust[ed]
    print("%-5s G 中心须上移 %+.3f%% (%.1f 个当前 u(G))；或中心不动须把 u(G) 撑大 %.0f× (至 u_r=%.1e)"
          % (ed, b["max_G_center_shift_percent"], b["max_G_center_shift_in_uG_current"],
             b["uG_inflation_needed_vs_current"], b["u_r_G_needed_for_overlap"]))

print("\n=== G 历史中心位移（ppm）===")
for k, v in out["G_center_shift"].items():
    print("  %-22s %+9.2f" % (k, v))

print("\n=== 对 2022 α 端的 G 未来压力测试 ===")
for k, v in press.items():
    print("  %-30s u_r(G)=%.2e disjoint=%-5s gap=%.4e = %.2f u(G)  σ=%.1f"
          % (k, v["u_r_G"], v["disjoint"], v["gap"], v["gap_in_uG"] or 0.0, v["sigma"]))

t = out["overlap_threshold_2022_alpha"]
print("\n=== 重合阈值（G 中心钉住 2022 值）===")
print("  u_r(G) 须达 %.0f ppm（当前 %.1f ppm 的 %.0f 倍）才重合；"
      "历史最大相邻位移 %.1f ppm，差 %.0f 倍"
      % (t["u_r_G_threshold_ppm"], rows["2022"]["u_r_G"] * 1e6,
         t["uG_inflation_vs_CODATA2022"], t["max_hist_adjacent_shift_ppm"],
         t["ratio_threshold_over_max_hist_shift"]))
