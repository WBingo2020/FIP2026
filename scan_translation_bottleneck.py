"""
scan_translation_bottleneck.py — (b) 翻译容量瓶颈敏感性扫描

目的：演示在更高表达水平（或更低翻译容量 K_TR）下，μ_het 如何被
"核糖体翻译容量"这一 RBA 定义性约束硬压低（乃至归零）。

机理回顾（见 calibrate_burden_fba_rba.py）：
    RIB_POOL = HOST_TRANS + HET_TRANS = (C_PROT + f)·μ
    RIB_POOL ≤ K_TR·μ            # 核糖体丰度 ∝ μ（密度受限）
  => μ>0 可行  ⇔  (C_PROT + f) ≤ K_TR  ⇔  rib_util ≤ 1
  其中 f = expression_level × FORM_FACTOR[form]，rib_util = (C_PROT+f)/K_TR

物理含义：
  * rib_util<1：核糖体充足，翻译非瓶颈，μ_het 由前体竞争主导（平滑下降）；
  * rib_util=1：临界点，核糖体被宿主+异源蛋白合成 100% 占满；
  * rib_util>1：核糖体供不应求，任意 μ>0 均不可行 → μ_het→0（翻译饥饿，
    细胞无法在维持任何生长的同时合成异源蛋白）。

本脚本做两维扫描：
  (A) 表达水平扫描：固定 K_TR=0.73，扫 expression_level，画 4 种形式 μ_het 曲线，
      标注各自的 rib_util=1 阈值（level* = (K_TR-C_PROT)/FF[form]）。
  (B) K_TR 扫描：固定可溶性形式 + level=0.30（默认 K_TR 下已超阈值），扫 K_TR，
      展示 μ_het 在 K_TR=C_PROT+f 处的硬跳变（0 → 前体竞争值）。

运行：
    cd fip && python -m scripts.scan_translation_bottleneck
"""
from __future__ import annotations

import json
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import cobra
from scripts.calibrate_burden_fba import FORM_FACTOR, _load_base
from fip.twins.process.mechanistic import EXPRESSION_FORMS as EF_PY

C_PROT = 0.55
K_TR = 0.73
AVG_RESIDUE_DA = 110.0
HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, "..", "docs", "scan_translation_bottleneck.png")


def _build_het_model(base, f: float, k_tr: float, c_prot: float = C_PROT):
    """在 base 模型上克隆并搭建平行-biomass + 核糖体容量结构，返回求解对象。"""
    model = base.copy()
    bm = next(r for r in model.reactions if "BIOMASS" in r.id.upper())
    mu_wt = float(model.optimize().objective_value)

    # ① 前体池竞争：平行 biomass
    bm_het = cobra.Reaction("BIOMASS_het_PROT")
    model.add_reactions([bm_het])
    for met, coef in bm.metabolites.items():
        bm_het.add_metabolites({met: coef * f})
    model.add_cons_vars(model.problem.Constraint(
        bm.flux_expression - bm_het.flux_expression, lb=0.0, ub=0.0))

    # ② 翻译容量：核糖体池（丰度∝μ）
    rib = cobra.Metabolite("RIBOSOME_c", compartment="c")
    model.add_metabolites([rib])
    rib_pool = cobra.Reaction("RIB_POOL")
    rib_pool.add_metabolites({rib: 1.0})
    host_trans = cobra.Reaction("HOST_TRANS")
    host_trans.add_metabolites({rib: -1.0})
    het_trans = cobra.Reaction("HET_TRANS")
    het_trans.add_metabolites({rib: -1.0})
    model.add_reactions([rib_pool, host_trans, het_trans])
    model.add_cons_vars(model.problem.Constraint(
        rib_pool.flux_expression - host_trans.flux_expression - het_trans.flux_expression,
        lb=0.0, ub=0.0))
    return model, bm, bm_het, host_trans, het_trans, rib_pool, mu_wt


def _bisect_mu(model, bm, bm_het, host_trans, het_trans, rib_pool,
               mu_wt, f, k_tr, c_prot=C_PROT):
    """二分法：最大可行 μ（μ 固定为候选值，检查网络 + 翻译容量是否可支撑）。"""
    def _fix(rxn, val):
        if val >= rxn.lower_bound:
            rxn.upper_bound = val
            rxn.lower_bound = val
        else:
            rxn.lower_bound = val
            rxn.upper_bound = val

    lo, hi = 0.0, mu_wt
    for _ in range(28):
        mid = (lo + hi) / 2.0
        _fix(bm, mid)
        _fix(bm_het, mid)
        _fix(host_trans, c_prot * mid)
        _fix(het_trans, f * mid)
        rib_pool.lower_bound = 0.0
        rib_pool.upper_bound = k_tr * mid
        sol = model.optimize()
        if sol.status == "optimal":
            lo = mid
        else:
            hi = mid
    return lo


def main():
    base = _load_base()
    mu_wt = float(base.optimize().objective_value)
    print(f"μ_wt (空载) = {mu_wt:.4f} h⁻¹ | C_PROT={C_PROT} K_TR={K_TR}")

    # ----- (A) 表达水平扫描（固定 K_TR=0.73）-----
    levels = np.round(np.linspace(0.02, 0.60, 30), 3)
    curves = {}
    thresholds = {}
    for form in EF_PY:
        ff = FORM_FACTOR.get(form, 1.0)
        muv = []
        for lv in levels:
            f = lv * ff
            mdl, bm, bm_het, host_trans, het_trans, rib_pool, _ = _build_het_model(base, f, K_TR)
            muv.append(_bisect_mu(model=mdl, bm=bm, bm_het=bm_het, host_trans=host_trans,
                                  het_trans=het_trans, rib_pool=rib_pool, mu_wt=mu_wt,
                                  f=f, k_tr=K_TR))
        curves[form] = np.array(muv)
        # 临界表达水平：lv* 使 (C_PROT + lv*·FF) = K_TR
        thresholds[form] = (K_TR - C_PROT) / ff

    # ----- (B) K_TR 扫描（可溶性形式, level=0.30，默认 K_TR 下已超阈值）-----
    scan_form = "intracellular_soluble"
    scan_level = 0.30
    ff = FORM_FACTOR[scan_form]
    f_scan = scan_level * ff
    ktrs = np.round(np.linspace(0.30, 1.20, 40), 3)
    mu_k = []
    for kt in ktrs:
        mdl, bm, bm_het, host_trans, het_trans, rib_pool, _ = _build_het_model(base, f_scan, kt)
        mu_k.append(_bisect_mu(model=mdl, bm=bm, bm_het=bm_het, host_trans=host_trans,
                               het_trans=het_trans, rib_pool=rib_pool, mu_wt=mu_wt,
                               f=f_scan, k_tr=kt))

    # ---- 图 ----
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(12.5, 5.2))
    colors = {"intracellular_soluble": "#d62728", "inclusion_body": "#1f77b4",
              "periplasmic": "#2ca02c", "secreted": "#9467bd"}
    for i, (form, arr) in enumerate(curves.items()):
        axA.plot(levels, arr, "-o", ms=3, color=colors.get(form, None), label=form)
        lv = thresholds[form]
        axA.axvline(lv, ls=":", lw=1.2, color=colors.get(form, None), alpha=0.7)
        axA.text(lv, 0.04 + i * 0.055, f" {form}\n lv*={lv:.2f}",
                 color=colors.get(form, None), fontsize=7, va="bottom",
                 ha="right" if lv > 0.25 else "left")
    axA.axhline(0, color="k", lw=0.8)
    axA.set_xlabel("expression_level (g/gDW)")
    axA.set_ylabel("mu_het (h^-1)")
    axA.set_title("(A) Expression-level sweep: mu_het vs level\n"
                  "dotted line = each form's rib_util=1 threshold")
    axA.set_ylim(-0.02, mu_wt * 1.05)
    axA.legend(fontsize=7, loc="upper right")
    axA.grid(alpha=0.3)

    axB.plot(ktrs, mu_k, "-o", ms=3, color="#d62728")
    kcrit = C_PROT + f_scan
    axB.axvline(kcrit, ls="--", color="k", lw=1)
    axB.text(kcrit, mu_wt * 0.5, f"  K_TR* = C_PROT+f\n  = {kcrit:.2f}\n  (ribosome limit)",
             fontsize=8, va="center")
    axB.axhspan(0, 0.01, color="#ffcccc", alpha=0.5)
    axB.set_xlabel("translation-capacity coeff K_TR (ribosome abundance / mu)")
    axB.set_ylabel("mu_het (h^-1)")
    axB.set_title(f"(B) K_TR sweep: soluble, level={scan_level}\n"
                  "hard step 0 -> precursor value at K_TR*=C_PROT+f")
    axB.set_ylim(-0.02, mu_wt * 1.05)
    axB.grid(alpha=0.3)

    fig.suptitle("Translation-capacity bottleneck scan: (b) FBA + ribosome resource-demand model\n"
                 "high expression / low K_TR -> mu_het hard-limited by ribosome pool",
                 fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(FIG, dpi=130)
    print(f"\n已保存图：{FIG}")

    # ---- 文本报告 ----
    print("\n=== (A) 表达水平扫描（K_TR=%.2f）===" % K_TR)
    print(f"{'形式':<22}{'FF':>6}{'lv*(rib=1)':>12}{'μ@lv=0.15':>11}{'μ@lv=0.30':>11}{'μ@lv=0.60':>11}")
    for form in EF_PY:
        arr = curves[form]
        def _mu_at(lv):
            i = int(np.argmin(np.abs(levels - lv)))
            return arr[i]
        print(f"{form:<22}{FORM_FACTOR[form]:>6.2f}{thresholds[form]:>12.3f}"
              f"{_mu_at(0.15):>11.4f}{_mu_at(0.30):>11.4f}{_mu_at(0.60):>11.4f}")

    print("\n=== (B) K_TR 扫描（可溶性, level=%.2f, f=%.3f）===" % (scan_level, f_scan))
    print(f"翻译容量临界 K_TR* = C_PROT + f = {kcrit:.3f}")
    print(f"  K_TR=0.50 (<*) → μ_het = {mu_k[int(np.argmin(np.abs(ktrs-0.50)))]:.4f}  (核糖体饥饿, 归零)")
    print(f"  K_TR=0.73 (=默认) → μ_het = {mu_k[int(np.argmin(np.abs(ktrs-0.73)))]:.4f}")
    print(f"  K_TR=0.95 (>*) → μ_het = {mu_k[int(np.argmin(np.abs(ktrs-0.95)))]:.4f}  (前体竞争值)")

    # 写回扫描结果到 JSON（附加段）
    json_path = os.path.join(HERE, "calibrated_burden.json")
    data = json.load(open(json_path, encoding="utf-8")) if os.path.exists(json_path) else {}
    data.setdefault("rba_approx", {})
    data["rba_approx"]["translation_bottleneck_scan"] = {
        "schema": "fip.stageH.burden_calibration.translation_bottleneck_scan.v1",
        "note": ("演示高表达/低K_TR下翻译容量硬瓶颈：(rib_util=(C_PROT+f)/K_TR>1 时任意μ>0"
                 "不可行→μ_het→0)。(A) 表达水平扫描各形式临界 lv*=(K_TR-C_PROT)/FF；(B) K_TR"
                 "扫描在 K_TR*=C_PROT+f 处 μ_het 硬跳变 0→前体竞争值。"),
        "parameters": {"C_PROT": C_PROT, "K_TR": K_TR, "mu_wt": round(mu_wt, 6)},
        "form_threshold_level": {k: round(v, 4) for k, v in thresholds.items()},
        "scan_A": {"levels": [float(x) for x in levels],
                   "mu_het": {k: [float(x) for x in arr] for k, arr in curves.items()}},
        "scan_B": {"form": scan_form, "level": scan_level, "f": round(f_scan, 4),
                   "ktrs": [float(x) for x in ktrs],
                   "mu_het": [float(x) for x in mu_k],
                   "ktr_crit": round(kcrit, 4)},
        "figure": "scan_translation_bottleneck.png",
    }
    with open(json_path, "w", encoding="utf-8") as fh:
        json.dump(data, fh, indent=2, ensure_ascii=False)
    print(f"\n已扩展写回：{json_path}  (新增 translation_bottleneck_scan 段)")


if __name__ == "__main__":
    main()
