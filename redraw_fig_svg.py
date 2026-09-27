"""redraw_fig_svg.py — 用 calibrated_burden.json 缓存的扫描数据重画成 SVG。

目的：原图是 PNG，被仓库 .gitignore 的 `*.png` 规则挡住（仅 `docs/**/*.png` 白名单放行）。
为彻底绕过忽略、并产出可缩放矢量图，这里从 JSON 中读回扫描结果重绘为 SVG。
SVG 不在 .gitignore 任何规则内，GitHub Desktop 必然能识别。

运行：cd fip && python3 -m scripts.redraw_fig_svg
"""
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
JSON_PATH = os.path.join(HERE, "calibrated_burden.json")
DOC_DIR = os.path.join(HERE, "..", "docs")
SVG = os.path.join(DOC_DIR, "scan_translation_bottleneck.svg")


def main():
    data = json.load(open(JSON_PATH, encoding="utf-8"))
    sc = data["rba_approx"]["translation_bottleneck_scan"]
    C_PROT = sc["parameters"]["C_PROT"]
    K_TR = sc["parameters"]["K_TR"]
    mu_wt = sc["parameters"]["mu_wt"]
    thresholds = sc["form_threshold_level"]
    levels = sc["scan_A"]["levels"]
    curves = sc["scan_A"]["mu_het"]
    sb = sc["scan_B"]
    ktrs = sb["ktrs"]
    mu_k = sb["mu_het"]
    kcrit = sb["ktr_crit"]

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
    axB.axvline(kcrit, ls="--", color="k", lw=1)
    axB.text(kcrit, mu_wt * 0.5, f"  K_TR* = C_PROT+f\n  = {kcrit:.2f}\n  (ribosome limit)",
             fontsize=8, va="center")
    axB.axhspan(0, 0.01, color="#ffcccc", alpha=0.5)
    axB.set_xlabel("translation-capacity coeff K_TR (ribosome abundance / mu)")
    axB.set_ylabel("mu_het (h^-1)")
    axB.set_title(f"(B) K_TR sweep: {sb['form']}, level={sb['level']}\n"
                  "hard step 0 -> precursor value at K_TR*=C_PROT+f")
    axB.set_ylim(-0.02, mu_wt * 1.05)
    axB.grid(alpha=0.3)

    fig.suptitle("Translation-capacity bottleneck scan: (b) FBA + ribosome resource-demand model\n"
                 "high expression / low K_TR -> mu_het hard-limited by ribosome pool",
                 fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    os.makedirs(DOC_DIR, exist_ok=True)
    fig.savefig(SVG, dpi=130)
    print(f"已保存 SVG：{SVG}")


if __name__ == "__main__":
    main()
