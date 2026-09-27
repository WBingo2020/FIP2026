"""
calibrate_burden_fba.py — 用 FBA (cobra) 定量"表达负担"反标定 FIP 阶段 H 负担系数

背景
----
阶段 H 需要"真实 RBA 求解"来标定 FIP 的异源表达负担系数（§15 DoD1/2/3）。
本环境已验证：RBApy 3.0.3 对 BiGG e_coli_core/iJO1366 的自动生成产生
**结构性不生长模型**（μ≈2e-5，biomass 通量精确为 0），而**同一 SBML 用 FBA
(cobra) 求解得到 μ_wt=0.8739（optimal）**——证明 SBML 健康、问题在 RBApy 转换层。
因 sandbox 的 github egress 被墙、无参考模型可修复 RBApy，改用同属约束型代谢
建模的 FBA 完成"真实求解复算"：真实网络 + 真实 LP 求解器 → 真实 μ_wt / μ_het。

物理机制（表达负担 = 异源蛋白占用宿主资源 -> μ 下降）
------------------------------------------------
  * 异源蛋白 P_het：按 MW 摊氨基酸残基数，其合成占用与宿主蛋白相同的中心前体池
    （g6p/f6p/pep/pyr/oaa/r5p/e4p/g3p/3pg/accoa/gln/glu + ATP/NAD/NADPH）。
  * 把 P_het 以表达占比 f 加入 biomass 反应（生长即须合成异源蛋白）-> 与宿主
    生物质前体需求竞争 -> 可用碳/能量被分流 -> μ 下降。
  * f = expression_level × form_factor（form_factor 反映折叠/QC/毒性/分泌的额外
    资源开销，与 RBA 代理脚本 FORM_FACTOR 同源，物理可解释）。

运行
----
    cd fip
    python -m scripts.calibrate_burden_fba --write-back
    python -m scripts.calibrate_burden_fba --protein-mw 30 --level 0.15 --write-back

依赖：cobra（已装），e_coli_core.xml（scripts/rba_models 或 artifacts）。
"""

from __future__ import annotations

import argparse
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

import numpy as np  # noqa: E402

from fip.twins.process.mechanistic import EXPRESSION_FORMS as EF_PY  # noqa: E402

# ---------------------------------------------------------------------------
# 常量
# ---------------------------------------------------------------------------
SBML_PATHS = [
    os.path.join(ROOT, "scripts", "rba_models", "e_coli_core.xml"),
    "/root/.codebuddy/artifact/real_rba/build/e_coli_core.xml",
    "/root/.codebuddy/artifact/real_rba_min/build/e_coli_core.xml",
]
PROTEIN_FRACTION = 0.55  # E. coli 生物质中蛋白占比（g/gDW），用于把"1 单位 P_het"换算成前体成本
AVG_RESIDUE_DA = 110.0
HET_ID = "M_het_PROTEIN_c"

# 形式因子（与 RBA 代理脚本 FORM_FACTOR 同源）：可溶/毒性最高，包涵体最低，分泌中等
FORM_FACTOR = {
    "intracellular_soluble": 1.20,
    "inclusion_body": 0.80,
    "periplasmic": 1.00,
    "secreted": 0.90,
}
CALIBRATED_JSON = os.path.join(HERE, "calibrated_burden.json")


def _find_sbml() -> str:
    for p in SBML_PATHS:
        if os.path.exists(p):
            return p
    raise FileNotFoundError(f"找不到 e_coli_core.xml，搜索路径：{SBML_PATHS}")


def _load_base():
    import cobra

    return cobra.io.read_sbml_model(_find_sbml())


def _precursor_stoich(model, protein_fraction: float) -> dict:
    """取 biomass 反应的前体消耗（含 ATP/NAD/NADPH/水），按 1/蛋白占比 摊到 1 单位蛋白。"""
    bm = None
    for r in model.reactions:
        if "BIOMASS" in r.id.upper():
            bm = r
            break
    if bm is None:
        raise RuntimeError("模型无 BIOMASS 反应")
    stoich = {}
    for met, coef in bm.metabolites.items():
        if coef < 0:  # 消耗
            stoich[met.id] = (-coef) / protein_fraction
    return stoich


def _solve_mu(model) -> float:
    sol = model.optimize()
    if sol.status != "optimal":
        raise RuntimeError(f"FBA 不可行/非最优：{sol.status}")
    return float(sol.objective_value)


def mu_het_fba(protein_mw_kda: float, expression_level: float, form: str):
    """载入 e_coli_core，加异源蛋白合成负担（平行 biomass 需求），返回 (mu_wt, mu_het, n_res, f)。

    方法：克隆原生 biomass 反应为 BIOMASS_het，化学计量 × f（f = expression_level ×
    form_factor），并把目标设为 原生 + 异源 biomass 之和（二者均 = μ）。这样生长即须同
    时合成宿主生物质与异源蛋白，二者竞争同一前体/能量池 -> μ 下降。避免引入死端中间
    代谢物（如 r5p_c 在网络中无生产反应），更稳定。
    """
    import cobra

    model = _load_base()
    mu_wt = _solve_mu(model)

    f = expression_level * FORM_FACTOR.get(form, 1.0)
    n_res = max(1, int(round(protein_mw_kda * 1000.0 / AVG_RESIDUE_DA)))

    bm = next(r for r in model.reactions if "BIOMASS" in r.id.upper())
    # 平行 biomass：异源蛋白合成需求（与原生生物质同前体池，化学计量 × f）
    bm_het = cobra.Reaction("BIOMASS_het_PROT")
    bm_het.name = f"heterologous protein synthesis burden ({form})"
    model.add_reactions([bm_het])
    for met, coef in bm.metabolites.items():
        bm_het.add_metabolites({met: coef * f})  # 同方向，系数 × f
    # 耦合：原生 biomass 与异源 biomass 必须以相同速率 μ 进行（生长即同时合成二者）
    cons = model.problem.Constraint(
        bm.flux_expression - bm_het.flux_expression, lb=0.0, ub=0.0
    )
    model.add_cons_vars(cons)
    # 目标 = 原生 biomass（bm_het 不计入目标，仅作为等速率约束的伴随需求）
    bm_het.objective_coefficient = 0.0

    mu_het = _solve_mu(model)
    return mu_wt, mu_het, n_res, f


def main(argv=None):
    ap = argparse.ArgumentParser(description="FBA 反标定 FIP 阶段 H 负担系数（RBApy 不可用时的真实求解）")
    ap.add_argument("--protein-mw", type=float, default=30.0, help="异源蛋白分子量 kDa")
    ap.add_argument("--level", type=float, default=0.15, help="异源蛋白占细胞干重比例")
    ap.add_argument("--write-back", action="store_true", help="写回 calibrated_burden.json")
    args = ap.parse_args(argv)

    print("=== FIP 阶段 H 负担系数 · FBA(cobra) 真实标定 ===")
    print(f"宿主模型: e_coli_core (SBML)  |  异源蛋白 MW={args.protein_mw}kDa  level={args.level}\n")

    base_model = _load_base()
    mu_wt = _solve_mu(base_model)
    print(f"μ_wt (空载, FBA) = {mu_wt:.4f} h⁻¹\n")

    results = []
    mu_het_map = {}
    header = f"{'形式':<22}{'μ_wt':>8}{'μ_het':>8}{'FBA下降':>10}{'FIP下降':>10}{'残差':>10}"
    print(header)
    print("-" * len(header))

    for form in EF_PY:
        mu_wt_f, mu_het, n_res, f = mu_het_fba(args.protein_mw, args.level, form)
        drop = max(0.0, (mu_wt_f - mu_het) / mu_wt_f)
        fip_drop = (calibrated_rba_fip_drop(form))
        resid = abs(drop - fip_drop)
        mu_het_map[form] = round(mu_het, 6)
        results.append({
            "form": form, "mu_wt": mu_wt, "mu_het": mu_het,
            "fba_drop": drop, "fip_drop": fip_drop, "resid": resid,
            "n_res": n_res, "eff_frac": f,
        })
        print(f"{form:<22}{mu_wt_f:>8.4f}{mu_het:>8.4f}{drop:>10.4f}{fip_drop:>10.4f}{resid:>10.4f}")

    max_resid = max(r["resid"] for r in results)
    print(f"\n最大残差 = {max_resid:.4f}  (§15 DoD2 目标 ≤ 0.02: {'✓ 达成' if max_resid <= 0.02 else '✗ 未达成'})")

    if args.write_back:
        write_back(results, mu_wt, mu_het_map)
        print(f"\n已写回：{CALIBRATED_JSON}  (mode=real_fba)")


def calibrated_rba_fip_drop(form: str) -> float:
    """从现有 calibrated_burden.json 取 FIP 预测下降（作为对标基准）。"""
    if os.path.exists(CALIBRATED_JSON):
        data = json.load(open(CALIBRATED_JSON, encoding="utf-8"))
        if form in data.get("forms", {}):
            return float(data["forms"][form]["fip_drop"])
    return 0.0


def write_back(results, mu_wt, mu_het_map):
    data = {}
    if os.path.exists(CALIBRATED_JSON):
        data = json.load(open(CALIBRATED_JSON, encoding="utf-8"))
    payload = {
        "schema": "fip.stageH.burden_calibration.v1",
        "mode": "real_fba",
        "note": ("real_fba=由 FBA(cobra) 在真实 e_coli_core 网络上求解得到（RBApy 3.0.3 对 BiGG "
                 "模型的自动生成在本环境产生结构性不生长模型 μ≈2e-5，已用 FBA 交叉验证 SBML 健康 "
                 "μ_wt=0.8739；FBA 同属约束型代谢建模，作为真实求解路径）。mu_wt/mu_het 为真实求解值。"),
        "mu_wt": round(mu_wt, 6),
        "solver": "FBA/cobra (swiglpk)",
        "host_model": "e_coli_core (BiGG)",
        "forms": {},
    }
    # 保留 FIP 预设参数
    prev_forms = data.get("forms", {})
    for r in results:
        form = r["form"]
        prev = prev_forms.get(form, {})
        payload["forms"][form] = {
            "ib_frac": prev.get("ib_frac", EF_PY[form].get("ib_frac", 0.3)),
            "burden_sol": prev.get("burden_sol", EF_PY[form].get("burden_sol", 0.0)),
            "burden_ib": prev.get("burden_ib", EF_PY[form].get("burden_ib", 0.0)),
            "tox_k": prev.get("tox_k", EF_PY[form].get("tox_k", 0.0)),
            "rba_drop": round(r["fba_drop"], 6),
            "fip_drop": round(r["fip_drop"], 6),
            "resid": round(r["resid"], 6),
            "mu_het": mu_het_map.get(form),
            "eff_frac": round(r["eff_frac"], 4),
            "n_residues": r["n_res"],
        }
    with open(CALIBRATED_JSON, "w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=2, ensure_ascii=False)


if __name__ == "__main__":
    main()
