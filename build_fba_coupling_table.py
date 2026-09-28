"""
build_fba_coupling_table.py — 预计算 FBA 两层耦合查表（供独立 HTML / JS 端口消费）

浏览器无法跑 FBA 求解器，故把阶段 FBA 耦合所需的"查表量"预先算好、随 HTML 下发：
  * drop_frac = μ_het / μ_wt   （FBA 表达负担降幅；o2=30 参考、无溢流）
  * y_ac_gg   = 乙酸产率 g/g   （FBA 氧化还原 yield；o2=8 DO 受限，动态仿真溢流时用）
维度：carbon(葡萄糖/甘油) × fba_temp × strain_factor × expression_form × expression_level。

μ 内禀上限 mu_cap(T)、宿主摄取 effective_uptake(carbon,T,sf)、生长温度因子
growth_temperature_factor(T) 为纯数学（cardinal + Q10），由 JS 端口直接移植、运行期计算，
无需进表；本脚本同时把相关常数写进 JSON 的 constants 段，保证 JS/Python 同源。

运行：
    cd fip && python3 -m scripts.build_fba_coupling_table
产物：
    scripts/fba_coupling_table.json
"""
from __future__ import annotations

import json
import os
import sys
from datetime import datetime, timezone

import cobra

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from scripts.calibrate_burden_process import (  # noqa: E402
    CARBON, _load_base, _add_burden, _apply_process,
)
from scripts.calibrate_burden_fba import FORM_FACTOR  # noqa: E402
from scripts.uptake_capacity import (  # noqa: E402
    effective_uptake, mu_cap, growth_temperature_factor,
    MU_REF_37, Q10_MU, T_MIN, T_OPT, T_MAX, T_REF,
)
from scripts.calibrate_burden_pichia import (  # noqa: E402
    metrics_pichia, effective_uptake_pichia, growth_temperature_factor_pichia,
    mu_cap_pichia, MU_REF_P_GLY, MU_REF_P_MEOM, Q10_MU_P, T_REF_P, T_OPT_P,
    T_MAX_P, T_MIN_P, V_LIT_P, Q10_UPTAKE_P, MM_P, YXS_P, FORM_FACTOR_PICHIA,
)
from fip.twins.process.mechanistic import EXPRESSION_FORMS as EF_PY  # noqa: E402

OUT = os.path.join(HERE, "fba_coupling_table.json")

# 查表网格（覆盖 Streamlit / HTML 控件范围）
FBA_TEMP_GRID = [25.0, 28.0, 31.0, 34.0, 37.0, 40.0, 42.0]
STRAIN_GRID = [0.5, 0.75, 1.0, 1.25, 1.5, 2.0]
LEVEL_GRID = [0.05, 0.10, 0.15, 0.20, 0.30]
FORMS = ["none"] + list(EF_PY.keys())            # none = 无表达形式（ff=1.0）
FF_MAP = {"none": 1.0}
FF_MAP.update(FORM_FACTOR)
O2_REF = 30.0       # 参考工况：氧充足、无溢流（算 μ_het/μ_wt 降幅）
O2_LIMITED = 8.0    # DO 受限工况：算乙酸氧化还原 yield
MM = {"glucose": 180.0, "glycerol": 92.0}

# 毕赤酵母网格（Pichia 专属温区，与 ferment_sim.js growthTemperatureFactorPichia 一致）
FBA_TEMP_GRID_P = [20.0, 24.0, 28.0, 30.0, 33.0, 36.0]
STRAIN_GRID_P = STRAIN_GRID
LEVEL_GRID_P = LEVEL_GRID
FORMS_P = FORMS                                    # 表达形式两种宿主共享
FF_MAP_P = FF_MAP
O2_REF_P = 30.0
MM_P_ALL = {"glycerol": 92.0, "methanol": 32.0}

_model_cache = {}


def get_model(carbon: str):
    if carbon not in _model_cache:
        _model_cache[carbon] = _load_base(carbon)
    return _model_cache[carbon]


def metrics_fast(carbon: str, sub_max, o2_max, f, T):
    """复用 FBA Oracle 的化学计量，但缓存模型避免重复读 SBML（提速）。"""
    base = get_model(carbon)
    model = base.copy()
    if f and f > 0:
        _add_burden(model, f)
    _apply_process(model, carbon, sub_max, o2_max)
    cap = mu_cap(T)
    bm = next(r for r in model.reactions if "BIOMASS" in r.id.upper())
    bm.upper_bound = cap
    sol = model.optimize()
    if sol.status != "optimal":
        return {"mu": None, "acetate_yield": 0.0}
    exch = CARBON[carbon]["exch"]
    sub = -sol.fluxes.get(exch, 0.0)
    ac = sol.fluxes.get("EX_ac_e", 0.0)
    y_ac = ac / sub if sub > 1e-9 else 0.0
    return {"mu": float(sol.objective_value), "acetate_yield": float(y_ac)}


def build_ecoli():
    data = {}
    for carbon in ("glucose", "glycerol"):
        print(f"[build:ecoli] carbon={carbon}")
        data[carbon] = []
        for t in FBA_TEMP_GRID:
            temp_row = []
            for sf in STRAIN_GRID:
                sub_max = effective_uptake(carbon, T=t, strain_factor=sf)
                mw = metrics_fast(carbon, sub_max, O2_REF, 0.0, t)   # wt，与形式/水平无关
                mu_wt = mw["mu"]
                strain_row = []
                for form in FORMS:
                    ff = FF_MAP.get(form, 1.0)
                    cell = []
                    for lv in LEVEL_GRID:
                        f = lv * ff
                        mh = metrics_fast(carbon, sub_max, O2_REF, f, t)
                        mh_ac = metrics_fast(carbon, sub_max, O2_LIMITED, f, t)
                        mu_het = mh["mu"]
                        feasible = (mu_wt is not None and mu_het is not None)
                        drop = (mu_het / mu_wt) if (feasible and mu_wt > 0) else 1.0
                        y_ac = mh_ac.get("acetate_yield", 0.0)
                        y_ac_gg = y_ac * 60.0 / MM[carbon]
                        cell.append({"drop": round(float(drop), 4), "yac": round(float(y_ac_gg), 4),
                                     "feasible": bool(feasible)})
                    strain_row.append(cell)
                temp_row.append(strain_row)
            data[carbon].append(temp_row)
    return data


def build_pichia():
    """毕赤酵母查表：data_pichia[carbon][temp][strain][form] = list[level 单元格]。

    关键优化（已实测验证）：平行 biomass 异源负担降幅 drop = μ_het/μ_wt ≡ 1/(1+f)，
    对温度(20–36°C)与菌株因子(0.5–2.0)严格不变（纯化学计量前体分流性质）。
    故仅在参考工况 (T=28°C, sf=1.0, o2=30) 对每 (form, level) 跑一次 FBA 取得 drop，
    再沿 temp×strain 广播；温度可行性由 Pichia cardinal 模型 gtf(T)>0 解析给出
    （不依赖逐格 FBA，既快又物理一致）。甘油/甲醇两碳源 drop 完全相同（负担与碳源无关），
    iPP668 统一校准；甲醇相绝对化学计量得率由动态层文献锚定。"""
    # 参考工况的 drop[form][level]
    ref_T, ref_sf = 28.0, 1.0
    ref_sub = effective_uptake_pichia("glycerol", T=ref_T, strain_factor=ref_sf)
    drop_ref = {}     # (form_idx, level_idx) -> drop
    for fi, form in enumerate(FORMS_P):
        ff = FF_MAP_P.get(form, 1.0)
        for li, lv in enumerate(LEVEL_GRID_P):
            f = lv * ff
            mh = metrics_pichia("glycerol", ref_sub, O2_REF_P, f, ref_T)
            mu_het = mh.get("mu_het")
            mu_wt = mh.get("mu_wt")
            drop = (mu_het / mu_wt) if (mu_wt and mu_wt > 0 and mu_het is not None) else 1.0
            drop_ref[(fi, li)] = round(float(drop), 4)
    # 沿 temp × strain 广播
    data = {}
    for carbon in ("glycerol", "methanol"):
        print(f"[build:pichia] carbon={carbon} (drop 由 iPP668 参考工况 FBA 取得并广播)")
        data[carbon] = []
        for t in FBA_TEMP_GRID_P:
            gtf = growth_temperature_factor_pichia(t)
            temp_row = []
            for sf in STRAIN_GRID_P:
                strain_row = []
                for fi, form in enumerate(FORMS_P):
                    cell = []
                    for li, lv in enumerate(LEVEL_GRID_P):
                        drop = drop_ref[(fi, li)]
                        feasible = gtf > 0
                        cell.append({"drop": drop, "yac": 0.0, "feasible": bool(feasible)})
                    strain_row.append(cell)
                temp_row.append(strain_row)
            data[carbon].append(temp_row)
    return data


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--host", choices=("ecoli", "pichia", "both"), default="both")
    args = ap.parse_args()

    out = {
        "schema": "fip.fba_coupling_table.v1",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "note": ("FBA 两层耦合查表：浏览器无法跑 FBA 求解器，故预计算 drop_frac(μ_het/μ_wt) 与 "
                 "y_ac_gg(乙酸 yield g/g) 供 HTML/JS 端口做'伪耦合'。μ 内禀上限 mu_cap(T)、宿主摄取 "
                 "effective_uptake、生长温度因子为纯数学，由 JS 端口同源移植（见 constants）。与 Python "
                 "mechanistic._fba_couple(_pichia) 逐位一致。\n"
                 "· E. coli：甘油/葡萄糖走 e_coli_core_glycerol.xml / e_coli_core.xml。\n"
                 "· Pichia：iPP668 (Chung2010, GS115) 统一校准两阶段异源负担（≡1/(1+f)，与碳源无关）；"
                 "iPP668 不含甲醇交换，甲醇相绝对化学计量得率由动态层文献锚定。"),
    }
    if args.host in ("ecoli", "both"):
        data_ecoli = build_ecoli()
        out["data"] = data_ecoli
        out["constants"] = {
            "MU_REF_37": MU_REF_37, "Q10_MU": Q10_MU,
            "T_MIN": T_MIN, "T_OPT": T_OPT, "T_MAX": T_MAX, "T_REF": T_REF,
            "V_LIT": {"glucose": 10.5, "glycerol": 2.5},
            "Q10_UPTAKE": {"glucose": 2.0, "glycerol": 1.9},
            "MM": MM,
            "FORM_FACTOR": {**FORM_FACTOR, "none": 1.0},
            "O2_REF": O2_REF, "O2_LIMITED": O2_LIMITED,
        }
        out["grid"] = {"fba_temp": FBA_TEMP_GRID, "strain": STRAIN_GRID, "level": LEVEL_GRID}
        out["meta"] = {
            "carbons": list(data_ecoli.keys()),
            "cells_per_carbon": len(FBA_TEMP_GRID) * len(STRAIN_GRID) * len(FORMS) * len(LEVEL_GRID),
            "infeasible_cells": int(sum(
                0 if data_ecoli[c][ti][si][fi][li].get("feasible", True) else 1
                for c in data_ecoli for ti in range(len(FBA_TEMP_GRID)) for si in range(len(STRAIN_GRID))
                for fi in range(len(FORMS)) for li in range(len(LEVEL_GRID))
            )),
            "infeasible_cause": "低温(25/28°C)×低菌株(0.5)下甘油/葡萄糖摄取不足，FBA 模型连野生型维持代谢都供不起 → 无生长物理边界；drop 退化为 1.0、y_ac 退化为 0.0。",
        }
    if args.host in ("pichia", "both"):
        data_pichia = build_pichia()
        out["data_pichia"] = data_pichia
        out["grid_pichia"] = {"fba_temp": FBA_TEMP_GRID_P, "strain": STRAIN_GRID_P, "level": LEVEL_GRID_P}
        out["constants_pichia"] = {
            "MU_REF_P_GLY": MU_REF_P_GLY, "MU_REF_P_MEOM": MU_REF_P_MEOM,
            "Q10_MU_P": Q10_MU_P, "T_REF_P": T_REF_P, "T_OPT_P": T_OPT_P,
            "T_MAX_P": T_MAX_P, "T_MIN_P": T_MIN_P,
            "V_LIT_P": V_LIT_P, "Q10_UPTAKE_P": Q10_UPTAKE_P,
            "MM_P": MM_P_ALL, "YXS_P": YXS_P,
            "FORM_FACTOR_PICHIA": {**FORM_FACTOR_PICHIA, "none": 1.0},
            "O2_REF_P": O2_REF_P,
            "model": "iPP668 (Chung2010, GS115) — 两阶段负担校准；甲醇相 stoich 文献锚定",
        }
        out["meta_pichia"] = {
            "carbons": list(data_pichia.keys()),
            "cells_per_carbon": len(FBA_TEMP_GRID_P) * len(STRAIN_GRID_P) * len(FORMS_P) * len(LEVEL_GRID_P),
            "infeasible_cells": int(sum(
                0 if data_pichia[c][ti][si][fi][li].get("feasible", True) else 1
                for c in data_pichia for ti in range(len(FBA_TEMP_GRID_P)) for si in range(len(STRAIN_GRID_P))
                for fi in range(len(FORMS_P)) for li in range(len(LEVEL_GRID_P))
            )),
            "infeasible_cause": "T<=15°C 或 T>=37°C（Pichia cardinal 温区外）无生长 → feasible=False，drop 退化为 1.0。",
        }
    out["forms"] = FORMS   # 两宿主共享表达形式表

    with open(OUT, "w", encoding="utf-8") as fh:
        json.dump(out, fh, ensure_ascii=False)
    print(f"[build] written {OUT}  (~{os.path.getsize(OUT)//1024} KB)")


if __name__ == "__main__":
    main()
