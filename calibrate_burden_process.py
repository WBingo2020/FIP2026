"""
calibrate_burden_process.py — 工艺约束下的乙酸溢出模块（阶段 H 补充，支持葡萄糖/甘油二选一）

为什么需要
----------
真实发酵工艺里两类工况会触发 E. coli 乙酸溢出（overflow metabolism）：
  (1) 碳源补料速率 > 宿主碳源摄取能力：胞内碳过剩；
  (2) 溶氧(DO)受限：呼吸链氧化容量被限，丙酮酸来不及全氧化 -> 以乙酸/甲酸/
      乙醇形式排出，同时 rebalance 还原力与 ATP。
标准 FBA 以"最大 biomass 得率"为目标，默认预测乙酸 ≡ 0（见 calibrate_burden_fba.py
的诊断）。本模块把"碳源摄取上限"与"溶氧(DO)氧化容量上限"作为**工艺约束**显式
加入，复现并量化乙酸溢出，并叠加异源表达负担看其如何放大溢出。

FIP 实际工艺：碳源为葡萄糖或甘油二选一
--------------------------------------
e_coli_core 原生只含葡萄糖交换，不含甘油代谢通路。本阶段新增"甘油对应的模型"
`e_coli_core_glycerol.xml`（在 e_coli_core 上补 GLYCt 转运 + GLYK 甘油激酶 +
GLYCK 甘油-3-P 脱氢酶，甘油->DHAP 接入糖酵解），使同一核心模型同时支持葡萄糖与
甘油两种碳源，二者μ可直接对比、且不破坏任何已有葡萄糖校准（葡萄糖-only μ 仍为
1.2917）。"二选一"在模块里通过封死另一碳源交换实现。

机理（已用两个模型实测验证，见脚本内 sweep）
  * 碳源摄取上限 sub_max：EX_<sub>_e LB = -sub_max（宿主摄取能力；补料超此值时稳态
    下细胞只吃到 sub_max，多余碳在胞内累积/溢出）。
  * 溶氧上限 o2_max：EX_o2_e LB = -o2_max（DO 受限的呼吸氧化容量）。
  * 当 o2_max 不足以氧化 sub_max 供来的碳时，FBA 在 biomass 目标下分泌乙酸(+甲酸+
    乙醇) -> 乙酸溢出；μ 被氧化容量限制。
  * 甘油单位碳还原度高于葡萄糖，DO 受限时同样溢出；甘油-only 最大 μ(0.746)低于葡萄糖
    (1.292)，符合其单位碳能量密度的物理预期。
  * 异源负担（平行 biomass，化学计量 ×f，与原生 biomass 等速率耦合）在受限营养/能量
    池里抢资源 -> 进一步压低 μ_het。乙酸通量由 O2/碳氧化还原平衡决定，平行 biomass 使
    "宿主+异源"总蛋白碳通量≈恒定，故负担主要改 μ 而非改乙酸产率（实测 Δ乙酸≈0）。

输出（无实测数据也可用的"正向预测"，待有工艺数据时作校验基准）
  * 乙酸产率 y_ac = v(EX_ac_e) / |v(EX_<sub>_e)|  (mol/mol 碳源)；
  * 不同 (sub_max, o2_max) 下的 μ 与乙酸盐流；
  * 4 种表达形式在"问题工况"(碳过剩+DO受限)下的 μ_het 与溢出放大（葡萄糖/甘油各一份）。

运行
    cd fip && python3 -m scripts.calibrate_burden_process --write-back
    cd fip && python3 -m scripts.calibrate_burden_process --carbon-source glycerol --write-back
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np
import cobra

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from fip.twins.process.mechanistic import EXPRESSION_FORMS as EF_PY  # noqa: E402

# 复用真实 FBA 求解的异源负担结构
from scripts.calibrate_burden_fba import (  # noqa: E402
    _load_base as _load_glucose, FORM_FACTOR, AVG_RESIDUE_DA,
)
# 文献基准 + 菌株改造 × 温度修正 的宿主碳源摄取能力；及生长率自身的温度依赖(μ~Q10)
from scripts.uptake_capacity import (  # noqa: E402
    effective_uptake, describe as describe_uptake,
    mu_cap, growth_temperature_factor, describe_growth,
    MU_REF_37 as MU_NOTE_REF, T_REF as T_REF, T_MIN as T_MIN, T_OPT as T_OPT, T_MAX as T_MAX,
)

RBA_MODELS = os.path.join(HERE, "rba_models")
PROTEIN_FRACTION = 0.55
CALIBRATED_JSON = os.path.join(HERE, "calibrated_burden.json")
DOC_DIR = os.path.join(ROOT, "docs")
FIG = os.path.join(DOC_DIR, "process_overflow_scan.svg")

# 碳源配置：FIP 工艺葡萄糖或甘油二选一
CARBON = {
    "glucose": {
        "model": "e_coli_core.xml",
        "exch": "EX_glc__D_e",
        "other": "EX_glyc_e",           # 二选一 -> 封死
        "host_uptake": None,            # 运行期由 uptake_capacity.effective_uptake 注入（文献基准×菌株×温度）
        "default_sub": 15.0,            # 碳过剩（补料超摄取能力）
        "problem_o2": 8.0,             # DO 受限的"问题工况"
        "carbon_per_mol": 6,
    },
    "glycerol": {
        "model": "e_coli_core_glycerol.xml",
        "exch": "EX_glyc_e",
        "other": "EX_glc__D_e",         # 二选一 -> 封死
        "host_uptake": None,            # 运行期注入
        "default_sub": 15.0,
        "problem_o2": 8.0,
        "carbon_per_mol": 3,
    },
}

# 工艺约束默认扫描范围（物理合理、非标定值；待实测数据回填）
SUB_MAXS = [3.0, 5.0, 7.0, 10.0, 12.0, 15.0, 18.0, 20.0]   # 碳源摄取上限 mmol/gDW/h
O2_MAXS = [2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 15.0, 20.0, 30.0]  # DO 受限的呼吸氧化容量
EXPRESSION_LEVEL = 0.15    # 默认异源表达水平
CASE_LEVELS = [0.05, 0.10, 0.15, 0.20, 0.30]   # case-by-case 扫描的表达水平网格
PROTEIN_MW_KDA = 30.0       # 默认异源蛋白分子量 kDa


def _load_base(carbon: str):
    """载入指定碳源的模型（葡萄糖=e_coli_core；甘油=e_coli_core_glycerol）。"""
    cfg = CARBON[carbon]
    return cobra.io.read_sbml_model(os.path.join(RBA_MODELS, cfg["model"]))


def _add_burden(model, f: float):
    """克隆原生 biomass 为平行 biomass（化学计量 ×f），二者等速率耦合。f>0 才加。"""
    import cobra
    bm = next(r for r in model.reactions if "BIOMASS" in r.id.upper())
    bm_het = cobra.Reaction("BIOMASS_het_PROT")
    bm_het.name = "heterologous protein synthesis burden"
    model.add_reactions([bm_het])
    for met, coef in bm.metabolites.items():
        bm_het.add_metabolites({met: coef * f})
    cons = model.problem.Constraint(bm.flux_expression - bm_het.flux_expression, lb=0.0, ub=0.0)
    model.add_cons_vars(cons)
    bm_het.objective_coefficient = 0.0
    return model


def _apply_process(model, carbon: str, sub_max: float, o2_max: float):
    cfg = CARBON[carbon]
    # 二选一：封死另一碳源
    if cfg["other"] in model.reactions:
        model.reactions.get_by_id(cfg["other"]).lower_bound = 0.0
    model.reactions.get_by_id(cfg["exch"]).lower_bound = -abs(sub_max)
    model.reactions.EX_o2_e.lower_bound = -abs(o2_max)
    return model


def _metrics(carbon: str, sub_max, o2_max, f=0.0, T=37.0):
    import cobra
    model = _load_base(carbon)
    exch = CARBON[carbon]["exch"]
    if f > 0:
        _add_burden(model, f)
    _apply_process(model, carbon, sub_max, o2_max)
    # 生长率温度依赖：把细胞内禀最大生长上限 μ_cap(T) 设为 FBA 目标(对数生长/biomass)的上界，
    # 使乙酸/μ 通量分配与受限后的 μ 自洽（min(FBA_μ, μ_cap(T)) 等价于 biomass 反应 ub=μ_cap）。
    cap = mu_cap(T)
    bm = next(r for r in model.reactions if "BIOMASS" in r.id.upper())
    bm.upper_bound = cap
    sol = model.optimize()
    if sol.status != "optimal":
        return {"feasible": False, "mu": None}
    sub = -sol.fluxes.get(exch, 0.0)
    ac = sol.fluxes.get("EX_ac_e", 0.0)
    fo = sol.fluxes.get("EX_for_e", 0.0)
    et = sol.fluxes.get("EX_etoh_e", 0.0)
    co2 = sol.fluxes.get("EX_co2_e", 0.0)
    y_ac = ac / sub if sub > 1e-9 else 0.0
    return {
        "feasible": True,
        "mu": round(float(sol.objective_value), 4),
        "mu_cap": round(float(cap), 4),          # 本温度下的内禀生长上限
        "substrate_consumed": round(float(sub), 4),
        "acetate": round(float(ac), 4),
        "formate": round(float(fo), 4),
        "ethanol": round(float(et), 4),
        "co2": round(float(co2), 4),
        "acetate_yield": round(float(y_ac), 4),  # mol acetate / mol 碳源
    }


def _sweep_carbon(carbon: str, T: float = 37.0):
    cfg = CARBON[carbon]
    name = "葡萄糖" if carbon == "glucose" else "甘油"
    print(f"\n=== 碳源={name} ({carbon}) @ {T:.1f}°C ===")
    print(f"  文献基准宿主摄取能力(经菌株×温度修正) host_uptake = {cfg['host_uptake']:.3f} mmol/gDW/h")
    print(f"  生长率温度上限 μ_cap({T:.0f}°C) = {mu_cap(T):.4f} /hr（FBA 目标上界）")
    # (A) DO 受限扫描（碳过剩 sub=default_sub）
    print(f"(A) {name}摄取={cfg['default_sub']} 固定，扫 DO 氧化容量 o2_max:")
    print(f"  {'o2_max':>7}{'μ':>8}{'μ_cap':>8}{'乙酸':>9}{'甲酸':>9}{'乙醇':>8}{'y_ac':>8}")
    sweep_o2 = []
    for o2 in O2_MAXS:
        m = _metrics(carbon, cfg["default_sub"], o2, T=T)
        sweep_o2.append({"o2_max": o2, **m})
        print(f"  {o2:>7.1f}{m['mu']:>8.4f}{m['mu_cap']:>8.4f}{m['acetate']:>9.2f}{m['formate']:>9.2f}"
              f"{m['ethanol']:>8.2f}{m['acetate_yield']:>8.3f}")
    # (B) 碳源摄取上限扫描（DO 受限 o2=problem_o2）
    print(f"(B) DO 氧化容量={cfg['problem_o2']} 固定，扫{name}摄取上限 sub_max:")
    print(f"  {'sub_max':>8}{'μ':>8}{'μ_cap':>8}{'乙酸':>9}{'甲酸':>9}{'y_ac':>8}")
    sweep_sub = []
    for sub in SUB_MAXS:
        m = _metrics(carbon, sub, cfg["problem_o2"], T=T)
        sweep_sub.append({"sub_max": sub, **m})
        print(f"  {sub:>8.1f}{m['mu']:>8.4f}{m['mu_cap']:>8.4f}{m['acetate']:>9.2f}{m['formate']:>9.2f}{m['acetate_yield']:>8.3f}")
    return sweep_o2, sweep_sub


def _heterologous_carbon(carbon: str, ideal: dict, T: float = 37.0,
                         level: float = EXPRESSION_LEVEL, protein_mw: float = PROTEIN_MW_KDA):
    cfg = CARBON[carbon]
    name = "葡萄糖" if carbon == "glucose" else "甘油"
    print(f"\n(C) 问题工况({name} sub={cfg['default_sub']}, o2={cfg['problem_o2']}) @ {T:.1f}°C "
          f"下，4 种形式的 μ 对比 (level={level}, MW={protein_mw}kDa):")
    print(f"  {'形式':<22}{'μ_理想FBA':>10}{'μ_wt工艺':>10}{'μ_het工艺':>10}{'乙酸_wt':>9}{'乙酸_het':>9}")
    out = {}
    for form in EF_PY:
        ff = FORM_FACTOR.get(form, 1.0)
        f = level * ff
        mw = _metrics(carbon, cfg["default_sub"], cfg["problem_o2"], f=0.0, T=T)
        mh = _metrics(carbon, cfg["default_sub"], cfg["problem_o2"], f=f, T=T)
        out[form] = {
            "f": round(f, 4),
            "protein_mw_kda": protein_mw,    # 记录用：FBA 平行-biomass 模型下 μ_het 由 f=level×form 决定，与 MW 无关
            "mu_ideal_fba": ideal.get(form),
            "mu_wt": mw["mu"], "mu_het": mh["mu"],
            "acetate_wt": mw["acetate"], "acetate_het": mh["acetate"],
            "d_ac": round(mh["acetate"] - mw["acetate"], 3),
        }
        print(f"  {form:<22}{ideal.get(form,0):>10.4f}{mw['mu']:>10.4f}{mh['mu']:>10.4f}"
              f"{mw['acetate']:>9.2f}{mh['acetate']:>9.2f}")
    return out


def _case_by_case(carbon: str, T: float = 37.0, levels=None,
                 protein_mw: float = PROTEIN_MW_KDA):
    """case-by-case μ_het 网格：对问题工况，扫 (形式 × 表达水平) 给出 μ_wt/μ_het/乙酸。

    注：FBA 平行-biomass 模型中负担 f = level×form_factor（占干重比例），故 μ_het 由
    (form, level) 决定；protein_mw 仅作记录（固定 level 下 μ_het 与 MW 无关，因负担按
    质量而非残基数计）。"""
    cfg = CARBON[carbon]
    name = "葡萄糖" if carbon == "glucose" else "甘油"
    levels = levels or CASE_LEVELS
    print(f"\n(CB) case-by-case μ_het 网格 ({name} @ {T:.1f}°C, 问题工况 sub={cfg['default_sub']}, o2={cfg['problem_o2']}):")
    print(f"  {'形式':<22}{'level':>7}{'f':>7}{'μ_wt':>9}{'μ_het':>9}{'下降':>8}{'乙酸_het':>9}")
    grid = []
    for form in EF_PY:
        ff = FORM_FACTOR.get(form, 1.0)
        row = {"form": form}
        for lv in levels:
            f = lv * ff
            mw = _metrics(carbon, cfg["default_sub"], cfg["problem_o2"], f=0.0, T=T)
            mh = _metrics(carbon, cfg["default_sub"], cfg["problem_o2"], f=f, T=T)
            mu_wt = mw["mu"]; mu_het = mh["mu"]
            drop = (mu_wt - mu_het) / mu_wt if mu_wt > 0 else 0.0
            row.setdefault("levels", {})[lv] = {
                "f": round(f, 4),
                "mu_wt": mu_wt, "mu_het": mu_het,
                "drop": round(drop, 4),
                "acetate_het": mh["acetate"], "acetate_yield_het": mh["acetate_yield"],
            }
            print(f"  {form:<22}{lv:>7.2f}{f:>7.3f}{mu_wt:>9.4f}{mu_het:>9.4f}{drop:>8.4f}{mh['acetate']:>9.2f}")
        grid.append(row)
    return grid


def main(argv=None):
    ap = argparse.ArgumentParser(description="工艺约束(碳源摄取/DO)下的乙酸溢出模块（葡萄糖/甘油二选一）")
    ap.add_argument("--carbon-source", choices=["glucose", "glycerol", "both"], default="both",
                    help="扫描的碳源：glucose / glycerol / 二者皆扫(both, 默认)")
    ap.add_argument("--temperature", type=float, default=37.0,
                    help="工艺温度 °C（默认 37，参考温度；按 Q10 修正宿主摄取能力）")
    ap.add_argument("--strain-factor", type=float, default=1.0,
                    help="菌株改造因子（野生型=1.0；>1 过表达上调，<1 敲除/弱化下调），同时作用于两碳源")
    ap.add_argument("--strain-factor-glucose", type=float, default=None,
                    help="仅作用于葡萄糖的菌株改造因子（覆盖 --strain-factor）")
    ap.add_argument("--strain-factor-glycerol", type=float, default=None,
                    help="仅作用于甘油的菌株改造因子（覆盖 --strain-factor）")
    ap.add_argument("--level", type=float, default=EXPRESSION_LEVEL,
                    help=f"异源表达水平（占细胞干重比例，默认 {EXPRESSION_LEVEL}）；用于问题工况 μ_het 对比")
    ap.add_argument("--protein-mw", type=float, default=PROTEIN_MW_KDA,
                    help=f"异源蛋白分子量 kDa（默认 {PROTEIN_MW_KDA}；记录用，FBA 模型下 μ_het 主要由 level×form 决定）")
    ap.add_argument("--write-back", action="store_true", help="写回 calibrated_burden.json")
    args = ap.parse_args(argv)

    carbons = ["glucose", "glycerol"] if args.carbon_source == "both" else [args.carbon_source]

    # 注入文献基准 × 菌株 × 温度 修正后的宿主摄取能力
    uptake_records = {}
    for carbon in carbons:
        sf = args.strain_factor
        if carbon == "glucose" and args.strain_factor_glucose is not None:
            sf = args.strain_factor_glucose
        if carbon == "glycerol" and args.strain_factor_glycerol is not None:
            sf = args.strain_factor_glycerol
        CARBON[carbon]["host_uptake"] = effective_uptake(carbon, T=args.temperature, strain_factor=sf)
        uptake_records[carbon] = describe_uptake(carbon, T=args.temperature, strain_factor=sf)

    ideal = {}
    if os.path.exists(CALIBRATED_JSON):
        jd = json.load(open(CALIBRATED_JSON, encoding="utf-8"))
        ideal = {k: v.get("mu_het") for k, v in jd.get("forms", {}).items()}

    per_carbon = {}
    for carbon in carbons:
        sweep_o2, sweep_sub = _sweep_carbon(carbon, T=args.temperature)
        het = _heterologous_carbon(carbon, ideal, T=args.temperature,
                                   level=args.level, protein_mw=args.protein_mw)
        grid = _case_by_case(carbon, T=args.temperature, protein_mw=args.protein_mw)
        per_carbon[carbon] = {
            "config": {k: CARBON[carbon][k] for k in ("exch", "host_uptake", "default_sub", "problem_o2", "carbon_per_mol")},
            "sweep_o2_limited": sweep_o2,
            "sweep_sub_limited": sweep_sub,
            "heterologous_under_process": het,
            "case_by_case": grid,     # (形式 × 表达水平) 网格 μ_het
        }

    if args.write_back:
        data = {}
        if os.path.exists(CALIBRATED_JSON):
            data = json.load(open(CALIBRATED_JSON, encoding="utf-8"))
        data["process_overflow"] = {
            "schema": "fip.stageH.process_overflow.v5",
            "note": ("工艺约束下的乙酸溢出正向模块（v5，支持葡萄糖/甘油二选一 + 文献基准摄取能力 + 生长率温度依赖 + "
                     "case-by-case 表达水平扫描）。碳源摄取上限(sub_max)与 DO 氧化容量上限(o2_max)为工艺约束；o2_max "
                     "不足以氧化 sub_max 供来的碳时 FBA 在 biomass 目标下分泌乙酸(+甲酸+乙醇)，μ 被氧化容量限制。宿主"
                     "摄取能力 host_uptake 由 scripts/uptake_capacity.py 给出：文献基准(葡萄糖 10.5 / 甘油 2.5 "
                     "mmol/gDW/h, 37°C) × 菌株改造因子 × Q10^((T-37)/10)。生长率自身温度依赖——内禀最大生长上限 μ_max(T)"
                     "作 FBA 目标上界（biomass ub=μ_cap(T)），>40°C 因热失活回落(cardinal 模型)。甘油模型为 e_coli_core"
                     "补 GLYCt+GLYK+GLYCK 的 e_coli_core_glycerol.xml（葡萄糖-only μ 仍为 1.2917）。本模块现已与动态"
                     "发酵仿真器 simulate_ecoli 打通（SimRecipe.fba_coupled=True）：把 FBA 的 μ_cap、宿主摄取、μ_het "
                     "降幅、乙酸氧化还原 yield 灌入动态仿真的内禀参数，关掉两层各算各的；动态仿真保留时序/诱导/DO 级联。"
                     "case_by_case 段给出 (形式 × 表达水平) 网格 μ_het（FBA 平行-biomass 下 μ_het 由 f=level×form 决定，"
                     "与蛋白 MW 无关）。无实测数据时为物理合理默认扫描，待工艺数据回填作校验基准。"),
            "expression_scan": {
                "level": args.level, "protein_mw_kda": args.protein_mw,
                "case_levels": CASE_LEVELS,
            },
            "uptake_capacity": {
                "temperature_C": args.temperature,
                "strain_factor": args.strain_factor,
                "strain_factor_glucose": uptake_records.get("glucose", {}).get("strain_factor"),
                "strain_factor_glycerol": uptake_records.get("glycerol", {}).get("strain_factor"),
                "sources": uptake_records,
            },
            "growth_temperature": {
                "temperature_C": args.temperature,
                "mu_ref_37": 1.30,
                "q10_mu": 2.2,
                "T_min": 4.0, "T_opt": 40.0, "T_max": 46.0,
                "growth_temperature_factor": growth_temperature_factor(args.temperature),
                "mu_cap": mu_cap(args.temperature),
                "model": ("μ_max(T) = μ_ref_37 × cardinal(T);  T<=T_opt: Q10^((T-37)/10);  T_opt<T<T_max: "
                          "线性衰减到 0 @ T_max;  T<=T_min 或 T>=T_max: μ=0（生长停止）"),
                "figure": "growth_temperature_dependence.svg",
                "source": ("Ratkowsky (1983) square-root / Arrhenius growth (Q10≈2–2.5, Ea≈50–70 kJ/mol); "
                           "Rosso et al. (1993) cardinal temperatures E. coli (Tmin≈4–8, Topt≈39–40, Tmax≈46–48 °C); "
                           "Cell Biol by Numbers BNID 100919."),
            },
            "carbon_sources": per_carbon,
            "figures": {
                "overflow_scan": "process_overflow_scan.svg",
                "growth_temperature": "growth_temperature_dependence.svg",
            },
        }
        with open(CALIBRATED_JSON, "w", encoding="utf-8") as fh:
            json.dump(data, fh, indent=2, ensure_ascii=False)
        print(f"\n已写回：{CALIBRATED_JSON} (process_overflow 段 v5, 含 {list(per_carbon)} 碳源)")

    _draw(per_carbon, T=args.temperature)
    _draw_growth_temperature(T=args.temperature)


def _draw(per_carbon: dict, T: float = 37.0):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=(13.5, 9.0))
    titles = {"glucose": "Glucose", "glycerol": "Glycerol"}

    panels = [("glucose", "o2", 0, 0), ("glucose", "sub", 0, 1),
              ("glycerol", "o2", 1, 0), ("glycerol", "sub", 1, 1)]
    for carbon, kind, r, c in panels:
        ax = axes[r][c]
        info = per_carbon.get(carbon)
        if info is None:
            ax.set_title(f"{titles.get(carbon, carbon)}: 未扫描")
            ax.axis("off"); continue
        cfg = CARBON[carbon]
        if kind == "o2":
            sweep = info["sweep_o2_limited"]
            x = [d["o2_max"] for d in sweep]
            sublbl = f"carbon uptake fixed = {cfg['default_sub']}"
            xlabel = "O2 uptake limit (DO capacity), mmol/gDW/h"
            titleA = f"(A-{carbon[:3].upper()}) {titles[carbon]}: {sublbl}\nlower DO -> acetate overflow, mu drops"
            ax.axvline(cfg["host_uptake"], ls=":", color="k", lw=1, alpha=0.5)
        else:
            sweep = info["sweep_sub_limited"]
            x = [d["sub_max"] for d in sweep]
            sublbl = f"DO capacity fixed = {cfg['problem_o2']}"
            xlabel = "carbon uptake limit (feed capacity), mmol/gDW/h"
            titleA = f"(B-{carbon[:3].upper()}) {titles[carbon]}: {sublbl}\nabove uptake capacity mu saturates, carbon overflows"
            ax.axvline(cfg["host_uptake"], ls=":", color="k", lw=1, alpha=0.6)
            ax.text(cfg["host_uptake"] + 0.3, max(d["acetate"] for d in sweep) * 0.6,
                    f"host uptake\ncapacity ~{cfg['host_uptake']:g}", fontsize=7)
        mu = [d["mu"] for d in sweep]
        ac = [d["acetate"] for d in sweep]
        ax.plot(x, mu, "-o", ms=3, color="#1f77b4", label="mu (h^-1)")
        ax.plot(x, ac, "-s", ms=3, color="#d62728", label="acetate flux")
        ax.set_xlabel(xlabel)
        ax.set_ylabel("flux / mu")
        ax.set_title(titleA, fontsize=10)
        ax.legend(fontsize=7); ax.grid(alpha=0.3)

    fig.suptitle(f"Acetate overflow under process constraints (glucose vs glycerol, e_coli_core FBA) @ {T:.1f}°C\n"
                 "carbon feed > host uptake OR DO-limited -> overflow; heterologous burden lowers mu_het; "
                 "mu capped by growth-temperature limit mu_cap(T)",
                 fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    os.makedirs(DOC_DIR, exist_ok=True)
    fig.savefig(FIG, dpi=130)
    print(f"\n已保存图：{FIG}")


def _draw_growth_temperature(T: float = 37.0):
    """专门可视化"生长率温度依赖"扩展：μ_max(T) 上限曲线 + 双碳源摄取(T) 对比。"""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    temps = np.linspace(4.0, 46.0, 100)
    mu_curve = [mu_cap(t) for t in temps]
    f_curve = [growth_temperature_factor(t) for t in temps]
    upt_glc = [effective_uptake("glucose", T=t) for t in temps]
    upt_gly = [effective_uptake("glycerol", T=t) for t in temps]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13.5, 5.2))

    # 左：μ_max(T) 内禀生长上限（cardinal + Q10）
    ax1.plot(temps, mu_curve, "-", color="#1f77b4", lw=2, label="mu_cap(T) = μ_max(T) /hr")
    ax1.axvline(T_OPT, ls="--", color="grey", lw=1, alpha=0.7)
    ax1.text(T_OPT + 0.3, max(mu_curve) * 0.5, f"T_opt={T_OPT:.0f}°C", fontsize=8, color="grey")
    ax1.axvspan(0, T_MIN, color="red", alpha=0.06)
    ax1.axvspan(T_MAX, 50, color="red", alpha=0.06)
    ax1.axvline(T, color="#d62728", lw=2, alpha=0.8)
    ax1.text(T + 0.3, mu_cap(T) + 0.05, f"operating {T:.0f}°C\nμ_cap={mu_cap(T):.3f}",
             fontsize=8, color="#d62728")
    ax1.set_xlabel("Temperature (°C)")
    ax1.set_ylabel("μ_max(T)  (h⁻¹)")
    ax1.set_title(f"(I) Growth-rate temperature limit  μ_max(T)  [{T_REF:.0f}°C ref={MU_NOTE_REF:.2f}/h]\n"
                  "Q10 below Topt, linear decay to Tmax (cardinal, not pure Q10 extrapolation)",
                  fontsize=10)
    ax1.legend(fontsize=8); ax1.grid(alpha=0.3)
    ax1.set_xlim(4, 46)

    # 右：双碳源摄取速率(T) vs 生长上限(T) 对比
    ax2.plot(temps, upt_glc, "-", color="#2ca02c", lw=2, label="glucose uptake v_eff(T)")
    ax2.plot(temps, upt_gly, "-", color="#9467bd", lw=2, label="glycerol uptake v_eff(T)")
    ax2b = ax2.twinx()
    ax2b.plot(temps, mu_curve, "--", color="#1f77b4", lw=1.6, label="μ_max(T) (growth cap)")
    ax2b.set_ylabel("μ_max(T) (h⁻¹)", color="#1f77b4")
    ax2b.tick_params(axis="y", labelcolor="#1f77b4")
    ax2.axvline(T, color="#d62728", lw=2, alpha=0.8)
    ax2.text(T + 0.3, max(upt_glc) * 0.6, f"operating {T:.0f}°C", fontsize=8, color="#d62728")
    ax2.set_xlabel("Temperature (°C)")
    ax2.set_ylabel("host carbon uptake (mmol/gDW/h)")
    ax2.set_title("(II) Uptake(T) vs Growth-cap(T): at >40°C uptake keeps rising\nbut growth μ_max falls (heat inactivation) — decoupled temperature effects",
                  fontsize=10)
    ax2.legend(loc="upper left", fontsize=8); ax2b.legend(loc="upper right", fontsize=8)
    ax2.grid(alpha=0.3); ax2.set_xlim(4, 46)

    fig.suptitle("Growth-rate temperature dependence (μ ~ Q10 + cardinal) — FIP process module extension\n"
                 "effective μ = min(FBA_μ, μ_cap(T)); biomass reaction upper_bound = μ_cap(T)",
                 fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    out = os.path.join(DOC_DIR, "growth_temperature_dependence.svg")
    os.makedirs(DOC_DIR, exist_ok=True)
    fig.savefig(out, dpi=130)
    print(f"已保存图：{out}")


if __name__ == "__main__":
    main()
