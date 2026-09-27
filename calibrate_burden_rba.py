"""
calibrate_burden_rba.py — 用 RBA / ME 定量"表达负担"反标定 FIP 阶段 H 负担系数

目标
----
让 FIP 前向仿真里"换表达形式导致的生长下降"与基因组尺度资源分配模型
（RBA / ME-model / 酶约束 GEM）预测的"异源蛋白合成成本导致的 μ 下降"对齐。

两条计算路径
------------
  * 代理路径（默认，无需 RBApy）：rba_growth_drop_proxy —— Goelzer & Fromion 2011
    思路的极简解析近似，用于无依赖时跑通骨架。
  * 真实路径（--real-rba）：rba_growth_drop_real —— 用 RBApy 载入宿主 RBA 模型，
    计算空载 μ_wt 与带异源表达的 μ_het，drop = (μ_wt - μ_het)/μ_wt。
    异源表达反应由 _add_heterologous_expression 真实定义（加异源蛋白 macromolecule +
    翻译/折叠/分泌过程输入 + 生产 target，使蛋白质量池与核糖体容量被占用 -> μ 下降）。

标定
----
  对每个 EXPRESSION_FORM，网格搜索 burden_sol/burden_ib 使 FIP 终点下降 ≈ RBA 下降；
  tox_k 取形式预设作起点（可扩展为 3D 拟合）。
  定量结果可写回：--write-back 生成 calibrated_burden.json（单一事实源）+ 打印
  可直接粘贴回 mechanistic.EXPRESSION_FORMS 的系数块（见 §15 DoD 第 3 条）。

运行
----
    cd fip
    python -m scripts.calibrate_burden_rba                     # 默认 E.coli + 代理 RBA
    python -m scripts.calibrate_burden_rba --host Pichia
    python -m scripts.calibrate_burden_rba --real-rba \
        --rba-model scripts/rba_models/ec_core_rba \           # generate-rba-model 的输出目录
        --protein-mw 30.0 --level 0.15 --write-back           # 真实 RBA 定量标定

依赖
----
    标准库 + numpy（FIP 已装）。真实 RBA 接入需可选依赖 RBApy（未装则 --real-rba 回退代理并告警）。
"""

from __future__ import annotations

import csv
import json
import os
import sys

# ---- 把仓库根（含 fip 包）加入 path，保证任何调用方式都能 import fip.* ----
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)  # .../fip
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

import numpy as np  # noqa: E402

from fip.twins.process.mechanistic import (  # noqa: E402
    SimRecipe,
    simulate_ecoli,
    EXPRESSION_FORMS as EF_PY,
)
from fip.twins.expression.model import expression_form_to_recipe  # noqa: E402

# Pichia 侧同款接口（单一事实源 resolve_expression_burden 已 import 自 mechanistic）
try:  # pragma: no cover - 防御性导入
    from fip.twins.process.pichia_sim import PichiaRecipe, simulate_pichia  # type: ignore
    _PICHIA_OK = True
except Exception:  # pragma: no cover
    _PICHIA_OK = False


# ===========================================================================
# 0. 常量与可选依赖守卫
# ===========================================================================
# RBA 异源蛋白的氨基酸组成（按平均 E. coli 组成；单字母 AA 即模型中的 component id）。
RBA_AA_FREQ = {
    "A": 0.094, "R": 0.055, "N": 0.042, "D": 0.058, "C": 0.009,
    "Q": 0.040, "E": 0.062, "G": 0.084, "H": 0.023, "I": 0.056,
    "L": 0.094, "K": 0.059, "M": 0.024, "F": 0.040, "P": 0.047,
    "S": 0.057, "T": 0.055, "W": 0.014, "Y": 0.032, "V": 0.072,
}
_AVG_RESIDUE_DA = 110.0  # 平均氨基酸残基分子量
HET_PROTEIN_ID = "M_het_PROTEIN"
HET_TARGET_GROUP = "heterologous_expression"
# 可溶/毒性蛋白因错误折叠形成"无效循环"（合成->降解->再合成）额外占用翻译容量；
# 用等效合成率提升表示毒性维持成本（与 FIP 的 tox_k 物理对应，但作为建模选择独立于 FIP 预设）。
TOX_FUTILE_CYCLE_GAIN = 0.10

# 标定产物落盘（单一事实源；被 .gitignore 收录？不——作为守卫产物提交，见 tests）。
CALIBRATED_JSON = os.path.join(HERE, "calibrated_burden.json")


def _real_rba_available() -> bool:
    """是否已安装 RBApy（模块名小写 rba）或 COBRAme（ME-model）。"""
    import importlib.util as u  # noqa: WPS433 (局部导入，便于可选依赖守卫)
    return u.find_spec("rba") is not None or u.find_spec("cobrame") is not None


class RBASolveError(RuntimeError):
    """RBA 模型已载入但求解失败（缺 LP 求解器 / 模型结构不一致 / 合成蛋白组缺 Uniprot）。"""


# ===========================================================================
# 1. 量出 FIP 预测的「生长下降」
# ===========================================================================
def _sim(host: str, recipe):
    """按宿主分发到对应仿真器。"""
    if host == "Pichia":
        if not _PICHIA_OK:
            raise RuntimeError("pichia_sim 不可用")
        return simulate_pichia(recipe)
    return simulate_ecoli(recipe)


def _recipe_cls(host: str):
    return PichiaRecipe if (host == "Pichia" and _PICHIA_OK) else SimRecipe


def fip_growth_drop(form: str | None, burden: dict, host: str = "E.coli") -> float:
    """给定表达形式 + 负担系数，返回终点生物量下降分数 (X_base - X_burdened)/X_base。

    ⚠ 必须 expression_form=None：否则 resolve_expression_burden 会用 EXPRESSION_FORMS
    预设覆盖 burden 字段，网格搜索将完全失效（见仓库 issue 复盘）。
    """
    Rc = _recipe_cls(host)
    base = Rc(expression_form=None, ib_frac=0.3, burden_sol=0.0, burden_ib=0.0, tox_k=0.0)
    bur = Rc(expression_form=None, **{k: burden[k] for k in ("ib_frac", "burden_sol", "burden_ib", "tox_k")})
    x_base = float(_sim(host, base).X[-1])
    x_bur = float(_sim(host, bur).X[-1])
    return max(0.0, (x_base - x_bur) / x_base)


# ===========================================================================
# 2. RBA 侧：定量「表达负担」（代理 + 真实接入点）
# ===========================================================================
# 形式因子：可溶/毒性蛋白额外占伴侣与质量控制资源 -> 负担更高；包涵体因隔离 -> 更低。
# 注：真实路径的形式差异来自 machinery 路由（见 _add_heterologous_expression），
# 此处 FORM_FACTOR 仅用于代理路径的近似。
FORM_FACTOR = {
    "intracellular_soluble": 1.20,  # 可溶 + 可能 toxic：QC/伴侣成本最高
    "inclusion_body": 0.80,         # 物理隔离：对生长直接负担最低
    "periplasmic": 1.00,            # 分泌通路竞争
    "secreted": 0.90,               # 分泌：折叠/转运成本中等
}


def rba_growth_drop_proxy(expression_level: float = 0.15, form_factor: float = 1.0) -> float:
    """资源分配解析代理（Goelzer & Fromion 2011 思路的极简版）。

    异源蛋白合成占用蛋白质量池，留给代谢酶的比例下降 -> μ 按比例下降：
        drop ≈ (form_factor * level) / (capacity + form_factor * level)
    level = 异源蛋白占细胞干重比例（典型强表达 0.10–0.30，默认 0.15）。
    这是代理，用于无 RBApy 时跑通骨架；真实值来自 rba_growth_drop_real。
    """
    capacity = 1.0
    eff = form_factor * expression_level
    return max(0.0, eff / (capacity + eff))


def _load_rba_model(model_path: str):
    """载入 RBApy 模型。

    model_path 应为 generate-rba-model 的输出【目录】（含 model_file_index.in）。
    若传入单个 .xml 或不存在的路径，给出清晰错误，避免静默失败。
    """
    import importlib.util as u  # noqa: WPS433
    if u.find_spec("rba") is None and u.find_spec("cobrame") is None:
        raise NotImplementedError(
            "未安装 RBApy/cobrame：先 `pip install RBApy`（见 fip/requirements_rba.txt）。"
            "或去掉 --real-rba 使用 rba_growth_drop_proxy 代理。"
        )
    if not model_path or not os.path.exists(model_path):
        raise FileNotFoundError(
            f"需提供宿主 RBA/ME 模型【目录】（如 scripts/rba_models/ec_core_rba，"
            f"含 model_file_index.in）：当前 model_path={model_path!r}"
        )
    # 允许直接传 model_file_index.in 文件 -> 用其父目录
    if os.path.isfile(model_path):
        if model_path.endswith("model_file_index.in"):
            model_path = os.path.dirname(model_path)
        else:
            raise FileNotFoundError(
                "RBApy 模型是【目录】而非单个 .xml。请传入 generate-rba-model 的输出目录"
                f"（当前传入文件：{model_path!r}）。"
            )
    if not os.path.exists(os.path.join(model_path, "model_file_index.in")):
        raise FileNotFoundError(
            f"目录 {model_path!r} 内找不到 model_file_index.in，不是合法 RBApy 模型目录。"
        )
    if u.find_spec("rba") is not None:
        from rba import RbaModel  # 模块名小写 rba（RBApy 3.x）
        return RbaModel.from_xml(model_path)
    # 备选：COBRAme（ME-model）——接口不同，仅占位守卫
    import cobrame  # type: ignore  # pragma: no cover
    return cobrame.load_mat_model(model_path)  # pragma: no cover


def _solve_mu(model) -> float:
    """求解 RBA 模型并返回最大比生长速率 μ（Results.mu_opt）。"""
    try:
        res = model.solve()
    except Exception as e:  # LP 求解器缺失 / 结构不一致 / 合成蛋白组缺 Uniprot
        raise RBASolveError(
            f"RBA 求解失败（{type(e).__name__}: {str(e)[:160]}）。常见原因：未装 LP 求解器"
            "(glpk/cplex)、或模型用合成蛋白组生成（sandbox 无 Uniprot 访问）导致结构不一致。"
            "请在有 Uniprot 网络的机器用 `generate-rba-model` 生成真实模型后重试。"
        ) from e
    # RBApy 的 solve() 返回 rba.utils.results.Results，最优 μ 存在 .mu_opt
    mu = getattr(res, "mu_opt", None)
    if mu is None:
        mu = getattr(res, "objective_value", None)
    if mu is None:
        raise RBASolveError("RBA 求解返回无 μ（模型不可行/无界）。")
    return float(mu)


def rba_growth_drop_real(model_path, form, protein_mw_kda, expression_level):
    """【真实接入点 · 可用】用 RBApy 计算带异源表达的 μ 下降。

    步骤：
      1. 载入宿主 RBA 模型（目录，含 model_file_index.in）。
      2. μ_wt = model.solve().objective_value（空载最大比生长速率）。
      3. _add_heterologous_expression：加异源蛋白合成反应 + 容量约束（按 form 路由 machinery）。
      4. μ_het = model.solve().objective_value；drop = (μ_wt - μ_het) / μ_wt。

    模型不可解时抛 RBASolveError（由 --real-rba 调用方回退代理并告警）。
    """
    model = _load_rba_model(model_path)
    mu_wt = _solve_mu(model)
    # 注意：base 模型 solve() 会改变其内部状态，必须用【全新加载】的副本构建异源模型，
    # 否则写盘重载后会继承被污染的状态而不可行。
    model_het = _add_heterologous_expression(_load_rba_model(model_path), protein_mw_kda, expression_level, form)
    mu_het = _solve_mu(model_het)
    return max(0.0, (mu_wt - mu_het) / mu_wt)


def _add_heterologous_expression(model, protein_mw_kda, expression_level, form):
    """【真实反应定义 · RBApy 3.x】在 RBA 模型中加入异源蛋白合成反应与容量约束。

    物理机制（表达负担 = 占用宿主资源 -> μ 下降）：
      * 加异源蛋白 macromolecule M_het_PROTEIN（氨基酸组成按 MW + 平均 E.coli 组成）。
      * 把它加入「翻译」过程的生产输入（可被核糖体合成）+「蛋白降解」过程的输入
        （稳态周转，使其持续占用翻译容量）。
      * 加 production target（值 = expression_level，g/gDW/h），强制细胞以该速率合成 ->
        占用蛋白质量池 + 核糖体容量 -> 留给生长相关酶的比例下降 -> μ 下降。
      * form 修正（machinery 路由，决定 RBA 下降的大小排序）：
          - intracellular_soluble：翻译 + 折叠（QC 成本）+ 毒性无效循环（合成率 +10%）。
          - inclusion_body      ：仅翻译（物理隔离，不占折叠容量 -> 负担最低）。
          - periplasmic/secreted：翻译 + 折叠 + 分泌（转运/折叠成本中等）。
    该定义对任意合法 RBApy 模型（宿主）结构一致；具体 μ 下降由求解给出。
    """
    from rba.xml import Macromolecule, TargetSpecies, TargetGroup  # 小写 rba

    # --- 1. 异源蛋白 macromolecule（按 MW 摊氨基酸残基数）---
    n_res = max(1, int(round(protein_mw_kda * 1000.0 / _AVG_RESIDUE_DA)))
    composition = {aa: max(1, int(round(freq * n_res))) for aa, freq in RBA_AA_FREQ.items()}
    het = Macromolecule(HET_PROTEIN_ID, "Cytoplasm", composition=composition)
    model.proteins.macromolecules.append(het)

    # --- 2. 翻译 + 降解过程输入（使其可被合成/周转）---
    tr = model.processes.processes.get_by_id("P_Translation")
    _append_processing_input(tr.processings.productions[0], HET_PROTEIN_ID)
    dg = model.processes.processes.get_by_id("P_PROTEINdeg_noMM")
    _append_processing_input(dg.processings.degradations[0], HET_PROTEIN_ID)

    # --- 3. form 路由：折叠 / 分泌 machinery ---
    if form in ("intracellular_soluble", "periplasmic", "secreted"):
        fo = model.processes.processes.get_by_id("P_Folding")
        _append_processing_input(fo.processings.productions[0], HET_PROTEIN_ID)
    if form in ("periplasmic", "secreted"):
        se = model.processes.processes.get_by_id("P_Secretion_noMM")
        _append_processing_input(se.processings.productions[0], HET_PROTEIN_ID)

    # --- 4. 生产 target（强制合成率 = 等效表达水平）---
    # 可溶/毒性：错误折叠无效循环 -> 等效合成率提升，吸收额外维持/QC 容量。
    eff_level = expression_level
    if form == "intracellular_soluble":
        eff_level = expression_level * (1.0 + TOX_FUTILE_CYCLE_GAIN)
    # RBApy 的 production target value 必须是【参数 ID】（查 parameters.xml 表），
    # 不能填裸数字 -> 加一个 constant 函数作等效常量，再在 target 引用其 id。
    fn_id = "HET_PROD_{}".format(int(round(eff_level * 1000)))
    from rba.xml import Function  # 小写 rba
    model.parameters.functions.append(
        Function(fn_id, "constant", {"CONSTANT": float(eff_level)}, variable="growth_rate")
    )
    tg = TargetGroup(HET_TARGET_GROUP)
    pf = TargetSpecies(HET_PROTEIN_ID)
    pf.value = fn_id
    tg.production_fluxes.append(pf)
    model.targets.target_groups.append(tg)

    # --- 5. 写盘重载：重建参数表（含新常数函数）后再求解，保证 value 可查表 ---
    import os as _os
    import tempfile as _tf
    _tmp = _tf.mkdtemp(prefix="rba_het_")
    model.write(_tmp)
    from rba import RbaModel  # 小写 rba
    return RbaModel.from_xml(_tmp)


def _append_processing_input(processing, species_id):
    """向某 Processing 的 inputs 列表追加一个 species（stoich=1.0），复用其元素类型。"""
    ref = processing.inputs[0] if len(processing.inputs) else None
    if ref is None:
        # 兜底：直接构造 SpeciesReference（rba.xml.common）
        from rba.xml.common import SpeciesReference  # noqa: WPS433
        processing.inputs.append(SpeciesReference(species_id, 1.0))
    else:
        processing.inputs.append(type(ref)(species_id, 1.0))


# ===========================================================================
# 3. 标定：让 FIP 下降 ≈ RBA 下降
# ===========================================================================
def calibrate_form(
    form: str,
    rba_drop: float,
    host: str = "E.coli",
    grid_sol=(0.0, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0),
    ib_scale=1.0,   # burden_ib = burden_sol * ib_scale（包涵体负担走 burden_ib）
    tox_from_preset=True,
):
    """对单一表达形式，网格搜索 burden_sol/burden_ib 使 FIP 终点下降尽量逼近 rba_drop。

    仅标定主负担（burden_sol + burden_ib），tox_k 取形式预设作起点（可扩展为 3D）。
    """
    spec = expression_form_to_recipe(form)
    ib_frac = spec.get("ib_frac", 0.3)
    tox_k = spec.get("tox_k", 0.0) if tox_from_preset else 0.0
    best = None
    for bs in grid_sol:
        bi = bs * ib_scale
        burden = {"ib_frac": ib_frac, "burden_sol": bs, "burden_ib": bi, "tox_k": tox_k}
        drop = fip_growth_drop(form, burden, host=host)
        resid = abs(drop - rba_drop)
        if best is None or resid < best["resid"]:
            best = {"burden": burden, "fip_drop": drop, "resid": resid}
    return {
        "form": form, "rba_drop": rba_drop,
        "fitted": best["burden"], "fip_drop": best["fip_drop"], "resid": best["resid"],
    }


# ===========================================================================
# 4. 写回：calibrated_burden.json（单一事实源）+ EXPRESSION_FORMS 系数块
# ===========================================================================
def write_calibration_json(results, mode: str, path: str = CALIBRATED_JSON,
                           mu_wt: float | None = None, mu_het: dict | None = None):
    """把标定结果落盘为 calibrated_burden.json（§15 DoD 第 3 条守卫产物）。

    mode: 'real'（RBApy 求解成功）/ 'proxy_unvalidated'（sandbox 无 Uniprot/求解器，
    仅代理近似，待用户机器用 --real-rba 复算）。
    """
    payload = {
        "schema": "fip.stageH.burden_calibration.v1",
        "mode": mode,
        "note": ("real=由 RBApy 求解得到；proxy_unvalidated=沙箱无 Uniprot/求解器，"
                 "由代理近似得到，需在有网络的机器用 --real-rba 复算后写回。"),
        "mu_wt": mu_wt,
        "forms": {},
    }
    for r in results:
        b = r["fitted"]
        payload["forms"][r["form"]] = {
            "ib_frac": b["ib_frac"], "burden_sol": b["burden_sol"],
            "burden_ib": b["burden_ib"], "tox_k": b["tox_k"],
            "rba_drop": r["rba_drop"], "fip_drop": r["fip_drop"], "resid": r["resid"],
            "mu_het": (mu_het or {}).get(r["form"]),
        }
    with open(path, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2, ensure_ascii=False)
    return payload


def read_calibration_json(path: str = CALIBRATED_JSON) -> dict:
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def merge_calibrated_into_expression_forms(path: str = CALIBRATED_JSON) -> dict:
    """把 calibrated_burden.json 的系数合并回 EXPRESSION_FORMS（运行时视图）。

    返回 {form: {ib_frac, burden_sol, burden_ib, tox_k}}，未标定的形式保留原预设。
    """
    data = read_calibration_json(path)
    merged = {k: dict(v) for k, v in EF_PY.items()}
    for form, coeff in data.get("forms", {}).items():
        if form in merged:
            merged[form] = {
                "ib_frac": coeff["ib_frac"], "burden_sol": coeff["burden_sol"],
                "burden_ib": coeff["burden_ib"], "tox_k": coeff["tox_k"],
            }
    return merged


def _expression_forms_override_block(results) -> str:
    """生成可直接粘贴回 mechanistic.EXPRESSION_FORMS 的 Python 字典块。"""
    lines = ["EXPRESSION_FORMS = {"]
    for r in results:
        b = r["fitted"]
        lines.append(
            f'    "{r["form"]}": {{"ib_frac": {b["ib_frac"]:.2f}, '
            f'"burden_sol": {b["burden_sol"]:.2f}, "burden_ib": {b["burden_ib"]:.2f}, '
            f'"tox_k": {b["tox_k"]:.4f}}},'
        )
    lines.append("}")
    return "\n".join(lines)


# ===========================================================================
# 5. 主流程
# ===========================================================================
def main(argv=None):
    import argparse
    ap = argparse.ArgumentParser(description="用 RBA/ME 反标定 FIP 阶段 H 负担系数")
    ap.add_argument("--host", default="E.coli", choices=["E.coli", "Pichia"])
    ap.add_argument("--real-rba", action="store_true",
                    help="尝试真实 RBApy（需 RBApy + 模型目录；不可解则回退代理并告警）")
    ap.add_argument("--rba-model", default=None,
                    help="宿主 RBA 模型【目录】（generate-rba-model 输出，含 model_file_index.in）")
    ap.add_argument("--protein-mw", type=float, default=30.0, help="异源蛋白分子量 kDa")
    ap.add_argument("--level", type=float, default=0.15, help="异源蛋白占细胞干重比例（强表达 0.1–0.3）")
    ap.add_argument("--write-back", action="store_true",
                    help="把标定结果写回 calibrated_burden.json（§15 DoD 第 3 条守卫产物）")
    args = ap.parse_args(argv)

    mode = "real" if args.real_rba else "proxy_unvalidated"
    print("=== FIP 阶段 H 负担系数 · RBA/ME 反标定 ===")
    print(f"宿主: {args.host}  |  RBA 源: "
          f"{'真实(RBApy, 模型=' + (args.rba_model or '?') + ')' if args.real_rba else '代理(资源分配解析近似, level=%.2f)' % args.level}\n")

    header = f"{'形式':<22}{'RBA下降':>10}{'FIP下降':>10}{'残差':>10}  建议 burden_sol / burden_ib / tox_k"
    print(header)
    print("-" * len(header))

    results = []
    mu_het = {}
    rba_unavailable = False
    for form in EF_PY:
        ff = FORM_FACTOR.get(form, 1.0)
        if args.real_rba:
            try:
                rba_drop = rba_growth_drop_real(args.rba_model, form,
                                                protein_mw_kda=args.protein_mw,
                                                expression_level=args.level)
            except (RBASolveError, NotImplementedError, FileNotFoundError) as e:
                print(f"  [warn] {form}: 真实 RBA 不可用（{type(e).__name__}: {str(e)[:90]}）-> 回退代理")
                rba_unavailable = True
                rba_drop = rba_growth_drop_proxy(args.level, ff)
        else:
            rba_drop = rba_growth_drop_proxy(args.level, ff)
        res = calibrate_form(form, rba_drop, host=args.host)
        results.append(res)
        b = res["fitted"]
        print(f"{form:<22}{rba_drop:>10.3f}{res['fip_drop']:>10.3f}{res['resid']:>10.3f}  "
              f"{b['burden_sol']:.2f} / {b['burden_ib']:.2f} / {b['tox_k']:.4f}")

    # CSV（与脚本同目录）
    out = os.path.join(HERE, "burden_calibration.csv")
    with open(out, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["form", "host", "rba_drop", "fip_drop", "resid",
                    "burden_sol", "burden_ib", "tox_k", "ib_frac"])
        for r in results:
            b = r["fitted"]
            w.writerow([r["form"], args.host, f"{r['rba_drop']:.4f}", f"{r['fip_drop']:.4f}",
                        f"{r['resid']:.4f}", b["burden_sol"], b["burden_ib"], b["tox_k"], b["ib_frac"]])
    print(f"\n已写出标定对照表：{out}")

    if args.write_back:
        eff_mode = "proxy_unvalidated" if rba_unavailable else mode
        write_calibration_json(results, eff_mode, mu_het=mu_het)
        print(f"已写回守卫产物：{CALIBRATED_JSON}  (mode={eff_mode})")
        if rba_unavailable:
            print("⚠ 本次为代理近似（沙箱无 Uniprot/求解器）。在有网络的机器用 "
                  "`--real-rba --rba-model <dir> --write-back` 复算以得到真实定量系数。")
        print("\n可粘贴回 mechanistic.EXPRESSION_FORMS 的系数块：\n" + _expression_forms_override_block(results))
    else:
        print("提示：加 --write-back 可生成 calibrated_burden.json（§15 DoD 第 3 条守卫产物）。")


if __name__ == "__main__":
    main()
