"""
阶段 H 负担系数反标定骨架的回归测试。

覆盖：
  * 代理 RBA 模式（无需 RBApy / cobra）：E.coli + Pichia 双侧标定可跑、残差小（<0.05）。
  * 零负担 -> 终点下降≈0；有负担 -> 下降严格更大（单调性）。
  * CSV 写出结构正确。
  * 真实 RBA 接入点守卫：缺依赖 / 缺模型显式报错；sandbox 合成模型不可解时抛 RBASolveError
    （便于 --real-rba 安全回退代理）。
  * _add_heterologous_expression 结构有效性（载入模型->加异源蛋白->写回重载无错）。
  * calibrated_burden.json 写回守卫：json 是单一事实源，合并回 EXPRESSION_FORMS 可往返。

真实 RBA 模式的【求解】不在本测试覆盖；由 CI 的 import / 模型 skip 守卫负责。
"""

import csv
import glob
import os
import sys

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # .../fip
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from fip.twins.process.mechanistic import EXPRESSION_FORMS  # noqa: E402
from scripts.calibrate_burden_rba import (  # noqa: E402
    FORM_FACTOR,
    RBASolveError,
    calibrate_form,
    fip_growth_drop,
    merge_calibrated_into_expression_forms,
    rba_growth_drop_proxy,
    rba_growth_drop_real,
    write_calibration_json,
)


def _generated_rba_dir():
    hits = glob.glob(os.path.join(ROOT, "scripts", "rba_models", "**", "model_file_index.in"),
                     recursive=True)
    return os.path.dirname(hits[0]) if hits else None


def _rba_installed():
    import importlib.util as u
    return u.find_spec("rba") is not None or u.find_spec("cobrame") is not None


# ---------------------------------------------------------------------------
# 代理路径
# ---------------------------------------------------------------------------
def test_zero_burden_drop_is_zero():
    drop = fip_growth_drop(
        None, {"ib_frac": 0.3, "burden_sol": 0.0, "burden_ib": 0.0, "tox_k": 0.0},
        host="E.coli",
    )
    assert drop < 1e-6


def test_burden_reduces_growth_monotonic():
    d0 = fip_growth_drop(
        None, {"ib_frac": 0.1, "burden_sol": 0.0, "burden_ib": 0.0, "tox_k": 0.0},
        host="E.coli",
    )
    d1 = fip_growth_drop(
        None, {"ib_frac": 0.1, "burden_sol": 0.8, "burden_ib": 0.4, "tox_k": 0.0},
        host="E.coli",
    )
    assert d1 > d0


@pytest.mark.parametrize("host", ["E.coli", "Pichia"])
def test_calibrate_proxy_residual_small(host):
    for form in EXPRESSION_FORMS:
        rba_drop = rba_growth_drop_proxy(0.15, FORM_FACTOR.get(form, 1.0))
        res = calibrate_form(form, rba_drop, host=host)
        assert res["resid"] < 0.05
        assert res["fitted"]["burden_sol"] >= 0.0
        assert res["fitted"]["ib_frac"] == EXPRESSION_FORMS[form]["ib_frac"]


def test_csv_structure(tmp_path):
    results = []
    for form in EXPRESSION_FORMS:
        rba_drop = rba_growth_drop_proxy(0.15, FORM_FACTOR.get(form, 1.0))
        results.append(calibrate_form(form, rba_drop, host="E.coli"))
    out = tmp_path / "burden_calibration.csv"
    with open(out, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["form", "host", "rba_drop", "fip_drop", "resid",
                    "burden_sol", "burden_ib", "tox_k", "ib_frac"])
        for r in results:
            b = r["fitted"]
            w.writerow([r["form"], "E.coli", f"{r['rba_drop']:.4f}", f"{r['fip_drop']:.4f}",
                        f"{r['resid']:.4f}", b["burden_sol"], b["burden_ib"], b["tox_k"], b["ib_frac"]])
    rows = list(csv.reader(open(out)))
    assert len(rows) == len(EXPRESSION_FORMS) + 1
    assert rows[0][0] == "form"
    assert float(rows[1][2]) >= 0.0  # rba_drop 列


# ---------------------------------------------------------------------------
# 真实 RBA 接入点守卫
# ---------------------------------------------------------------------------
def test_real_rba_missing_model_raises():
    # 模型路径缺失/为 None -> 显式 FileNotFoundError（区分"未装依赖"与"缺模型"）。
    with pytest.raises((NotImplementedError, FileNotFoundError)):
        rba_growth_drop_real(None, "intracellular_soluble", 30.0, 0.15)


def test_real_rba_single_xml_rejected():
    # 传单个 .xml（非模型目录）应给出清晰错误，避免静默失败。
    with pytest.raises((NotImplementedError, FileNotFoundError)):
        rba_growth_drop_real("/no/such/model.xml", "intracellular_soluble", 30.0, 0.15)


@pytest.mark.skipif(not _rba_installed(), reason="RBApy/cobrame 未安装（可选依赖）")
@pytest.mark.skipif(_generated_rba_dir() is None, reason="无已生成的 RBA 模型目录（需 generate-rba-model）")
def test_real_rba_unsolvable_sandbox_model_raises_rbasolveerror():
    # sandbox 合成蛋白组模型结构不一致（缺 Uniprot）-> 求解抛 RBASolveError，
    # 使 --real-rba 能安全回退代理而非"看似成功"。
    with pytest.raises(RBASolveError):
        rba_growth_drop_real(_generated_rba_dir(), "intracellular_soluble", 30.0, 0.15)


@pytest.mark.skipif(not _rba_installed(), reason="RBApy/cobrame 未安装（可选依赖）")
@pytest.mark.skipif(_generated_rba_dir() is None, reason="无已生成的 RBA 模型目录（需 generate-rba-model）")
def test_add_heterologous_expression_structural():
    # _add_heterologous_expression 在真实模型结构上可运行、写回重载无结构性错误。
    from rba import RbaModel
    mdir = _generated_rba_dir()
    m = RbaModel.from_xml(mdir)
    n0 = len(m.proteins.macromolecules)
    from scripts.calibrate_burden_rba import _add_heterologous_expression
    _add_heterologous_expression(m, 30.0, 0.15, "intracellular_soluble")
    assert len(m.proteins.macromolecules) == n0 + 1
    import tempfile
    td = tempfile.mkdtemp()
    m.write(td)
    m2 = RbaModel.from_xml(td)  # 重载无错 = 结构有效
    assert len(m2.proteins.macromolecules) == n0 + 1


# ---------------------------------------------------------------------------
# 写回守卫（§15 DoD 第 3 条）
# ---------------------------------------------------------------------------
def test_calibration_json_roundtrip_guard(tmp_path):
    results = []
    for form in EXPRESSION_FORMS:
        rba_drop = rba_growth_drop_proxy(0.15, FORM_FACTOR.get(form, 1.0))
        results.append(calibrate_form(form, rba_drop, host="E.coli"))
    json_path = tmp_path / "calibrated_burden.json"
    write_calibration_json(results, "proxy_unvalidated", path=str(json_path))
    merged = merge_calibrated_into_expression_forms(str(json_path))
    assert set(merged) == set(EXPRESSION_FORMS)
    for form, coeff in merged.items():
        assert {"ib_frac", "burden_sol", "burden_ib", "tox_k"} <= set(coeff)
        # 未标定形式保留原预设；标定形式取 json 值
        assert coeff["ib_frac"] == EXPRESSION_FORMS[form]["ib_frac"]
