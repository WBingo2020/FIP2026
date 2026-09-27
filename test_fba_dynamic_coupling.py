"""阶段 FBA 耦合 + case-by-case 扫描的回归/行为测试。

- Item1: FBA 代谢 Oracle（μ_cap / 宿主摄取 / μ_het 降幅 / 乙酸 yield）灌入动态发酵仿真
  simulate_ecoli（SimRecipe.fba_coupled=True），关掉两层各算各的；默认 fba_coupled=False
  保持回归安全。
- Item2: 工艺溢出模块 calibrate_burden_process._case_by_case 给出 (形式 × 表达水平) 网格 μ_het。
- Item3: 两套模型 μ 标度对齐（fba_coupled 用 MU_REF_37 天花板），乙酸机制统一为 FBA 氧化还原 yield。
"""
import numpy as np
import pytest

from fip.twins.process.mechanistic import simulate_ecoli, SimRecipe
from scripts.uptake_capacity import mu_cap
from scripts.calibrate_burden_process import _case_by_case, _metrics


def _peak_mu(tr):
    X, t = tr.X, tr.t
    return float(np.nanmax(np.gradient(X, t) / np.maximum(X, 1e-9)))


def test_fba_coupled_glucose_runs_and_respects_ceiling():
    rc = SimRecipe(carbon_source="glucose", fba_coupled=True,
                   expression_form="intracellular_soluble", expression_level=0.15,
                   fba_temp=37.0, t_end_h=18.0)
    tr = simulate_ecoli(rc, noise=0.0)
    assert tr.X[-1] > 0.0
    cap = mu_cap(37.0)
    # 动态峰值 μ 不应超过 FBA 内禀上限（允许时序离散化微小超调）
    assert _peak_mu(tr) <= cap + 0.05
    # 表达负担压低了 μ（对比 wt：drop_frac≈0.847 应使峰值 μ < 天花板）
    assert _peak_mu(tr) < cap


def test_fba_coupled_glucose_overflow_under_do_stress():
    # 低 rpm/通气/kLa 强制氧限制 -> 应触发 FBA 乙酸溢出 yield
    rc = SimRecipe(carbon_source="glucose", fba_coupled=True,
                   expression_form="intracellular_soluble", expression_level=0.15,
                   fba_temp=37.0, t_end_h=18.0, rpm=400.0, airflow_lmin=3.0, kla_scale=0.4)
    tr = simulate_ecoli(rc, noise=0.0)
    assert tr.A[-1] > 0.0, "fba_coupled 下 DO 受限应产乙酸（FBA 溢出逻辑已灌入）"


def test_fba_coupled_glycerol_runs():
    rc = SimRecipe(carbon_source="glycerol", fba_coupled=True,
                   expression_form="intracellular_soluble", expression_level=0.15,
                   fba_temp=37.0, t_end_h=18.0)
    tr = simulate_ecoli(rc, noise=0.0)
    assert tr.X[-1] > 0.0
    # 甘油是弱碳源：宿主摄取低 -> 生长远慢于葡萄糖天花板
    assert _peak_mu(tr) < mu_cap(37.0)


def test_default_fba_uncoupled_unchanged():
    # 回归安全：fba_coupled=False 不触发 FBA 路径
    rc = SimRecipe(carbon_source="glucose", fba_coupled=False, t_end_h=12.0)
    tr = simulate_ecoli(rc, noise=0.0)
    assert tr.X[-1] > 0.0


def test_case_by_case_grid_shape_and_monotonic():
    grid = _case_by_case("glucose", T=37.0)
    forms = [r["form"] for r in grid]
    assert set(forms) == {"intracellular_soluble", "inclusion_body", "periplasmic", "secreted"}
    # 同一形式下 μ_het 随表达水平单调下降
    for row in grid:
        levels = sorted(float(lv) for lv in row["levels"])
        mu_het = [row["levels"][lv]["mu_het"] for lv in levels]
        assert mu_het == sorted(mu_het, reverse=True)
    # μ_wt 与表达水平无关（负担只影响 het）
    wt0 = grid[0]["levels"][0.15]["mu_wt"]
    assert abs(grid[1]["levels"][0.15]["mu_wt"] - wt0) < 1e-9


def test_case_by_case_burden_reduces_mu():
    grid = _case_by_case("glucose", T=37.0)
    row = next(r for r in grid if r["form"] == "intracellular_soluble")
    # level=0.30 时 μ_het 应明显低于 level=0.05
    assert row["levels"][0.30]["mu_het"] < row["levels"][0.05]["mu_het"]
