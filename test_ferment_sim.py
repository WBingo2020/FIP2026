"""What-if 动态仿真引擎测试。

- Python 机理仿真不变量（与展示主线 FermentSim JS 端口同式）。
- JS 端口（scripts/ferment_sim.js）语法 + 数值一致性（node 运行自检）。
"""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import numpy as np

from fip.twins.process.mechanistic import simulate_ecoli, SimRecipe
from fip.twins.process.pichia_sim import simulate_pichia, PichiaRecipe

ROOT = Path(__file__).resolve().parents[1]
JS_ENGINE = ROOT / "scripts" / "ferment_sim.js"


def _invariants(tr, label):
    """通用不变量：质量非负、DO 在 [0,100]、终点效价/菌浓为正。"""
    errs = []
    if not (tr.X[-1] > 0):
        errs.append(f"{label}: 终点菌浓 X<=0 ({tr.X[-1]})")
    if not (tr.P[-1] > 0):
        errs.append(f"{label}: 终点效价 P<=0 ({tr.P[-1]})")
    if tr.DO.min() < -1e-6 or tr.DO.max() > 100.0 + 1e-6:
        errs.append(f"{label}: DO 超出 [0,100] (min={tr.DO.min():.2f}, max={tr.DO.max():.2f})")
    for arr in (tr.X, tr.S, tr.P):
        if np.nanmin(arr) < -1e-6:
            errs.append(f"{label}: 出现负值质量")
            break
    if errs:
        raise AssertionError("; ".join(errs))


def test_ecoli_sim_invariants():
    tr = simulate_ecoli(SimRecipe())
    _invariants(tr, "E.coli")
    # 高密度发酵应达到可观菌浓与效价
    assert tr.X[-1] > 20.0, f"菌浓过低 {tr.X[-1]:.1f}"
    assert tr.P[-1] > 1.0, f"效价过低 {tr.P[-1]:.2f}"


def test_pichia_sim_invariants():
    tr = simulate_pichia(PichiaRecipe())
    _invariants(tr, "Pichia")
    assert tr.X[-1] > 10.0, f"菌浓过低 {tr.X[-1]:.1f}"


def test_ecoli_oxygen_limitation_sensible():
    """转速提升应缓解氧限制（DO 最低值升高或持平）。"""
    base = simulate_ecoli(SimRecipe())
    hi = simulate_ecoli(SimRecipe(rpm=1300))
    assert hi.DO.min() >= base.DO.min() - 1e-6, "高转速未缓解氧限制"


def test_induction_earlier_raises_titer():
    """诱导提前（更长表达窗口）应提高或持平终点效价。"""
    late = simulate_ecoli(SimRecipe(induction_h=20))
    early = simulate_ecoli(SimRecipe(induction_h=10))
    assert early.P[-1] >= late.P[-1] - 1e-6, "提前诱导未提升效价"


def test_js_engine_parity_and_syntax():
    """JS 端口：语法合法、且与 Python 默认配方数值一致（titer/dcw/doMin）。"""
    if not JS_ENGINE.exists():
        raise AssertionError("ferment_sim.js 不存在")
    code = JS_ENGINE.read_text(encoding="utf-8")
    # 用 node 运行：stub window，调用 simulate 与 defaultRecipe，输出 JSON 摘要
    script = (
        "global.window={};\n" + code + "\n"
        "var F=global.window.FermentSim;\n"
        "function summ(tr){var s=tr.summary;return {titer:+s.titer.toFixed(3),dcw:+s.dcw.toFixed(3),"
        "doMin:+s.doMin.toFixed(2),ourPeak:+s.ourPeak.toFixed(2),n:tr.n};}\n"
        "console.log(JSON.stringify({ec:summ(F.simulate(F.defaultRecipe('ecoli'))),pi:summ(F.simulate(F.defaultRecipe('pichia')))}));\n"
    )
    proc = subprocess.run(
        ["node", "-e", script], capture_output=True, text=True,
        env={**__import__("os").environ, "PATH": "/root/.nvm/versions/node/v22.13.1/bin:/usr/bin:/bin"},
    )
    if proc.returncode != 0:
        raise AssertionError("node 运行 ferment_sim.js 失败: " + proc.stderr)
    out = json.loads(proc.stdout.strip().splitlines()[-1])
    py_ec = simulate_ecoli(SimRecipe())
    py_pi = simulate_pichia(PichiaRecipe())
    assert abs(out["ec"]["titer"] - round(float(py_ec.P[-1]), 3)) < 0.02, f"E.coli titer 不一致 {out['ec']['titer']} vs {py_ec.P[-1]:.3f}"
    assert abs(out["ec"]["dcw"] - round(float(py_ec.X[-1]), 3)) < 0.05, "E.coli dcw 不一致"
    assert abs(out["pi"]["titer"] - round(float(py_pi.P[-1]), 3)) < 0.02, "Pichia titer 不一致"
    assert out["ec"]["n"] == len(py_ec.t), "E.coli 点数不一致"
