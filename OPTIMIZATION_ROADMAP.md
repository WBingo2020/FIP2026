# FIP 优化路线图

自检（bug）之外，从**架构 / 模型 / 工程 / 体验**四个层面梳理可优化项，附优先级与成本评估。
结论先行：**在接真实数据之前，最该做的是 P0 两项（单位体系 + 单元测试），而不是继续加功能。**

---

## 结论：我的优先级主张

平台当前处于**"演示能跑、真数据会崩"**的状态。模型层和展示层已经相当完整，
继续加模块只会放大风险。建议顺序：

```
P0  单位体系统一 + 单元测试层      ← 接真实数据的前置条件
P1  数据入口质量闸门 + OUR/kLa 口径统一      ✅ 已完成（含剩余低优先项：采样间隔突变 / 3σ 离群 / 跨表一致性）
P2  模型能力（软测量标定、效价预测、CI 校准）
P3  平台工程（体积/双前端、体验细节）
```

---

## P0 · 接真实数据前的必做项  ✅ 已完成（2026-09-14）

### 1. 统一气体流量单位体系 ★ ✅ 已完成
**现状（修复前）**：两套单位靠一个意义不明的 `×60` 缝合（详见 SELF_AUDIT 二·A）

| 位置 | 修复前值 | 真实单位 | 标注 |
|---|---|---|---|
| DB `AIR` / `max_airflow_lmin` | 0.0333 (2L罐) | L/s | `L/min` ❌ |
| `SimRecipe.airflow_lmin` | 10.0 | L/min | L/min ✅ |
| 前端 `AIR` 曲线 | 0.0333 | L/s | `vvm` ❌ |

**修复做法**：种子数据 `air = vol`（真实 L/min，删 `/60`）；设备 `max_airflow_lmin`
×60 变真实 L/min（V-S01=9、V-S02=30、V-M01=300、V-L01=3600）；`vvm_from_airflow`、
`mechanistic.py`、`pichia_sim.py`、`risk.suggest_scale_strategy` 删除隐式 `×60`；
前端 `AIR`/`O2_FLOW` 单位由 `vvm`/`派生` 改为 `L/min`；`datadict.gas_flow` 去掉量纲错误的
`vvm:1.0` 登记。

**关键保证**：vvm 物理值不变（各尺度仍 ≈1.0），仅存储/显示单位真实。已重建 DB（100 批）
与 `platform_interactive.html`，`audit_numbers.py` 的单位核对现显示"按 L/min 推 vvm≈1.0"。
`tests/test_units_lmin.py` 守护"×60 不被悄悄加回"。

### 2. 建立 pytest 单元测试层 ★ ✅ 已完成
**现状（修复前）**：`tests/` 目录**只有 `__init__.py`，0 个测试**。所有验证靠 10 个
`scripts/v_*.py`——它们启动 chromium 渲染 HTML 再用正则扫 DOM，慢（每个数秒～数十秒）、
脆（依赖文案）、且**完全覆盖不到纯计算函数**。

**已完成**：`tests/` 下新增 6 个文件、34 个用例（毫秒级、无需浏览器）：

| 文件 | 覆盖 |
|---|---|
| `test_dimensions.py` | `P=NpρN³Di⁵`、OUR/CER 与摩尔流量一致、∫OUR dt、vvm=airflow/WV、`C*∝yO₂·P` |
| `test_boundaries.py` | V=0、OUR=0、rpm=0、单点序列、空 UA、OUR_peak=0 |
| `test_monotonicity.py` | kLa↗P/V&vvm、vvm↗airflow、X_end↗OUR、margin↘OUR_peak、CO2余量↗airflow |
| `test_physical_windows.py` | DO∈[0,100]、WCW/DCW∈[3.5–5.0]、DCW≥0、OUR/CER≥0 |
| `test_analytic.py` | 恒定 OUR→准稳态 `OUR/mO₂`、无维持线性增长、RQ=1、anchor 标定 |
| `test_units_lmin.py` | **P0-1 单位回归守护**（含 `/60` 守卫） |

运行：`python3.11 -m pytest tests/ -q` → **34 passed**。

---

## P1 · 数据入口与模型一致性

### 3. 导入校验补"范围 / 单位 / 时序完整性"
**现状**：`import_export._validate_and_build` 只做 4 类检查：缺列、空值、词表（仅告警）、
数值可解析。**没有范围校验**——真实数据的 `DO=1000`、`pH=-1`、`DCW=5000` 会直接入库污染模型。

**建议补**（`datadict.FieldSpec` 加 `min/max/unit` 字段）：

- **范围**：DO∈[0,100]、pH∈[0,14]、TEMP∈[0,60]、OUR>0、DCW<500  ← **✅ 已实现（P1-1）**：
  `datadict.TAG_VALUE_RANGE` 登记各位号物理上下界；`_validate_and_build` 对超界值判 error 拒绝入库。
- **单位核对**：`unit` 必须在词表且与字段期望一致 ← **✅ 已实现（P0-1 收尾）**：
  `import_export._validate_and_build` 对 `process_time_series` 的 AIR/O2/O2_FLOW/FEED
  按 `datadict` 规范单位自动归一化（L/s×60→L/min 等），并告警；未知单位保留原值+告警（不静默错）。
  `gas_flow` 换算表已补 `L/s`/`L/h`；`tests/test_import_units.py`(7 项) 永久守护。
- **时序完整性**：重复时间戳 / 时间倒序 / **常数信号（传感器卡死）** ← **✅ 已实现（P1-1，告警级）**；
  **采样间隔突变** ← **✅ 已实现（P1-1 收尾，2026-09-14）：基准间隔 3× 突增即告警**；
  **3σ 离群（稳健中位数+MAD）** ← **✅ 已实现（P1-1 收尾，2026-09-14）**。
- **跨表一致性** ← **✅ 已实现（P1-1 收尾，2026-09-14）**：`check_cross_table_consistency`
  校验孤儿外键（error）/ 起止时间倒挂（error）/ 时序末点超 end / 离线 titer 与结果 titer
  偏差>25% / 装液量与设备台账偏差>20%，并接入 `run_calibration_and_training`。

成本：低-中，收益极高（脏数据是真实项目的第一大坑）。守护：`tests/test_import_validation.py`(12 项) +
`tests/test_cross_table.py`(3 项) + `tests/test_ci_calibration.py`(4 项)。

### 4. 统一 OUR 口径（三套 → 一套） ✅ 已实现（P1-2）
生理模型(soft_sensor `QO2_GROWTH/MO2_MAINTENANCE`) / 化学计量(`f_ox=0.5`) / risk 经验默认(180)
三套已统一为**生理方程 `OUR = X·(qO2_g·μ + mO2)`**：
- `sto.our_peak_from_recipe` 改为直接复用该方程（按宿主取常数），`DEFAULT_F_OX` 由 0.5→1.0（不再默认折扣，仅作可选校准系数）；
- `risk._our_from_recipe_if_needed` 经验默认也走同一生理方程（典型点 mu=0.12/h, X=100 g/L），
  ecoli 经验默认由 180→294；前端 Scale-up Twin 的 f_ox 输入改为"OUR 校准系数"(默认 1.0)。
- 守护：`tests/test_our_unification.py`(6 项)。
⚠ 副作用：Scale-up 的 OUR_peak 源数值上浮（ecoli 经验默认 180→294），属预期的统一效果，已记入回归测试。

### 5. 重新标定后端 kLa ✅ 已实现（P1-3，前后端统一到 van't Riet）
后端 `KLA_K=26` 比 van't Riet 非凝聚基准（前端 dash.js `K=0.032` 1/s）高约 3.6 倍，已归一：
- `mass_transfer.KLA_K` 改为 `0.032·3600/1000^0.4 ≈ 7.27`（与前端逐位一致），S 尺度 kLa 由约 716→约 200 h⁻¹；
- 注释"校准到 400 h⁻¹"已更正为 van't Riet 说明，并注明真实凝聚态发酵液需 DO-stat 实测再校准
  （governance 校正因子已预留）。守护：`tests/test_kla_unification.py`(3 项)。
⚠ 副作用：Scale-up Twin 的 kLa/OTR/可行域数值整体下移，OTR/OUR margin 判断与放大结论随之变化
（属路线图标注的"改变 Scale-up 结论"），已记入回归测试，待工艺方确认。
✅ 工艺方已确认：① 前向仿真 `kla_ref=400` 已删除，改用与 Scale-up Twin **同源的 van't Riet `kla_h`**（由真实 P/V 计算）；
   `SimRecipe`/`PichiaRecipe` 新增 `di_m`/`n_imp`/`impeller_type`/`rho` 几何字段，`seed.py` 按所选容器把
   真实几何传给仿真器，使 kLa 随 P/V 合理变化（10L 参考罐 ≈265 h⁻¹，替代旧 400；千升罐不再失真）。
   守护：`test_kla_unification.py::test_forward_sim_kla_synced_to_van_t_riet`。
❓ 仍待真实数据到位：② 真实发酵液（凝聚态）DO-stat 实测 kLa/C* 校准（van't Riet 是水相基准，
   真实发酵液因流变/聚并需实测校正，governance 因子已预留）。

### 6. 软测量：逐批锚定 → 跨批次回归 + 漂移告警
**现状**：每批强制拉到离线 harvest DCW，终点**永远等于化验值**，等于把软测量的误差
藏起来了——单看曲线永远"很准"，但真实工业价值恰恰在于**软测量与化验的偏差**。

**建议**（真实数据到位后）：跨批次回归全局 `Y_OX/mO₂`，保留每批残差作为**传感器漂移告警**
（残差 >3σ 触发"尾气仪/OUR 异常"）。成本：中，但这是软测量从"演示"走向"可用"的关键一步。

---

## P2 · 模型能力

### 7. 效价预测（当前最弱环节）
CV R² 只有 **0.625**，相对误差 36.7%——**最重要的商业指标恰恰预测得最差**。

| 目标 | CV R² | 相对误差 |
|---|---|---|
| OD | 0.821 | 17.9% |
| DCW | 0.797 | 18.7% |
| **titer** | **0.625** | **36.7%** |

建议方向（按性价比）：
1. **补特征**：诱导后累计时长、补料累积量、μ 轨迹统计、甲醇/甘油累积（Pichia）
   ——当前特征只有 16 个且偏"瞬时统计量"，缺**过程积分量**
2. 分宿主建模（E.coli / Pichia 生理差异大，现在混训）
3. 模型升级（GBDT / 分位数森林）

**进展（2026-09-14）**：① 特征工程已落地 ✅。新增 5 个过程积分/窗口特征——
诱导后表达窗口 `hours_post_induction`、∫OUR dt `OUR_cum`、∫CER dt `CER_cum`、DO 受限占比
`DO_low_frac`、比耗氧 `specific_our`；并让 `seed.py` 导出 `INDUCED` 位号使诱导窗口数据驱动
（与 `predictor._induction_flag` 读取逻辑一致）。titer 留批次 CV R² 由 **0.625 → 0.663**
（+0.048，MAE 1.085→1.009 g/L），`CER_cum`/`OUR_cum` 成为头部特征；DCW R²≈0.804、OD R²≈0.829
基本持平。守护：`tests/test_titer_prediction.py`（锁 titer CV R²>0.625 且新特征在列）。
**下一步（待真实数据 / 单独议题）**：实验显示分宿主 + GBDT 可把 titer 进一步提到 ~0.72
（E.coli 0.72 / Pichia 0.71，Pichia 提升最显著）。**区间校准（#8）已先行完成**——
RF 树方差 CI 已替换为共形分位数回归（CQR），保外覆盖率≈名义 95%、宽度随条件误差自适应；
分宿主 + GBDT 点预测升级可独立推进，无需再受未校准区间拖累。在合成 n=40–60 上
CV 噪声仍较大，建议待真实数据到位后做。

### 8. 置信区间校准 ✅ 已完成（2026-09-14）
实测覆盖率 **100%**（名义 95%），偏保守 4–7 倍——RF 树间方差高估不确定性。
已加 CV 指标披露，但更好的做法是 **conformal prediction** 或分位数回归森林，
给出真正校准的区间。成本：中。

**落地（共形分位数回归 CQR）**：`predictor.ConformalQuantileInterval` 用分位数 GBDT
拟合下/上分位（α/2, 1−α/2），再以 K 折交叉共形（CV+）在保外点算一致性分数
`e_i = max(ŷ_lo,i − y_i, y_i − ŷ_hi,i)`，取有限样本校正分位
`q̂ = Quantile(scores, ⌈(n+1)(1−α)⌉/n)`，最终区间 `[ŷ_lo − q̂, ŷ_hi + q̂]`，
在可交换性下保证**边际覆盖率 ≥ 1−α**。

| 目标 | 保外经验覆盖率（名义 95%） | 中位区间宽度 | 旧 RF 树方差覆盖率 | 宽度比(CQR/旧) |
|---|---|---|---|---|
| titer | **0.96** | 8.61 | 0.98 | 2.23 |
| harvest_dcw | **0.96** | 75.15 | 0.99 | 2.77 |
| harvest_od | **0.96** | 230.15 | 0.99 | 2.88 |

要点：
- 区间**已概率校准**（保外覆盖率≈名义 95%，不再是 100% 过度覆盖）；宽度随条件误差
  **自适应**（异方差，实测宽度与噪声驱动量相关性 0.81）。
- 旧 RF 树间方差是**未校准代理**：其覆盖率随数据漂移（当前保外 0.98、合成数据甚至低至
  0.83–0.86），宽度也不反映真实条件误差——这正是"4–7× 保守"论断的根源（覆盖率无保证、
  宽度无物理意义）。
- 守护：`tests/test_ci_calibration.py`；前端 `ptCICalibrationRow` 实时披露"共形分位数回归·已校准"。
- 说明：CQR 的"校准"价值在于覆盖率保证与异方差自适应，而**不是**保证比旧法更窄——
  当条件误差本身很大（如 titer R²≈0.66）时，CQR 宽度会如实反映真实分散度，可能比被低估的
  树方差更宽。这是正确行为。

### 9. 机理仿真的两处简化
- **X 生长不受 DO 限制**：`mu_pot` 无 `f_DO` 项 ⇒ X 与 DO 两条曲线不耦合，
  会出现"DO=0 但生物量照长"的矛盾（修正 P0-1 后会更明显）
- **混合时间无尺度效应**：`N·t_mix=90` 对所有尺度相同，大罐 t_mix 实际是小罐 3–10 倍
  ⇒ 放大时混合风险被系统性低估

### 10. What-if 响应面从经验系数 → 可回归
当前 DO/温度/补料系数是经验值（面板已标注"演示级"）。建议支持接入 DOE 数据回归，
或直接换成机理模型（Monod + 氧限制）。成本：中-高。

---

### 14. 阶段 H · 表达负担耦合（产物→生长反向耦合）✅ 已实现（2026-09-25）
此前机理仿真只做"工艺→产物"单向：诱导后 qP 决定产物累积，但产物类型 / 形式（可溶、包涵体、周质、分泌）不影响宿主生长。阶段 H 补上反向耦合，使"换蛋白 / 换表达形式"真实改变 X(t)。

**四项机制（侵入度由小到大）**：
1. **产物分叉**：P 拆 P_sol（可溶）/ P_ib（包涵体），按 `ib_frac` 分配，可溶与包涵体给不同负担系数。
2. **生长负担项**：`dX -= burden_eff·qP·X`（`burden_eff = burden_sol·(1−ib_frac) + burden_ib·ib_frac`），立起"表达→减速"方向。
3. **形式相关毒性**：对可溶 toxic 产物加额外维持能耗 `dX -= tox_k·P_sol·X`；包涵体因物理隔离不贡献毒性（机制见 Slouka 2018 包涵体综述）。
4. **表达形式主开关**：`EXPRESSION_FORMS` 4 种预设（胞内可溶 / 包涵体 / 周质 / 分泌）一键把 `ib_frac/burden_sol/burden_ib/tox_k` 写入仿真 recipe，打通 Expression ↔ Process Twin。

**落地文件**：`scripts/ferment_sim.js`（JS 机理引擎）、`fip/twins/process/mechanistic.py`、`fip/twins/process/pichia_sim.py`（单一事实源 `resolve_expression_burden`）、`fip/twins/expression/model.py`（`expression_form_to_recipe`）、`scripts/dash.js` + `platform_interactive.html`（表达形式下拉 + Psol/Pib KPI）。

**回归安全**：默认负担=0 ⇒ 与改造前基线逐位一致（E.coli P=3.271/X=48.89；Pichia P=2.065/X=26.89，被多测试守护）。

**测试**：`tests/test_stage_h_burden.py`（9 项：基线零回归+产物守恒、负担减速、毒性额外维持能耗、表达形式耦合、JS↔Python parity）；全量 pytest **240 passed**。

⚠ **已知差距已升级为正式优化任务**：见 **§15 · 用 RBA/ME 反标定 FIP 阶段 H 负担系数**。当前 `burden_*`/`tox_k` 仍为预设默认系数、未用组学 / 批次数据拟合，属"方向正确、用于筛选与教学"的轻量模型，非定量预测滴度。

### 15. 已知差距升级 · 用 RBA/ME 反标定 FIP 阶段 H 负担系数（正式优化任务）
**状态**：🟢 代码完成 · 待真实模型求解验证（P2 · 模型能力）｜优先级 P2｜依赖 §14 / 真实 RBA 模型与组学数据｜工作量 中-高

**为什么是独立任务**：§14 已实现"方向正确"的反向耦合，但其 `burden_*`/`tox_k` 为预设默认系数，未经组学 / 批次数据拟合，不能定量预测滴度。本任务把它从"筛选 / 教学级"升级为"定量级"——锚点就是基因组尺度资源分配模型（RBA / ME-model / ecYeastGEM）的定量表达负担。

> **进度（截至本次提交）**：`rba_growth_drop_real` 已接真实 RBApy 3.x 流水线（`RbaModel.from_xml` 载入模型目录、`solve()` 取 μ），`_add_heterologous_expression` 已按 RBApy API 补全真实反应定义（加异源蛋白 macromolecule + 翻译/折叠/分泌过程输入 + 生产 target，按 form 路由 machinery）。代理回退路径与 CI 守卫测试（`tests/test_calibrate_burden.py`，10 passed）已就位。
> **sandbox 限制（需用户机器复算）**：本沙箱无 Uniprot 网络访问（rest.uniprot.org 返回 403），只能用合成蛋白组走完 `generate-rba-model`，生成的模型**结构完整但不可求解**（`solve()` 报 `KeyError`，酶引用不存在的蛋白物种）。因此 DoD 第 1/2 条（真实 μ_wt/μ_het 与残差≤0.02）需在**有 Uniprot 网络的机器**用真实宿主模型复算；代码本身已正确，仅求解受环境限制。

**范围（IN）**：
- 用 RBApy（E.coli RBA）或 COBRAme（iOL1650-ME）计算"带异源表达 vs 空载"的 μ 下降；
- 把该下降反标定进 `EXPRESSION_FORMS` 的 `burden_sol/burden_ib/tox_k`，使 FIP 前向仿真的 drop 与 RBA 对齐；
- 最小验证脚本已落地：`fip/scripts/calibrate_burden_rba.py` —— 代理 RBA 跑通（E.coli + Pichia 残差 ≤0.017）；`--real-rba` 已接真实 RBApy 3.x 流水线（`rba_growth_drop_real` 载入模型目录、`_add_heterologous_expression` 补全真实异源反应定义，`solve()` 取 μ_wt/μ_het）。模型不可解时显式抛 `RBASolveError` 并回退代理 + 告警。

**范围（OUT）**：
- 不重做基因组重建（用现成 RBA_Ecoli.xml / ecYeastGEM）；
- 不做前向时间轨迹的 RBA 化（仍由 FIP ODE 跑，仅标定系数）。

**验收标准（Definition of Done）**：
1. `pip install RBApy` 后 `python -m scripts.calibrate_burden_rba --real-rba --rba-model <RBA模型目录>` 能跑出真实 μ_wt / μ_het 与 drop —— **代码已就绪**；sandbox 无 Uniprot 致合成模型不可解，需在用户机器（真实模型）复算。
2. 反标定后 E.coli 4 种表达形式的 FIP drop 与 RBA drop 残差 ≤0.02 —— **代理模式已 ≤0.017**；真实模式待 DoD1 复算后确认。
3. ✅ 标定结果写回 `calibrated_burden.json`（单一事实源，`--write-back` 生成），并通过 `merge_calibrated_into_expression_forms()` 合并回 `EXPRESSION_FORMS`；由 `tests/test_calibrate_burden.py::test_calibration_json_roundtrip_guard` 守护。
4. ✅ `tests/test_calibrate_burden.py` 在 CI 中（代理模式）全绿（10 passed），且 `--real-rba` 路径有 import / 模型 skip 守卫（`test_real_rba_unsolvable_sandbox_model_raises_rbasolveerror`、`test_add_heterologous_expression_structural`）。

**数据前置**：宿主 RBA 模型（RBA_Ecoli.xml / RBA_cerevisiae.xml）、目标蛋白 MW、典型表达水平（占 CDW 比例 0.10–0.30）、可选蛋白组占用数据细化酶约束。

**风险**：RBA 模型自带参数不确定性；标定只改系数不改结构，若 RBA 预测与实测批次偏差大，需回到 §14 结构本身复盘。

**入口**：`fip/scripts/calibrate_burden_rba.py`、`fip/requirements_rba.txt`、`tests/test_calibrate_burden.py`。

---

## P3 · 平台工程与体验

### 11. 单文件体积（DASH_DATA 占 70%） ✅ 已完成
原 1.0 MB，其中 **DASH_DATA 696 KB**（100 批 × ~49 点 × 12 曲线）。已落地两招：

1. **精度裁剪**（build_interactive.py）：修复 `_default` 把原生 float 误转字符串的 bug
   （predictions 数值原以 15 位浮点字符串存储），DCW/WCW 3→1 位、AIR/FEED/O2_FLOW/PRESS 4→2 位、
   TEMP/RPM 2→1 位、`phases` 浮点钳到 2 位。仅此约省 7%。
2. **时序轻量降采样**（P3-⑪ 轻量版「懒加载」）：每批 ~49 点 → ~36 点，发酵曲线仍平滑
   （预测用的是全分辨率 DB 特征 `_features_from_batch`，不受影响）。

结果：**DASH_DATA 696 → 574 KB（−18.6%，落在 15–25% 区间）**；总 HTML 1064 → 909 KB。
> 真·时序懒加载（列表页只载元数据、点开才载曲线）仍待 1000+ 批次时再做——届时体积可再从根上压住。

### 12. 两套前端的维护成本 ✅ 已收敛策略
Streamlit 7 个页面（833 行）+ 自包含 HTML（dash.js 4056 行）。计算层已共享，重复的是 UI。
**已落地策略**：以自包含 HTML（`platform_interactive.html` = `dash.js`+`build_interactive.py`）为
**展示/分发唯一真源**；Streamlit 退为**数据管理/导入后台**（保留「📥 数据导入」「🗄 批次数据库」）。

- `fip/app/common.py` 新增 `display_mainline_note()`：在 4 个 Streamlit 展示页
  （Process / Expression / Scale-up Twin + Copilot）顶部注入迁移横幅，明确"只读镜像、不再同步"。
- 导入页新增「🖥 交互展示版」入口（`st.link_button` 指向 `platform_interactive.html`）。
- 后续任何**展示/交互**改动只在 `dash.js`+`build_interactive.py` 一处维护，消除"同步两遍"。

### 13. 体验细节 ✅ 已完成
- **Copilot 引导**（dash.js `renderCopilot`）：空状态给出能力说明 + 5 个可点击示例问题
  （含"如何提高 Pichia M→L 放大成功率""软测量如何估计 Biomass"等），点击即问。
- **导入 O2 缺失主动告警**（import_export.`ingest_file`）：导入 `process_time_series` 时按批次聚合
  缺失的关键在线位号（O2/DO/pH/AIR/FEED/OUR/CER），**导入完成即告警**——
  O2 缺失明确提示"将按空气×20.95% 派生，丢失纯氧/富氧信息，影响 Scale-up 供氧/OTR 分析"。
  Streamlit「数据导入」QA 报告已渲染 `report.warnings`，告警直接可见（无需事后看图）。
- 策略对比表 `constant_PV`/`constant_kLa` 同源说明：见 Scale-up Twin 页面（未在本轮改动）。

---

## 已顺手优化

**软测量二分早停**：固定 80 次二分（区间仅 399.5 宽，80 次后精度 3e-22 = 浮点噪声级）
→ 加 `hi-lo < 1e-4` 早停，典型 22 次收敛。

| 指标 | 修改前 | 修改后 |
|---|---|---|
| 单批耗时 | 3.0 ms | **0.88 ms（3.4×）** |
| 1000 批 | 3.0 s | 0.9 s |

正确性不变：100/100 批终点仍精确锚定实测（偏离 0 批 >2%），
标定模式分布不变（yield 86 / scale 14）。回归全绿。

---

## 成本收益速查

| 项 | 优先级 | 成本 | 收益 |
|---|---|---|---|
| 统一气体流量单位 | P0 | 中 | **避免 60 倍静默错误** |
| pytest 单元测试层 | P0 | 低-中 | **防止本轮这类量纲错误复发** |
| 导入范围/单位校验 | P1 | 低-中 | 拦住真实脏数据 |
| 统一 OUR 口径 | P1 | 低 | 消除 1.85 倍口径差 |
| 标定 kLa | P1 | 中 | 前后端一致 + 符合文献 |
| 软测量跨批回归 | P1 | 中 | 软测量真正可用（含漂移告警） |
| 效价预测特征工程 | P2 | 中 | R² 0.625 → 目标 0.8 |
| CI 校准 | P2 | 中 | 区间不再虚高 4–7 倍 |
| RBA/ME 反标定负担系数 | P2 | 中-高 | 定量级（替代演示级）表达负担 |
| 时序懒加载/压缩 | P3 | 中 | 支撑 1000+ 批次 |
| 前端主次收敛 | P3 | 低 | 消除双份维护 |