"""合成数据生成器：E. coli 高密度补料分批发酵（机理仿真 + 噪声）。

模型（对应《5.2 模型结构》）：
- 生长：Monod + 乙酸抑制；生长受碳平衡约束（不可消耗不存在的底物）
- 溢流代谢：摄取超过生长需求 -> 乙酸积累 -> 抑制
- 产物：诱导后非生长偶联表达（T7）
- 气体：OUR = X·(qO2_growth·μ + mO2)；CER≈1.05·OUR（呼吸商修正）；RQ=CER/OUR
- DO：OTR = kLa·(C*−CL)，DO 低时 RPM/富氧级联；生长受 DO 限制（f_DO=CL/(ki_do+CL)）
  耦合 X 与 DO（氧限制时 μ 被压低，原模型 X 与 DO 两条曲线不耦合）
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

from fip.twins.scaleup.mass_transfer import kla_h
from fip.twins.scaleup.power import power_per_volume_w_m3
from fip.twins.process.media import params_for_medium, recipe_for_medium


# ---------------------------------------------------------------------------
# 阶段 H · 表达负担耦合：产物分叉（可溶/包涵体）+ 生长负担 + 可溶毒性
# 与 scripts/ferment_sim.js 的 EXPRESSION_FORMS / resolveBurden 逐位一致。
#   ib_frac     : 包涵体占产物比例（0..1）
#   burden_sol  : 可溶产物合成对生长的负担系数 λsol (gX/gP)
#   burden_ib   : 包涵体合成对生长的负担系数 λib (gX/gP)
#   tox_k       : 可溶（toxic）产物的额外维持/死亡能耗系数 1/(g/L·h)（IB 因隔离不贡献毒性）
# 默认（未选形式、未调滑块）负担全 0 ⇒ 与改造前行为逐位一致（回归安全）。
# ---------------------------------------------------------------------------
EXPRESSION_FORMS = {
    "intracellular_soluble": {"ib_frac": 0.10, "burden_sol": 0.80, "burden_ib": 0.20, "tox_k": 0.004},
    "inclusion_body":        {"ib_frac": 0.85, "burden_sol": 0.20, "burden_ib": 0.15, "tox_k": 0.000},
    "periplasmic":           {"ib_frac": 0.15, "burden_sol": 0.40, "burden_ib": 0.15, "tox_k": 0.001},
    "secreted":              {"ib_frac": 0.05, "burden_sol": 0.25, "burden_ib": 0.10, "tox_k": 0.000},
}
DEFAULT_IB_FRAC = 0.3
DEFAULT_BURDEN_SOL = 0.0
DEFAULT_BURDEN_IB = 0.0
DEFAULT_TOX_K = 0.0


def resolve_expression_burden(rc):
    """解析有效负担参数：给定表达形式用其预设（主开关）；否则读字段（默认 0 ⇒ 回归安全）。"""
    form = getattr(rc, "expression_form", None)
    preset = EXPRESSION_FORMS.get(form) if form else None
    if preset:
        return preset["ib_frac"], preset["burden_sol"], preset["burden_ib"], preset["tox_k"]
    return (
        rc.ib_frac if rc.ib_frac is not None else DEFAULT_IB_FRAC,
        rc.burden_sol if rc.burden_sol is not None else DEFAULT_BURDEN_SOL,
        rc.burden_ib if rc.burden_ib is not None else DEFAULT_BURDEN_IB,
        rc.tox_k if rc.tox_k is not None else DEFAULT_TOX_K,
    )


def _fba_couple(rc: "SimRecipe", p: "EColiParams") -> dict:
    """阶段 FBA 耦合：用 scripts/calibrate_burden_process 的 FBA 代谢 Oracle 校准动态仿真的
    内禀参数，关掉"两层各算各的"。惰性导入避免与 mechanistic 的循环依赖。

    返回灌入量：drop_frac(FBA μ_het/μ_wt 降幅)、y_ac_gg(FBA 乙酸产率 g/g)、mu_ceiling(FBA
    内禀生长上限)、并就地改写 p.mu_max=MU_REF_37、p.q_s_max=宿主摄取(转 gS/gX/h)。
    """
    from scripts.calibrate_burden_process import _metrics  # noqa: E402
    from scripts.uptake_capacity import (  # noqa: E402
        effective_uptake, mu_cap, MU_REF_37, growth_temperature_factor,
    )
    from scripts.calibrate_burden_fba import FORM_FACTOR  # noqa: E402
    carbon = rc.carbon_source
    T = rc.fba_temp
    sf = rc.strain_factor
    # FBA 代表性工况：宿主摄取能力为上限、氧充足（用于隔离表达负担的 μ 降幅）
    sub_max = effective_uptake(carbon, T=T, strain_factor=sf)
    o2_max = 30.0
    o2_limited = 8.0                              # DO 受限工况（取乙酸溢出 yield，用于动态仿真溢流）
    form = rc.expression_form
    ff = FORM_FACTOR.get(form, 1.0)
    f = rc.expression_level * ff
    mw = _metrics(carbon, sub_max, o2_max, f=0.0, T=T)
    mh = _metrics(carbon, sub_max, o2_max, f=f, T=T)
    mu_wt = mw.get("mu") or 0.0
    mu_het = mh.get("mu") or 0.0
    drop_frac = (mu_het / mu_wt) if mu_wt > 0 else 1.0
    # 乙酸产率取 DO 受限工况（动态仿真在 DO 受限/补料超摄取时才会溢流，参考点 o2=30 无溢出）
    mh_ac = _metrics(carbon, sub_max, o2_limited, f=f, T=T)
    y_ac = mh_ac.get("acetate_yield", 0.0)        # mol 乙酸 / mol 碳源
    mm = 180.0 if carbon == "glucose" else 92.0  # 碳源摩尔质量 g/mol
    y_ac_gg = y_ac * 60.0 / mm                    # g 乙酸 / g 碳源
    # 内禀参数对齐到 FBA 物理基准（动态仿真每步再乘 cardinal 温度因子 -> μ_cap(步温)）
    p.mu_max = MU_REF_37
    p.q_s_max = sub_max * mm / 1000.0             # mmol/gDW/h -> gS/gX/h
    return {
        "carbon": carbon, "T": T, "strain_factor": sf,
        "drop_frac": drop_frac, "y_ac_gg": y_ac_gg,
        "mu_ceiling": mu_cap(T), "q_s_max": p.q_s_max,
        "gtf": growth_temperature_factor,         # FBA cardinal 温度因子（替代内生 exp 形式）
    }


@dataclass
class EColiParams:
    # μ 标度说明（阶段 FBA 对齐）：mu_max=0.55 是"实现值"——高密度补料分批中受底物/DO 限制后
    # 的实际最大比生长（与观测批次数据拟合）。FBA 内禀上限 μ_cap(37°C)=1.30/h 是碳/氧不限时的
    # 理论天花板。二者不矛盾：0.55 < 1.30 正是"补料受限"的体现。fba_coupled=True 时 mu_max 被
    # 替换为 MU_REF_37(=μ_cap 参考)，使动态仿真与 FBA 共享同一内禀物理基准（每步再乘 cardinal
    # 温度因子 -> μ_cap(步温)）。详见 _fba_couple。
    mu_max: float = 0.55          # 1/h（实现值；fba_coupled 时改为 FBA 内禀上限 MU_REF_37）
    ks: float = 0.05              # g/L
    yxs: float = 0.45             # gX/gS
    q_s_max: float = 1.2          # gS/gX/h
    q_p_max: float = 0.012        # gP/gX/h 诱导后
    ki_acetate: float = 4.0       # g/L
    qo2_growth: float = 12.0      # mmol O2/gX（生长相关，折合每单位 μ）
    mo2: float = 1.5              # mmol O2/gX/h 维持
    c_star_mm: float = 0.25       # mmol/L 饱和 DO（空气,1bar,37C 近似）
    ki_do: float = 0.0002         # mmol/L 生长对 DO 的半饱和常数（≈0.08% 空气饱和度；正常操作点 f_DO≈1）


@dataclass
class SimRecipe:
    batch_volume_l: float = 10.0
    t_end_h: float = 24.0
    s0_g_l: float = 20.0
    x0_g_l: float = 0.3
    induction_h: float = 16.0
    temp_pre: float = 37.0
    temp_post: float = 30.0
    feed_start_h: float = 7.0
    mu_set: float = 0.11          # 1/h 指数补料（高密度发酵典型 0.09-0.14）
    feed_s_g_l: float = 500.0
    feed_s_g_l_post: float = 0.0   # 诱导期补料液碳源浓度 g/L；<=0 表示沿用 feed_s_g_l（实验设计按阶段注入）
    ph_set: float = 7.0
    rpm: float = 800.0
    airflow_lmin: float = 10.0    # vvm≈1
    press_bar: float = 1.0
    kla_scale: float = 1.0        # 大罐传质残余折扣（van't Riet 之上的治理因子）
    # 叶轮几何：用于按 van't Riet 由 P/V 计算 kLa，与 Scale-up Twin 共用同一模型
    di_m: float = 0.08            # 叶轮直径 m（10L 罐 rushton 典型）
    n_imp: int = 2
    impeller_type: str = "rushton"
    rho: float = 1030.0           # 发酵液密度 kg/m³
    medium: str | None = None     # 培养基预设键（见 fip.twins.process.media）；命中后覆盖机理参数
    ph_curve: list | None = None   # 分段 pH：[[t_h, value], ...]（按 t 升序）；缺省退化为 ph_set（供实验设计阶段注入）
    temp_curve: list | None = None # 分段温度：[[t_h, value], ...]；缺省退化为 temp_pre/temp_post
    # 阶段 H · 表达负担耦合
    ib_frac: float = 0.3           # 包涵体占产物比例（0..1）
    burden_sol: float = 0.0        # 可溶产物合成对生长的负担系数 λsol (gX/gP)
    burden_ib: float = 0.0         # 包涵体合成对生长的负担系数 λib (gX/gP)
    tox_k: float = 0.0             # 可溶（toxic）产物 μ 抑制常数 1/(g/L)
    expression_form: str | None = None  # 表达形式（intracellular_soluble/inclusion_body/periplasmic/secreted）；设定后覆盖上述负担字段
    # 阶段 FBA 耦合（将 FBA 代谢 Oracle 灌入动态仿真；默认 False 保持回归安全）
    carbon_source: str = "glucose"       # 动态仿真使用的碳源（glucose/glycerol），供 FBA 摄取/溢出 Oracle
    strain_factor: float = 1.0            # 菌株改造因子（野生型=1.0），作用于 FBA 摄取 Oracle
    protein_mw_kda: float = 30.0          # 异源蛋白分子量 kDa（FBA 负担标定记录用；μ_het 主要由 level×form 决定）
    expression_level: float = 0.15       # 异源蛋白占细胞干重比例（FBA 负担 f = level×form_factor）
    fba_coupled: bool = False             # 用 FBA Oracle 替代内生 mu_max/q_s_max/溢出系数（关掉两层各算各的）
    fba_temp: float = 37.0                # FBA Oracle 的代表性工艺温度（°C）


@dataclass
class BatchTrajectory:
    t: np.ndarray = field(default_factory=lambda: np.array([]))
    X: np.ndarray = field(default_factory=lambda: np.array([]))    # g/L DCW
    S: np.ndarray = field(default_factory=lambda: np.array([]))    # g/L 葡萄糖
    P: np.ndarray = field(default_factory=lambda: np.array([]))    # g/L 产物
    Psol: np.ndarray = field(default_factory=lambda: np.array([])) # g/L 可溶产物
    Pib: np.ndarray = field(default_factory=lambda: np.array([]))  # g/L 包涵体
    A: np.ndarray = field(default_factory=lambda: np.array([]))    # g/L 乙酸
    DO: np.ndarray = field(default_factory=lambda: np.array([]))   # %
    OUR: np.ndarray = field(default_factory=lambda: np.array([]))  # mmol/L/h
    CER: np.ndarray = field(default_factory=lambda: np.array([]))
    RQ: np.ndarray = field(default_factory=lambda: np.array([]))
    feed: np.ndarray = field(default_factory=lambda: np.array([]))  # L/h
    rpm: np.ndarray = field(default_factory=lambda: np.array([]))
    air: np.ndarray = field(default_factory=lambda: np.array([]))
    o2_enrich: np.ndarray = field(default_factory=lambda: np.array([]))
    temp: np.ndarray = field(default_factory=lambda: np.array([]))
    ph: np.ndarray = field(default_factory=lambda: np.array([]))
    induced: np.ndarray = field(default_factory=lambda: np.array([]))


def _interp_curve(curve, t, default):
    """分段曲线插值：curve 为 [[t_h, value], ...]（按 t 升序）；为空或越界返回 default。"""
    if not curve:
        return default
    if t <= curve[0][0]:
        return float(curve[0][1])
    for j in range(len(curve) - 1):
        t0, v0 = curve[j]; t1, v1 = curve[j + 1]
        if t0 <= t <= t1:
            if t1 == t0:
                return float(v1)
            return float(v0 + (v1 - v0) * (t - t0) / (t1 - t0))
    return float(curve[-1][1])


def simulate_ecoli(
    recipe: SimRecipe,
    params: EColiParams | None = None,
    dt_h: float = 0.1,
    noise: float = 0.0,
    rng: np.random.Generator | None = None,
) -> BatchTrajectory:
    p = params or EColiParams()
    rc = recipe
    # 阶段 H · 表达负担耦合：解析有效负担参数（表达形式为预设主开关）
    ib_frac, burden_sol, burden_ib, tox_k = resolve_expression_burden(rc)
    # 培养基预设：配方字段仅在仍为默认值（用户/UI 未显式设定）时应用；机理参数始终覆盖
    _rc_def = SimRecipe()
    for _k, _v in recipe_for_medium("ecoli", rc.medium).items():
        if hasattr(rc, _k) and getattr(rc, _k) == getattr(_rc_def, _k, None):
            setattr(rc, _k, _v)
    # 培养基预设：命中则用预设机理参数覆盖（如复合培养基 ki_acetate↑ 缓解溢流乙酸）
    for _k, _v in params_for_medium("ecoli", rc.medium).items():
        setattr(p, _k, _v)

    # 阶段 FBA 耦合：用 FBA 代谢 Oracle 校准内禀参数（仅 fba_coupled；默认 None 保持回归安全）
    fba = _fba_couple(rc, p) if rc.fba_coupled else None
    if fba is not None:
        print(f"[FBA-coupled] carbon={fba['carbon']} T={fba['T']} strain={fba['strain_factor']}: "
              f"mu_ceiling={fba['mu_ceiling']:.3f}/h, q_s_max={fba['q_s_max']:.3f} gS/gX/h, "
              f"burden drop_frac={fba['drop_frac']:.3f}, y_ac={fba['y_ac_gg']:.4f} g/g")

    rng = rng or np.random.default_rng(42)
    n = int(rc.t_end_h / dt_h) + 1
    tr = BatchTrajectory()
    tr.t = np.arange(n) * dt_h

    X, S, P, A, Psol, Pib = rc.x0_g_l, rc.s0_g_l, 0.0, 0.0, 0.0, 0.0
    vol = rc.batch_volume_l

    arrs = {k: np.zeros(n) for k in
            ("X", "S", "P", "Psol", "Pib", "A", "DO", "OUR", "CER", "RQ", "feed",
             "rpm", "air", "o2_enrich", "temp", "ph", "induced")}

    vvm = rc.airflow_lmin / rc.batch_volume_l

    for i, t in enumerate(tr.t):
        induced = t >= rc.induction_h
        temp = _interp_curve(rc.temp_curve, t, rc.temp_post if induced else rc.temp_pre)
        ph = _interp_curve(rc.ph_curve, t, rc.ph_set)
        # 对称惩罚：高温（热应激）与低温（代谢减慢）两侧都应降低 μ；37 ℃ ±1 ℃ 视为最优区
        # FBA 耦合时改用 FBA 的 cardinal 温度因子（与 μ_cap(T) 一致，非纯 exp 近似）
        if fba is not None:
            f_temp = fba["gtf"](temp)
        else:
            f_temp = 1.0 if abs(temp - 37.0) < 1.0 else math.exp(-0.06 * abs(37.0 - temp))
        f_ph = 1.0 if abs(ph - rc.ph_set) < 0.3 else math.exp(-0.25 * abs(ph - rc.ph_set))

        # ---- 补料（g/L/h 葡萄糖） ----
        # 诱导期可切换补料液浓度（实验设计 Induction 阶段碳源浓度）；质量流率由 mu_set 决定，
        # 浓度只改变体积流率 → 影响稀释与最终装液量（浓度越高、稀释越小）
        feed_c = rc.feed_s_g_l_post if (induced and rc.feed_s_g_l_post > 0) else rc.feed_s_g_l
        if t < rc.feed_start_h or X < 0.5:
            F_lh = 0.0
        else:
            F_lh = (rc.mu_set * X * vol / p.yxs) / feed_c
        feed_g = F_lh * feed_c / vol

        # ---- 潜在速率（不含 DO 限制） ----
        inhibition = 1.0 / (1.0 + A / p.ki_acetate)
        monod = S / (p.ks + S) if S > 0 else 0.0
        mu_pot_raw = p.mu_max * f_temp * f_ph * monod * inhibition    # 1/h

        # ---- 碳平衡（先用未限 DO 的潜在 μ 估 realized OUR，供 DO 级联） ----
        avail = S / dt_h + feed_g                                # g/L/h 可用
        uptake_pot = p.q_s_max * f_temp * monod * X              # g/L/h
        uptake = min(uptake_pot, avail)
        growth_raw = min(mu_pot_raw * X, uptake * p.yxs)         # g/L/h（无 DO 限制）

        # ---- DO 限制（X 生长耦合 DO）：用 realized OUR 估 CL，再 Monod 限生长 ----
        # 级联用与改动前一致的 realized OUR，正常操作点 DO 处于 25% 以上、f_DO≈1；
        # 仅当 OUR 超过供氧能力（真正氧限制）时 DO 下探、μ 被压低。
        our_tent = X * (p.qo2_growth * (growth_raw / X) + p.mo2)
        rpm, air, o2_frac = rc.rpm, rc.airflow_lmin, 0.21
        kla_last = c_star_last = 0.0
        for _ in range(16):  # 级联迭代（足够收敛到执行器上限，使 DO 在供氧能力内被维持）
            pv_w = power_per_volume_w_m3(rpm, rc.di_m, rc.impeller_type, vol, rc.n_imp, rc.rho)
            kla_last = rc.kla_scale * kla_h(pv_w, vvm)  # van't Riet，与 Scale-up Twin 一致
            c_star_last = p.c_star_mm * (o2_frac / 0.21) * rc.press_bar
            c_l = max(c_star_last - our_tent / kla_last, 0.0)
            do_sat = min(max(c_l / c_star_last * 100.0, 0.0), 100.0)
            if do_sat >= 25.0 or rpm >= rc.rpm * 1.6:
                break
            rpm = min(rpm * 1.15, rc.rpm * 1.6)
            o2_frac = min(o2_frac + 0.12, 0.95)
        f_do = c_l / (p.ki_do + c_l)   # DO 限制项（c_l 单位 mmol/L；ki_do 同单位）
        mu_pot = mu_pot_raw * f_do     # 氧不足时 μ 被压低 → X 与 DO 耦合

        uptake_pot = p.q_s_max * f_temp * monod * X              # g/L/h

        # ---- 碳平衡约束：摄取不能超过（存量+补料） ----
        avail = S / dt_h + feed_g                                # g/L/h 可用
        uptake = min(uptake_pot, avail)
        growth = min(mu_pot * X, uptake * p.yxs)                 # g/L/h
        # 溢流：摄取超过生长需求 -> 乙酸
        # FBA 耦合时乙酸产率改用 FBA 氧化还原 yield（g/g），否则用内生现象学系数 0.7
        overflow_carbon = max(0.0, uptake - growth / p.yxs)
        acetate_gen = overflow_carbon * (fba["y_ac_gg"] if fba is not None else 0.7)
        acetate_cons = 0.08 * A * X / (1.0 + S)
        q_p = p.q_p_max * (1.0 if induced else 0.0) * inhibition * f_ph

        # 阶段 H · 产物分叉（可溶/包涵体）+ 表达负担 + 可溶毒性（额外维持能耗）对生长的耦合
        # FBA 耦合时：表达负担用 FBA μ_het 降幅 drop_frac 直接压生长（关掉内生 burden_eff*q_p）；
        # 可溶毒性 tox_k 仍独立保留（FBA 不建模毒性死亡）。产物分叉 ib_frac 两种模式都用。
        d_psol = (1.0 - ib_frac) * q_p * X
        d_pib = ib_frac * q_p * X
        if fba is not None:
            dX = growth * fba["drop_frac"] - 0.01 * X - tox_k * Psol * X
        else:
            burden_eff = burden_sol * (1.0 - ib_frac) + burden_ib * ib_frac
            dX = growth - 0.01 * X - burden_eff * q_p * X - tox_k * Psol * X
        dS = feed_g - uptake
        dA = acetate_gen - acetate_cons
        dP = d_psol + d_pib

        X = max(X + dX * dt_h, 0.0)
        S = max(S + dS * dt_h, 0.0)
        A = max(A + dA * dt_h, 0.0)
        P = max(P + dP * dt_h, 0.0)
        Psol = max(Psol + d_psol * dt_h, 0.0)
        Pib = max(Pib + d_pib * dt_h, 0.0)
        vol += F_lh * dt_h

        # ---- 呼吸：OUR = X·(qo2_g·μ_eff + mO2) ----
        mu_eff = growth / X if X > 1e-9 else 0.0
        our = X * (p.qo2_growth * mu_eff + p.mo2)
        # 乙酸同化略降 CER；葡萄糖呼吸 RQ≈1.05
        cer = max(our * 1.05 - 0.5 * acetate_cons, 0.0)
        rq = cer / our if our > 1e-6 else 1.0

        # ---- DO 最终值：用实生长 OUR 与已确定执行器重算（≥ 级联 tentative） ----
        c_l_final = max(c_star_last - our / kla_last, 0.0) if kla_last > 0 else 0.0
        do_sat = (min(max(c_l_final / c_star_last * 100.0, 0.0), 100.0)
                  if c_star_last > 0 else 0.0)

        def nz(v: float, rel: float | None = None) -> float:
            r = noise if rel is None else rel
            return float(v * (1 + rng.normal(0, r))) if r > 0 else float(v)

        arrs["X"][i], arrs["S"][i], arrs["P"][i], arrs["Psol"][i], arrs["Pib"][i], arrs["A"][i] = X, S, P, Psol, Pib, A
        arrs["DO"][i] = max(min(nz(do_sat, noise * 0.3), 100.0), 0.0)
        arrs["OUR"][i], arrs["CER"][i], arrs["RQ"][i] = nz(our), nz(cer), nz(rq, noise * 0.5)
        arrs["feed"][i], arrs["rpm"][i] = F_lh, rpm
        arrs["air"][i], arrs["o2_enrich"][i] = air, o2_frac
        arrs["temp"][i], arrs["ph"][i] = nz(temp), nz(ph)
        arrs["induced"][i] = 1.0 if induced else 0.0

    for k, v in arrs.items():
        setattr(tr, k, v)
    return tr


def od600_from_dcw(dcw_g_l: float, host: str = "ecoli") -> float:
    """OD600 ≈ DCW/0.35（E. coli 线性近似）。"""
    return dcw_g_l / 0.35


# ===========================================================================
# 毕赤酵母（Pichia pastoris）高密度补料分批：甘油生长 + 组成型/诱导型表达
# 与 scripts/ferment_sim.js 的 simulatePichia 逐位一致；FBA 耦合用
# scripts/calibrate_burden_pichia 的 iPP668 平行 biomass Oracle（drop ≡1/(1+f)）。
# ===========================================================================
C_STAR_AIR = 0.25   # mmol/L 饱和 DO（空气, 1bar, ~28-37°C 近似）；与 E. coli 同基准


@dataclass
class PichiaParams:
    # 文献锚定（FBA 耦合时 μ_max 被替换为 MU_REF_P，Yxs 对齐 FBA 参考值）
    mu_max_gly: float = 0.28        # 1/h 甘油野生型（实现值；FBA 时改 MU_REF_P_GLY=0.20）
    mu_max_meoh: float = 0.14       # 1/h 甲醇野生型（FBA 时改 MU_REF_P_MEOM=0.12）
    ks_gly: float = 0.08            # g/L
    ks_meoh: float = 0.25           # g/L
    yxs_gly: float = 0.50           # gX/gS
    yxs_meoh: float = 0.42          # gX/gS
    q_p_max: float = 0.004          # gP/gX/h
    q_meoh_maint: float = 0.03      # 1/h AOX 维持甲醇耗
    ki_meoh: float = 4.0            # g/L 甲醇抑制常数
    qo2_growth: float = 10.0        # mmol O2/gX 生长相关
    mo2: float = 1.2                # mmol O2/gX/h 维持
    o2_per_g_meoh: float = 46.9     # mmol O2 / g 甲醇（甲醇完全氧化高耗氧）
    c_star_mm: float = C_STAR_AIR


@dataclass
class PichiaRecipe:
    host: str = "pichia"
    batch_volume_l: float = 10.0
    t_end_h: float = 96.0
    s0_g_l: float = 20.0            # 甘油批初始 g/L
    x0_g_l: float = 0.5
    gly_feed_start_h: float = 18.0
    gly_feed_end_h: float = 32.0
    mu_set: float = 0.08            # 1/h 甘油补料设定
    feed_s_g_l: float = 500.0       # 甘油补料液 g/L
    meoh_start_h: float = 34.0
    meoh_feed_rate_lh: float = 0.012  # L/h 甲醇补料
    meoh_density: float = 0.792
    temp_pre: float = 28.0
    ph_set: float = 5.5
    rpm: float = 900.0
    airflow_lmin: float = 10.0
    press_bar: float = 1.0
    kla_scale: float = 1.0
    di_m: float = 0.08
    n_imp: int = 2
    impeller_type: str = "rushton"
    rho: float = 1030.0
    medium: str | None = None
    ph_curve: list | None = None
    temp_curve: list | None = None
    # 阶段 H · 表达负担耦合（与 E. coli 同语义）
    ib_frac: float = 0.3
    burden_sol: float = 0.0
    burden_ib: float = 0.0
    tox_k: float = 0.0
    expression_form: str | None = None
    # 表达模式 / 诱导型补料模式
    expression_mode: str = "inducible"   # constitutive / inducible
    meoh_feed_mode: str = "methanol"     # methanol（纯甲醇）/ mixed（甲醇+甘油混合共利用）
    mixed_gly_feed_g_l: float = 200.0    # 混合补料甘油浓度 g/L
    mixed_gly_start_h: float = 0.0       # 0 ⇒ 沿用 meoh_start_h
    mixed_gly_end_h: float = 0.0         # 0 ⇒ 沿用 t_end_h
    # 阶段 FBA 耦合
    strain_factor: float = 1.0
    protein_mw_kda: float = 30.0
    expression_level: float = 0.15
    fba_coupled: bool = False
    fba_temp: float = 28.0


@dataclass
class PichiaTrajectory:
    t: np.ndarray = field(default_factory=lambda: np.array([]))
    X: np.ndarray = field(default_factory=lambda: np.array([]))    # g/L DCW
    S: np.ndarray = field(default_factory=lambda: np.array([]))    # g/L 甘油
    P: np.ndarray = field(default_factory=lambda: np.array([]))    # g/L 产物
    M: np.ndarray = field(default_factory=lambda: np.array([]))    # g/L 甲醇
    Psol: np.ndarray = field(default_factory=lambda: np.array([]))
    Pib: np.ndarray = field(default_factory=lambda: np.array([]))
    DO: np.ndarray = field(default_factory=lambda: np.array([]))
    OUR: np.ndarray = field(default_factory=lambda: np.array([]))
    CER: np.ndarray = field(default_factory=lambda: np.array([]))
    RQ: np.ndarray = field(default_factory=lambda: np.array([]))
    feed: np.ndarray = field(default_factory=lambda: np.array([]))
    rpm: np.ndarray = field(default_factory=lambda: np.array([]))
    air: np.ndarray = field(default_factory=lambda: np.array([]))
    temp: np.ndarray = field(default_factory=lambda: np.array([]))
    ph: np.ndarray = field(default_factory=lambda: np.array([]))
    induced: np.ndarray = field(default_factory=lambda: np.array([]))
    aox: np.ndarray = field(default_factory=lambda: np.array([]))   # 醇氧化酶诱导水平 0..1
    mode: np.ndarray = field(default_factory=lambda: np.array([]))  # 操作阶段标注


def _fba_couple_pichia(rc: "PichiaRecipe", p: "PichiaParams") -> dict:
    """阶段 FBA 耦合（Pichia）：用 calibrate_burden_pichia 的 iPP668 平行 biomass Oracle 校准。
    返回 drop_frac(FBA μ_het/μ_wt 降幅)、gtf(Pichia cardinal 温度因子)、feasible，
    并就地改写 p.mu_max_gly/mu_max_meoh=MU_REF_P、p.yxs=YXS_P。"""
    from scripts.calibrate_burden_pichia import (  # noqa: E402
        metrics_pichia, effective_uptake_pichia, growth_temperature_factor_pichia,
        mu_cap_pichia, MU_REF_P_GLY, MU_REF_P_MEOM, YXS_P, FORM_FACTOR_PICHIA,
    )
    carbon = "methanol" if rc.expression_mode == "inducible" else "glycerol"
    T = rc.fba_temp
    sf = rc.strain_factor
    form = rc.expression_form
    ff = FORM_FACTOR_PICHIA.get(form, 1.0)
    f = rc.expression_level * ff
    sub_max = effective_uptake_pichia(carbon, T=T, strain_factor=sf)
    m = metrics_pichia(carbon, sub_max, 30.0, f=f, T=T)
    mu_wt = m.get("mu_wt") or 0.0
    mu_het = m.get("mu_het")
    drop_frac = (mu_het / mu_wt) if (mu_wt > 0 and mu_het is not None) else 1.0
    # 内禀参数对齐到 FBA 物理基准（动态仿真每步再乘 cardinal 温度因子 -> μ_cap(步温)）
    p.mu_max_gly = MU_REF_P_GLY
    p.mu_max_meoh = MU_REF_P_MEOM
    p.yxs_gly = YXS_P["glycerol"]
    p.yxs_meoh = YXS_P["methanol"]
    return {
        "carbon": carbon, "T": T, "strain_factor": sf,
        "drop_frac": drop_frac, "mu_ceiling": mu_cap_pichia(T, carbon),
        "feasible": bool(m.get("feasible", False)),
        "gtf": growth_temperature_factor_pichia,   # Pichia cardinal 温度因子
    }


def simulate_pichia(
    recipe: PichiaRecipe,
    params: PichiaParams | None = None,
    dt_h: float = 0.25,
    noise: float = 0.0,
    rng: np.random.Generator | None = None,
) -> PichiaTrajectory:
    p = params or PichiaParams()
    rc = recipe
    from scripts.calibrate_burden_pichia import growth_temperature_factor_pichia  # noqa: E402
    ib_frac, burden_sol, burden_ib, tox_k = resolve_expression_burden(rc)
    # 培养基预设（配方仅在仍为默认值时覆盖；机理参数始终覆盖）
    if rc.medium:
        _def = PichiaRecipe()
        for _k, _v in recipe_for_medium("pichia", rc.medium).items():
            if hasattr(rc, _k) and getattr(rc, _k) == getattr(_def, _k, None):
                setattr(rc, _k, _v)
        for _k, _v in params_for_medium("pichia", rc.medium).items():
            setattr(p, _k, _v)

    expr_mode = rc.expression_mode or "inducible"
    meoh_mode = rc.meoh_feed_mode or "methanol"

    fba = _fba_couple_pichia(rc, p) if rc.fba_coupled else None
    FBA_DROP = 1.0
    if fba is not None:
        FBA_DROP = fba["drop_frac"]
        print(f"[FBA-coupled:Pichia] carbon={fba['carbon']} mode={expr_mode} T={fba['T']} "
              f"strain={fba['strain_factor']}: drop_frac={FBA_DROP:.3f}, "
              f"mu_ceiling={fba['mu_ceiling']:.3f}/h, feasible={fba['feasible']}")

    rng = rng or np.random.default_rng(7)
    n = int(rc.t_end_h / dt_h) + 1
    tr = PichiaTrajectory()
    tr.t = np.arange(n) * dt_h

    X = rc.x0_g_l
    S = rc.s0_g_l
    M = 0.0
    P = Psol = Pib = 0.0
    vol = rc.batch_volume_l
    aox = 0.0

    arrs = {k: np.zeros(n) for k in
            ("X", "S", "P", "M", "Psol", "Pib", "DO", "OUR", "CER", "RQ",
             "feed", "rpm", "air", "temp", "ph", "induced", "aox", "mode")}

    vvm = rc.airflow_lmin / rc.batch_volume_l

    for i, t in enumerate(tr.t):
        on_meoh = (expr_mode == "inducible") and (t >= rc.meoh_start_h)
        temp = _interp_curve(rc.temp_curve, t, rc.temp_pre)
        ph = _interp_curve(rc.ph_curve, t, rc.ph_set)
        # Pichia cardinal 温度因子（FBA 与否同基准）；T 出温区 -> 0
        f_temp = growth_temperature_factor_pichia(temp) if fba is not None else _pichia_temp_factor(temp)
        f_ph = 1.0 if abs(ph - rc.ph_set) < 0.3 else math.exp(-0.30 * abs(ph - rc.ph_set))
        if on_meoh and M > 0.1:
            aox = min(aox + 0.12 * dt_h, 1.0)
        elif not on_meoh:
            aox = max(aox - 0.05 * dt_h, 0.0)

        # 甘油补料：组成型全程；诱导型仅 meoh 前
        gly_feed_active = (t >= rc.gly_feed_start_h and t < rc.gly_feed_end_h and X > 0.5)
        f_gly = (rc.mu_set * X * vol / p.yxs_gly) / rc.feed_s_g_l if gly_feed_active else 0.0
        feed_gly_g = f_gly * rc.feed_s_g_l / vol

        # 甲醇补料（仅诱导型）
        f_meoh = rc.meoh_feed_rate_lh if on_meoh else 0.0
        feed_meoh_g = f_meoh * rc.meoh_density * 1000.0 / vol

        # 诱导型 · 甲醇+甘油混合补料（共利用）
        if expr_mode == "inducible" and meoh_mode == "mixed" and on_meoh:
            mg_start = rc.mixed_gly_start_h if rc.mixed_gly_start_h > 0 else rc.meoh_start_h
            mg_end = rc.mixed_gly_end_h if rc.mixed_gly_end_h > 0 else rc.t_end_h
            if t >= mg_start and t < mg_end and X > 0.5:
                f_gly_mix = (rc.mu_set * 0.6 * X * vol / p.yxs_gly) / rc.mixed_gly_feed_g_l
                feed_gly_g += f_gly_mix * rc.mixed_gly_feed_g_l / vol

        monod_g = S / (p.ks_gly + S) if S > 0 else 0.0
        avail_g = S / dt_h + feed_gly_g
        growth_g = min(p.mu_max_gly * f_temp * f_ph * monod_g * X, avail_g * p.yxs_gly)
        uptake_g = growth_g / p.yxs_gly + 0.015 * X

        meoh_inhib = 1.0 / (1.0 + M / p.ki_meoh)
        monod_m = M / (p.ks_meoh + M) if M > 0 else 0.0
        avail_m = M / dt_h + feed_meoh_g
        growth_m = min(p.mu_max_meoh * f_temp * f_ph * aox * monod_m * meoh_inhib * X, avail_m * p.yxs_meoh)
        maint_meoh = p.q_meoh_maint * aox * X if on_meoh else 0.0
        cons_meoh = min(growth_m / p.yxs_meoh + maint_meoh, avail_m)

        growth = growth_g + growth_m
        expr_active = 1.0 if (expr_mode == "constitutive" or on_meoh) else 0.0
        q_p = p.q_p_max * f_temp * f_ph * (aox if expr_mode == "inducible" else 1.0) * expr_active

        # 阶段 H · 产物分叉 + 表达负担 + 可溶毒性
        burden_eff = burden_sol * (1.0 - ib_frac) + burden_ib * ib_frac
        d_psol = (1.0 - ib_frac) * q_p * X
        d_pib = ib_frac * q_p * X
        # FBA 耦合：表达负担用 FBA μ_het 降幅 drop_frac 直接压表达相生长（关掉内生负担）；
        # 甲醇 DO 受限溢出（甲醛毒性）由动态层 ki_meoh 处理（FBA 表 yac=0）。
        growth_eff = growth * FBA_DROP if (fba is not None and expr_active) else growth
        dX = growth_eff - 0.008 * X - (0.0 if fba is not None else burden_eff * q_p * X) - tox_k * Psol * X
        dS = feed_gly_g - uptake_g
        dM = feed_meoh_g - cons_meoh
        dP = d_psol + d_pib

        X = max(X + dX * dt_h, 0.0)
        S = max(S + dS * dt_h, 0.0)
        M = max(M + dM * dt_h, 0.0)
        P = max(P + dP * dt_h, 0.0)
        Psol = max(Psol + d_psol * dt_h, 0.0)
        Pib = max(Pib + d_pib * dt_h, 0.0)
        vol += (f_gly + f_meoh) * dt_h

        mu_eff = growth / X if X > 1e-9 else 0.0
        our = X * (p.qo2_growth * mu_eff + p.mo2)
        our += (cons_meoh - growth_m / p.yxs_meoh) * p.o2_per_g_meoh
        rq = (0.75 - 0.05 * aox) if on_meoh else 1.0
        cer = max(our * rq, 0.0)

        rpm, air, o2_frac = rc.rpm, rc.airflow_lmin, 0.21
        do_sat = 100.0
        for _ in range(3):
            pv_w = power_per_volume_w_m3(rpm, rc.di_m, rc.impeller_type, vol, rc.n_imp, rc.rho)
            kla = rc.kla_scale * kla_h(pv_w, vvm)
            c_star = p.c_star_mm * (o2_frac / 0.21) * rc.press_bar
            c_l = max(c_star - our / kla, 0.0)
            do_sat = min(max(c_l / c_star * 100.0, 0.0), 100.0)
            if do_sat >= 25.0 or rpm >= rc.rpm * 1.6:
                break
            rpm = min(rpm * 1.2, rc.rpm * 1.6)
            o2_frac = min(o2_frac + 0.15, 0.95)

        # 操作阶段标注（mode 数组为数值编码，与 JS arrs.MODE 一致）
        if expr_mode == "constitutive":
            phase = "gly_bat" if t < rc.gly_feed_start_h else "gly_feed"
        elif t < rc.meoh_start_h:
            phase = "gly_feed" if t < rc.gly_feed_end_h else "starv"
        else:
            phase = "mix_feed" if meoh_mode == "mixed" else "meoh_induce"
        _PHASE_CODE = {"gly_bat": 0.0, "gly_feed": 1.0, "starv": 2.0,
                       "meoh_induce": 3.0, "mix_feed": 4.0}
        phase_code = _PHASE_CODE.get(phase, 1.0)

        def nz(v: float, rel: float | None = None) -> float:
            r = noise if rel is None else rel
            return float(v * (1 + rng.normal(0, r))) if r > 0 else float(v)

        arrs["X"][i], arrs["S"][i], arrs["P"][i] = X, S, P
        arrs["M"][i], arrs["Psol"][i], arrs["Pib"][i] = M, Psol, Pib
        arrs["DO"][i] = max(min(nz(do_sat, noise * 0.3), 100.0), 0.0)
        arrs["OUR"][i], arrs["CER"][i], arrs["RQ"][i] = nz(our), nz(cer, noise * 0.5), nz(rq, noise * 0.5)
        arrs["feed"][i], arrs["rpm"][i], arrs["air"][i] = f_meoh, rpm, air
        arrs["temp"][i], arrs["ph"][i] = nz(temp), nz(ph)
        arrs["induced"][i], arrs["aox"][i], arrs["mode"][i] = (1.0 if on_meoh else 0.0), aox, phase_code

    for k, v in arrs.items():
        setattr(tr, k, v)
    return tr


def _pichia_temp_factor(T: float) -> float:
    """非 FBA 耦合时的 Pichia 温度因子（对称惩罚，±2°C 内等效，超出 exp 衰减）。"""
    if abs(T - 28.0) <= 2.0:
        return 1.0
    return math.exp(-0.10 * (abs(T - 28.0) - 2.0))


def simulate(recipe, params=None, dt_h=0.25, noise=0.0, rng=None):
    """host 派发：PichiaRecipe -> simulate_pichia；SimRecipe -> simulate_ecoli。"""
    if isinstance(recipe, PichiaRecipe):
        return simulate_pichia(recipe, params, dt_h, noise, rng)
    return simulate_ecoli(recipe, params, dt_h, noise, rng)
