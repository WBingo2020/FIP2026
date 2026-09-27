"""
uptake_capacity.py — E. coli 碳源摄取速率（文献基准 + 菌株改造 × 温度修正）

FIP 工艺乙酸溢出模块需要"宿主碳源摄取能力"作为工艺约束。此前该值（glucose=10、
glycerol=8）是手设的物理合理默认、无文献出处。本模块把它换成**有文献依据的基准值**，
并支持两类用户明确要求的变异：
  (1) 菌株宿主被改造 -> 摄取速率变化（strain_factor，野生型=1）；
  (2) 不同温度条件 -> 摄取速率差异（Q10 / Arrhenius，参考温度 37°C）。

此外，用户要求**生长率本身也随温度变**（μ~Q10）。此前只把温度修正作用在"摄取速率"；
现扩展为细胞内禀最大生长率上限 μ_max(T) 也受温度控制（见文末"生长率温度依赖"段），
并作为 FBA 求解的上界（model.objective upper_bound = μ_cap(T)），使乙酸/μ 通量分配与
受限后的 μ 自洽。

有效摄取速率：
    v_eff = v_lit[carbon] × strain_factor[carbon] × Q10[carbon] ^ ((T - T_ref)/10)

--- 文献基准 v_lit (mmol/gDW/h, 37°C, 野生型) -----------------------------------------
葡萄糖 glucose:
  * 10.5  Varma & Palsson (1994) 好氧最大葡萄糖摄取 (被 DySEEP 2026 / Processes 14:1659 转引:
           "Varma (1994) experimentally determined the maximum glucose and oxygen uptake
            rates for E. coli under aerobic conditions as 10.5 and 15 mmol/CDW·h")
  * 15    过量葡萄糖 batch 的"典型生理"值 (PNAS 118:e2013836118 dETFL: "initial uptake rate
            of glucose is set to 15 mmol·gDW⁻¹·h⁻¹ ... characteristic of a typical physiology
            for E. coli growing on glucose with excess oxygen")
  * 6.5–8.5  葡萄糖受限稳态实测 (Wiley, Metabolite flux profiling, 3 株野生型 MG1655/W3110/JM101)
  => 取最大能力 10.5 作为"宿主摄取能力"基准（乙酸溢出由补料超此能力触发）。

甘油 glycerol:
  * 1.5–2.4  GlpF/GlpK 转运 Vmax 25–40 nmol/min/mgDW (Metabeng 引 Orjuela et al. 2020,
           "Glycerol GlpF/GlpK 10–20 µM Km, 25–40 nmol/min/mg dw") -> 换算 1.5–2.4 mmol/gDW/h
  * 甘油是 E. coli 弱碳源，摄取显著慢于葡萄糖（约 1/4–1/5），故基准取 2.5。
  （注：某 PLOS ONE 动力学拟合给 vmaxGly=174.86 mmol/g·h，量纲/模型特异性存疑，未采用。）

--- 温度 Q10 (参考 T_ref = 37°C) ----------------------------------------------------------
  * 通用生物/酶促速率 Q10≈2 (Ea≈50 kJ/mol)；E. coli 生长 17→27→37°C 每段 ~×2.5 (有效 Ea≈60 kJ/mol,
    Cell Biology by the Numbers, BNID 100919)。
  * 葡萄糖 PTS 为酶促转运 -> Q10_glc = 2.0。
  * 甘油：GlpF 孔道为"孔扩散"机制、Ea 仅 4.5 kcal/mol (≈Q10 1.2, Maurel 1994 JBC)，但胞内
    甘油激酶 GLYK 为酶促 (Ea≈50 kJ/mol) 且常限速 -> 整体 Q10_gly = 1.9（略低于纯酶促）。

--- 菌株改造因子 strain_factor (野生型=1.0) --------------------------------------------
  * ↑ 过表达 PTS / GlpF-GlpK：葡萄糖可推至 15–20 (Orders of Magnitude, masspy)；甘油上调同理。
  * ↓ 敲除 pgi / pfkA：葡萄糖摄取降至 2.5 / 1.4 (Wiley flux profiling)；其他通路改造同理。
  * 区间约 0.2–2.0+，默认 1.0 = 未改造野生型。

用法：
    from scripts.uptake_capacity import effective_uptake
    v = effective_uptake("glucose", T=30.0, strain_factor=1.2)   # 30°C, PTS 过表达 1.2×
    v = effective_uptake("glycerol", T=37.0, strain_factor=0.5)  # 甘油激酶弱化株
"""
from __future__ import annotations

# 参考温度（E. coli 标准好氧生长温度）
T_REF = 37.0

# 文献基准摄取速率 (mmol/gDW/h, T_REF, 野生型)
LIT_UPTAKE = {
    "glucose": 10.5,   # Varma & Palsson 1994 (aerobic max)
    "glycerol": 2.5,   # GlpF/GlpK Vmax 1.5–2.4 mmol/gDW/h (Orjuela 2020 via Metabeng)
}

# 温度 Q10（参考 T_REF）
Q10 = {
    "glucose": 2.0,    # PTS 酶促, Ea≈50 kJ/mol
    "glycerol": 1.9,   # GLYK 限速; GlpF 孔道近温度不敏感 (Ea 4.5 kcal/mol)
}

# 文献出处（透明可查）
SOURCES = {
    "glucose": {
        "v_lit": 10.5,
        "ref": "Varma & Palsson (1994), aerobic max glucose uptake; recast in Processes 14:1659 (2026)",
        "alt": "15 mmol/gDW/h typical excess-glucose batch (PNAS 118:e2013836118); "
               "6.5–8.5 steady-state glucose-limited (Wiley Metabolite flux profiling)",
        "q10": 2.0, "q10_ref": "Ea≈50 kJ/mol enzymatic PTS; Q10≈2 generic (Cell Biol by Numbers)",
    },
    "glycerol": {
        "v_lit": 2.5,
        "ref": "GlpF/GlpK Vmax 25–40 nmol/min/mgDW = 1.5–2.4 mmol/gDW/h (Orjuela et al. 2020, via Metabeng)",
        "alt": "glycerol weak carbon source, ~1/4–1/5 of glucose capacity",
        "q10": 1.9,
        "q10_ref": "GLYK kinase-limited (Ea≈50 kJ/mol); GlpF pore Ea 4.5 kcal/mol (Maurel 1994 JBC)",
    },
}


# --- 生长率温度依赖 (μ ~ Q10, 含 cardinal 上下界) -----------------------------------------
# 用户明确要求：生长率本身也随温度变。此前仅把温度修正作用在"摄取速率"；现扩展为细胞内禀
# 最大生长率上限 μ_max(T) 也受温度控制。
#
# 模型：低于最适 T_opt 用 Q10（生长随温升提速）；高于 T_opt 用线性衰减到 T_max（蛋白/酶热
# 失活）-> 即 cardinal 温度模型，避免纯 Q10 在 >37°C 仍上涨的不物理外推（细菌在 40°C 以上
# 生长反而下降）。
#   * Q10_mu ≈ 2.2：细菌生长 25→37°C 普遍 ≈×2（Ratkowsky 1983 平方根模型 / Arrhenius Ea≈50–70
#     kJ/mol；Cell Biology by the Numbers BNID 100919：E. coli μ(25°C)≈0.5, μ(37°C)≈1.0–1.2，
#     等效 Q10≈2–2.5）。
#   * T_opt≈40°C, T_max≈46°C, T_min≈4°C：E. coli cardinal 温度（Rosso et al. 1993；BNID 100919）。
#     低于 T_min / 高于 T_max 生长停止（μ=0）。
#   * MU_REF_37 = 1.30 /hr：e_coli_core 在 37°C、碳与氧均不受限时 FBA 求得的最大 μ（与本项目
#     实测基线 1.2917 一致），作为细胞内禀最大生长能力基准。工艺有效 μ = min(FBA_μ, μ_cap(T))。
MU_REF_37 = 1.30     # /hr, 37°C 内禀最大生长（碳/氧不限时 e_coli_core FBA）
Q10_MU = 2.2         # 生长率 Q10（<T_opt 区段）
T_MIN = 4.0          # °C, 生长停止下限
T_OPT = 40.0         # °C, 最适生长温度
T_MAX = 46.0         # °C, 生长停止上限（热失活）

GROWTH_SOURCE = ("Ratkowsky (1983) square-root / Arrhenius for microbial growth (Q10≈2–2.5, "
                 "Ea≈50–70 kJ/mol); Rosso et al. (1993) cardinal temperatures for E. coli "
                 "(Tmin≈4–8, Topt≈39–40, Tmax≈46–48 °C); Cell Biol by Numbers BNID 100919.")


def growth_temperature_factor(T: float = T_REF) -> float:
    """生长速率相对 37°C 的温度因子 μ_max(T) / μ_max(37)（cardinal + Q10）。"""
    if T <= T_MIN or T >= T_MAX:
        return 0.0
    if T <= T_OPT:
        return round(Q10_MU ** ((float(T) - T_REF) / 10.0), 4)
    # 高于最适：从 μ(T_opt) 线性衰减到 0 @ T_max
    f_opt = Q10_MU ** ((T_OPT - T_REF) / 10.0)
    return round(f_opt * (T_MAX - float(T)) / (T_MAX - T_OPT), 4)


def mu_cap(T: float = T_REF, mu_ref: float = MU_REF_37) -> float:
    """细胞内禀最大生长率上限 μ_max(T)（/hr）。工艺有效 μ 受此上界约束（FBA objective ub）。"""
    return round(mu_ref * growth_temperature_factor(T), 4)


def describe_growth(T: float = T_REF) -> dict:
    """返回生长率温度依赖的透明字典（基准/卡片参数/因子/上限/文献），便于写回 JSON 与审计。"""
    return {
        "T_ref_C": T_REF,
        "T_used_C": T,
        "mu_ref_37": MU_REF_37,
        "q10_mu": Q10_MU,
        "T_min": T_MIN, "T_opt": T_OPT, "T_max": T_MAX,
        "growth_temperature_factor": growth_temperature_factor(T),
        "mu_cap": mu_cap(T),
        "source": GROWTH_SOURCE,
    }


def effective_uptake(carbon: str, T: float = T_REF, strain_factor: float = 1.0) -> float:
    """有效宿主碳源摄取能力 (mmol/gDW/h)。

    carbon: "glucose" | "glycerol"
    T: 工艺温度 °C（默认 37）
    strain_factor: 菌株改造因子（1.0=野生型；>1 过表达上调，<1 敲除/弱化下调）
    """
    if carbon not in LIT_UPTAKE:
        raise ValueError(f"未知碳源 {carbon!r}，可选 {list(LIT_UPTAKE)}")
    v = LIT_UPTAKE[carbon] * float(strain_factor) * (Q10[carbon] ** ((float(T) - T_REF) / 10.0))
    return round(v, 4)


def describe(carbon: str, T: float = T_REF, strain_factor: float = 1.0) -> dict:
    """返回透明字典：文献基准 / 温度因子 / 菌株因子 / 有效值，便于写回 JSON 与审计。"""
    lit = LIT_UPTAKE[carbon]
    q = Q10[carbon]
    tfac = q ** ((T - T_REF) / 10.0)
    eff = lit * strain_factor * tfac
    return {
        "carbon": carbon,
        "T_ref_C": T_REF,
        "T_used_C": T,
        "v_lit": lit,
        "q10": q,
        "temp_factor": round(tfac, 4),
        "strain_factor": strain_factor,
        "effective_uptake": round(eff, 4),
        "source": SOURCES[carbon],
    }


if __name__ == "__main__":
    for c in ("glucose", "glycerol"):
        for T in (30.0, 37.0, 42.0):
            print(f"{c:>8} T={T:>4.0f}°C sf=1.0 -> v_eff={effective_uptake(c, T):.3f} mmol/gDW/h")
    print()
    print("PTS 过表达 1.3× @30°C:", effective_uptake("glucose", T=30.0, strain_factor=1.3))
    print("甘油激酶弱化 0.5× @37°C:", effective_uptake("glycerol", T=37.0, strain_factor=0.5))
    print()
    print("--- 生长率温度依赖 μ_max(T) ---")
    for T in (20.0, 25.0, 30.0, 37.0, 40.0, 42.0, 45.0):
        print(f"  T={T:>4.0f}°C  factor={growth_temperature_factor(T):.3f}  μ_cap={mu_cap(T):.3f} /hr")
