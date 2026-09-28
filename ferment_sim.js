/* ============================================================================
 * FermentSim · 发酵机理前向仿真引擎（纯前端、零依赖）
 * 移植自 fip/twins/process/mechanistic.py (simulate_ecoli) 与
 *        fip/twins/process/pichia_sim.py  (simulate_pichia)
 * 与 Scale-up Twin 共用同一套 van't Riet kLa / P/V 关联式，保证
 * 展示/放大的传氧链路数值一致。
 *
 * 用法：
 *   var traj = FermentSim.simulate({ host:"ecoli", rpm:800, airflow_lmin:10, ... });
 *   traj.t / traj.X / traj.P / traj.DO / traj.OUR / traj.CER ... 均为等长数组（小时轴对齐）
 *   traj.summary = { titer, dcw, acetatePeak, doMin, ourPeak, endH }
 * ==========================================================================*/
(function () {
  "use strict";

  // ---- van't Riet 型 kLa（与 Python mass_transfer.kla_h 逐位一致）----
  // K = 0.032 (1/s, P/V 用 kW/m³, vvm 用 1/min)；P/V 入参为 W/m³ = kW/m³×1000，
  // 故等价 K(1/h) = 0.032·3600 / 1000^0.4 ≈ 7.2698。
  var KLA_K = 0.032 * 3600.0 / Math.pow(1000.0, 0.4);
  var C_STAR_AIR = 0.25; // mmol/L（空气, 1bar, 37/28C 近似）

  // ---- 阶段 FBA 耦合：浏览器侧"伪耦合"所需的常量与纯数学（与 Python uptake_capacity / cardinal 逐位一致）----
  // 浏览器无法跑 FBA 求解器；drop_frac / y_ac_gg 由预计算查表 window.FBA_TABLE 提供（见 lookupFba），
  // 以下常量为同源移植，供运行期计算 μ 内禀上限 / 宿主摄取 / 生长温度因子。
  var MU_REF_37 = 1.30, Q10_MU = 2.2, T_MIN = 4.0, T_OPT = 40.0, T_MAX = 46.0, T_REF = 37.0;
  var V_LIT = { glucose: 10.5, glycerol: 2.5 };      // 文献基准摄取 mmol/gDW/h, 37°C, 野生型
  var Q10_UPTAKE = { glucose: 2.0, glycerol: 1.9 };  // 摄取温度 Q10
  var MM_CARBON = { glucose: 180, glycerol: 92 };     // 碳源摩尔质量 g/mol（mmol->g 换算）
  var FORM_FACTOR = { intracellular_soluble: 1.20, inclusion_body: 0.80, periplasmic: 1.00, secreted: 0.90, none: 1.0 };

  // cardinal 温度模型：T<=T_opt 用 Q10；T_opt<T<T_max 线性衰减到 0；越界 μ=0（热失活/停止）
  function growthTemperatureFactor(T) {
    T = (T == null) ? T_REF : T;
    if (T <= T_MIN || T >= T_MAX) return 0.0;
    if (T <= T_OPT) return Math.pow(Q10_MU, (T - T_REF) / 10.0);
    var fOpt = Math.pow(Q10_MU, (T_OPT - T_REF) / 10.0);
    return fOpt * (T_MAX - T) / (T_MAX - T_OPT);
  }
  function muCap(T) { return MU_REF_37 * growthTemperatureFactor(T); }
  function effectiveUptake(carbon, T, sf) {
    if (!V_LIT[carbon]) return 0.0;
    var s = (sf == null) ? 1.0 : sf;
    return V_LIT[carbon] * s * Math.pow(Q10_UPTAKE[carbon], (T - T_REF) / 10.0);
  }
  // 从预计算查表取 (drop_frac, y_ac_gg)；表未注入时返回 null（调用方退化为手调默认）
  function lookupFba(carbon, fbaTemp, strain, form, level) {
    var T = (typeof window !== "undefined") ? window.FBA_TABLE : null;
    if (!T || !T.data || !T.data[carbon]) return null;
    var tg = T.grid.fba_temp, sg = T.grid.strain, lg = T.grid.level;
    var fi = T.forms.indexOf(form); if (fi < 0) fi = 0;
    function bracket(arr, v) {
      if (v <= arr[0]) return [0, 0, 0];
      if (v >= arr[arr.length - 1]) return [arr.length - 1, arr.length - 1, 1];
      for (var i = 0; i < arr.length - 1; i++) {
        if (v >= arr[i] && v <= arr[i + 1]) {
          var t = (arr[i + 1] === arr[i]) ? 0 : (v - arr[i]) / (arr[i + 1] - arr[i]);
          return [i, i + 1, t];
        }
      }
      return [0, 0, 0];
    }
    var tb = bracket(tg, fbaTemp), sb = bracket(sg, strain);
    var ti0 = tb[0], ti1 = tb[1], wt = tb[2];
    var si0 = sb[0], si1 = sb[1], ws = sb[1] != null ? sb[2] : 0;
    function cell(ti, si) {
      var arr = T.data[carbon][ti][si][fi];
      if (level <= lg[0]) return arr[0];
      if (level >= lg[lg.length - 1]) return arr[arr.length - 1];
      for (var k = 0; k < lg.length - 1; k++) {
        if (level >= lg[k] && level <= lg[k + 1]) {
          var tl = (lg[k + 1] === lg[k]) ? 0 : (level - lg[k]) / (lg[k + 1] - lg[k]);
          var a = arr[k], b = arr[k + 1];
          return { drop: a.drop + (b.drop - a.drop) * tl, yac: a.yac + (b.yac - a.yac) * tl };
        }
      }
      return arr[0];
    }
    var c00 = cell(ti0, si0), c10 = cell(ti1, si0), c01 = cell(ti0, si1), c11 = cell(ti1, si1);
    return {
      drop: c00.drop * (1 - wt) * (1 - ws) + c10.drop * wt * (1 - ws) + c01.drop * (1 - wt) * ws + c11.drop * wt * ws,
      yac: c00.yac * (1 - wt) * (1 - ws) + c10.yac * wt * (1 - ws) + c01.yac * (1 - wt) * ws + c11.yac * wt * ws
    };
  }

  // ---- 毕赤酵母 Pichia 两层 FBA 耦合：宿主专属代谢常量与查表（与 Python _fba_couple_pichia 逐位一致）----
  // Pichia 生长温区与 E. coli 不同：最优 ~28°C（T_REF_P）、上限 ~37°C（T_MAX_P）、下限 ~15°C（T_MIN_P）。
  // FBA 内禀上限（文献/iPP668 参考值，待 GSM 上传后由 FBA Oracle 重写）：甘油 μ_max≈0.20/h、
  // 甲醇 μ_max≈0.12/h（@T_REF_P）。甲醇单位碳完全氧化需氧高（o2_per_g_meoh≈46.9 mmol O2/g 甲醇）。
  var MU_REF_P_GLY = 0.20, MU_REF_P_MEOM = 0.12;
  var Q10_MU_P = 2.0, T_REF_P = 28.0, T_OPT_P = 30.0, T_MAX_P = 37.0, T_MIN_P = 15.0;
  var V_LIT_P = { glycerol: 2.8, methanol: 1.6 };   // 文献基准宿主摄取 mmol/gDW/h @ T_REF_P, 野生型
  var Q10_UPTAKE_P = { glycerol: 1.9, methanol: 1.7 };
  var MM_P = { glycerol: 92.0, methanol: 32.0 };    // 碳源摩尔质量 g/mol（mmol->g 换算）
  var YXS_P = { glycerol: 0.50, methanol: 0.40 };   // 生物质得率 gX/gS（FBA biomass 得率参考）

  // cardinal 温度模型（Pichia 专属）：T<=T_min 或 T>=T_max → μ=0（停长/热失活）；T<=T_opt → Q10；T_opt<T<T_max → 线性衰减
  function growthTemperatureFactorPichia(T) {
    T = (T == null) ? T_REF_P : T;
    if (T <= T_MIN_P || T >= T_MAX_P) return 0.0;
    if (T <= T_OPT_P) return Math.pow(Q10_MU_P, (T - T_REF_P) / 10.0);
    var fOpt = Math.pow(Q10_MU_P, (T_OPT_P - T_REF_P) / 10.0);
    return fOpt * (T_MAX_P - T) / (T_MAX_P - T_OPT_P);
  }
  function muCapPichia(T, carbon) {
    var muRef = (carbon === "methanol") ? MU_REF_P_MEOM : MU_REF_P_GLY;
    return muRef * growthTemperatureFactorPichia(T);
  }
  function effectiveUptakePichia(carbon, T, sf) {
    if (!V_LIT_P[carbon]) return 0.0;
    var s = (sf == null) ? 1.0 : sf;
    return V_LIT_P[carbon] * s * Math.pow(Q10_UPTAKE_P[carbon], (T - T_REF_P) / 10.0);
  }
  // 毕赤酵母 FBA 查表：读 window.FBA_TABLE.data_pichia[carbon]（碳源 glycerol/methanol）。
  // 表未注入（GSM 未生成）时返回 null，simulatePichia 退化为手调默认（关掉两层耦合）。
  function lookupFbaPichia(carbon, fbaTemp, strain, form, level) {
    var T = (typeof window !== "undefined") ? window.FBA_TABLE : null;
    if (!T || !T.data_pichia || !T.data_pichia[carbon]) return null;
    var tg = T.grid_pichia.fba_temp, sg = T.grid_pichia.strain, lg = T.grid_pichia.level;
    var fi = T.forms.indexOf(form); if (fi < 0) fi = 0;
    function bracket(arr, v) {
      if (v <= arr[0]) return [0, 0, 0];
      if (v >= arr[arr.length - 1]) return [arr.length - 1, arr.length - 1, 1];
      for (var i = 0; i < arr.length - 1; i++) {
        if (v >= arr[i] && v <= arr[i + 1]) {
          var tt = (arr[i + 1] === arr[i]) ? 0 : (v - arr[i]) / (arr[i + 1] - arr[i]);
          return [i, i + 1, tt];
        }
      }
      return [0, 0, 0];
    }
    function cell(ti, si) {
      var arr = T.data_pichia[carbon][ti][si][fi];
      if (level <= lg[0]) return arr[0];
      if (level >= lg[lg.length - 1]) return arr[arr.length - 1];
      for (var k = 0; k < lg.length - 1; k++) {
        if (level >= lg[k] && level <= lg[k + 1]) {
          var tl = (lg[k + 1] === lg[k]) ? 0 : (level - lg[k]) / (lg[k + 1] - lg[k]);
          var a = arr[k], b = arr[k + 1];
          return { drop: a.drop + (b.drop - a.drop) * tl, yac: a.yac + (b.yac - a.yac) * tl };
        }
      }
      return arr[0];
    }
    var tb = bracket(tg, fbaTemp), sb = bracket(sg, strain);
    var ti0 = tb[0], ti1 = tb[1], wt = tb[2];
    var si0 = sb[0], si1 = sb[1], ws = sb[2];
    var c00 = cell(ti0, si0), c10 = cell(ti1, si0), c01 = cell(ti0, si1), c11 = cell(ti1, si1);
    return {
      drop: c00.drop * (1 - wt) * (1 - ws) + c10.drop * wt * (1 - ws) + c01.drop * (1 - wt) * ws + c11.drop * wt * ws,
      yac: c00.yac * (1 - wt) * (1 - ws) + c10.yac * wt * (1 - ws) + c01.yac * (1 - wt) * ws + c11.yac * wt * ws
    };
  }

  // ---- 高密表达培养基预设（与 Python fip.twins.process.media.MEDIA_PRESETS 逐位镜像）----
  // 每个预设：label / composition / cn_ratio / recipe(覆盖配方字段) / params(覆盖机理参数) / note
  var MEDIA = {
    ecoli: {
      complex_yeast: {
        label: "复合培养基（葡萄糖+酵母粉）",
        composition: "葡萄糖（碳源/补料）+ 酵母浸出物（有机氮·维生素·生长因子）+ 蛋白胨 + K₂HPO₄/KH₂PO₄ + MgSO₄ + 微量元素",
        cn_ratio: "~4.5（酵母粉提供约 40% 有机氮，碳主要来自葡萄糖）",
        recipe: { x0_g_l: 0.5 },
        params: { ki_acetate: 8.0 },
        c_source: "葡萄糖", n_source: "酵母浸出物 + 蛋白胨（有机氮）",
        c_source_by_phase: { batch: "葡萄糖", fedbatch: "葡萄糖", induction: "葡萄糖" },
        note: "酵母粉/蛋白胨提供氨基酸前体，富培养基支持更高接种密度且对溢流乙酸耐受更高（ki_acetate 4→8）；在胁迫工况（高 μset / 低 DO）下乙酸抑制更轻，模型检测时效价差异放大；高密度补料分批可达 50–120 gDCW/L。"
      },
      defined: {
        label: "化学限定培养基（CDM）",
        composition: "葡萄糖 + (NH₄)₂SO₄（无机氮）+ K₂HPO₄ + MgSO₄ + 柠檬酸 + 微量元素 + 维生素",
        cn_ratio: "~5.0（碳完全来自葡萄糖）",
        recipe: { s0_g_l: 20.0, feed_s_g_l: 500.0, x0_g_l: 0.3 },
        params: { ki_acetate: 4.0, yxs: 0.45 },
        c_source: "葡萄糖", n_source: "(NH₄)₂SO₄（无机氮）",
        c_source_by_phase: { batch: "葡萄糖", fedbatch: "葡萄糖", induction: "葡萄糖" },
        note: "组分完全确定、可重复、计量明确；无机氮下葡萄糖溢流更易生乙酸；常用于工艺开发与放大研究（= 当前模型默认基线）。"
      }
    },
    pichia: {
      bsm_methanol: {
        label: "BSM 基础盐培养基 + 甘油/甲醇诱导",
        composition: "BSM（H₃PO₄ / K₂SO₄ / MgSO₄ / CaSO₄ / (NH₄)₂SO₄）+ PTM1 痕量盐 + 生物素；批式甘油 → 甘油流加 → 甲醇诱导",
        cn_ratio: "甲醇 C:N≈3.4（完全氧化）；甘油期 C:N≈3.0",
        recipe: {},
        params: {},
        c_source: "甘油 → 甲醇（诱导）", n_source: "(NH₄)₂SO₄（BSM 无机氮）",
        c_source_by_phase: { batch: "甘油", fedbatch: "甘油", induction: "甲醇" },
        note: "Invitrogen 标准 BSM；甘油批式长到 OD≈150，甘油流加维持 μ≈0.08，碳饥饿后甲醇诱导 AOX1；甲醇完全氧化需氧高（RQ<1）。= 当前模型默认基线。"
      }
    }
  };
  var DEFAULT_MEDIUM = { ecoli: "defined", pichia: "bsm_methanol" };
  function hostKeyOf(host) { var h = (host || "ecoli").toLowerCase(); return (h.indexOf("pichia") >= 0 || h.indexOf("pastoris") >= 0 || h.indexOf("yeast") >= 0) ? "pichia" : "ecoli"; }
  function mediaParams(host, medium) {
    var grp = MEDIA[hostKeyOf(host)]; if (!grp || !medium) return {};
    var m = grp[medium]; return (m && m.params) ? m.params : {};
  }
  function mediaRecipe(host, medium) {
    var grp = MEDIA[hostKeyOf(host)]; if (!grp || !medium) return {};
    var m = grp[medium]; return (m && m.recipe) ? m.recipe : {};
  }

  var NP_MAP = { rushton: 5.0, pitched_blade: 1.7, marine: 1.5, hydrofoil: 1.6 };
  function powerNumber(impeller) {
    var t = (impeller || "").toLowerCase();
    for (var k in NP_MAP) { if (t.indexOf(k) >= 0) return NP_MAP[k]; }
    return 1.7;
  }
  function powerPerVolumeWm3(rpm, di, impeller, volL, nImp, rho) {
    var v = volL / 1000.0;
    if (!(v > 0)) return 0.0;
    var nRps = rpm / 60.0;
    var np = powerNumber(impeller);
    var p = nImp * np * rho * Math.pow(nRps, 3) * Math.pow(di, 5);
    return p / v;
  }
  function klaH(pv, vvm) {
    if (vvm <= 0 || pv <= 0) return 0.0;
    return KLA_K * Math.pow(pv, 0.4) * Math.sqrt(vvm);
  }

  function clamp(v, lo, hi) { return v < lo ? lo : (v > hi ? hi : v); }

  // -------------------------------------------------------------------------
  // 阶段 H · 表达负担耦合：产物分叉（可溶/包涵体）+ 生长负担 + 可溶毒性
  // -------------------------------------------------------------------------
  // 表达形式预设（宿主无关，描述产物定位带来的固有负担）。与 Python
  // fip/twins/expression/model.py.EXPRESSION_FORMS 逐位一致。
  //   ib_frac      : 包涵体占产物比例（0..1）
  //   burden_sol  : 可溶产物合成对生长的负担系数 λsol (gX/gP)
  //   burden_ib   : 包涵体合成对生长的负担系数 λib (gX/gP)
  //   tox_k       : 可溶（toxic）产物的额外维持/死亡能耗系数 1/(g/L·h)（IB 因隔离不贡献毒性）
  var EXPRESSION_FORMS = {
    intracellular_soluble: { ib_frac: 0.10, burden_sol: 0.80, burden_ib: 0.20, tox_k: 0.004, label: "胞内可溶" },
    inclusion_body:        { ib_frac: 0.85, burden_sol: 0.20, burden_ib: 0.15, tox_k: 0.000, label: "包涵体" },
    periplasmic:           { ib_frac: 0.15, burden_sol: 0.40, burden_ib: 0.15, tox_k: 0.001, label: "周质空间" },
    secreted:              { ib_frac: 0.05, burden_sol: 0.25, burden_ib: 0.10, tox_k: 0.000, label: "分泌表达" }
  };
  // 默认（未选形式、未调滑块）：负担全 0 ⇒ 与改造前行为逐位一致（回归安全）
  var DEFAULT_IB_FRAC = 0.3, DEFAULT_BURDEN_SOL = 0.0, DEFAULT_BURDEN_IB = 0.0, DEFAULT_TOX_K = 0.0;

  // 解析有效负担参数：给定表达形式则采用其预设（表达形式=耦合主开关）；否则读滑块/字段（默认 0）。
  function resolveBurden(rec) {
    rec = rec || {};
    var form = rec.expression_form;
    var preset = (form && EXPRESSION_FORMS[form]) ? EXPRESSION_FORMS[form] : null;
    if (preset) {
      return { ib: preset.ib_frac, bs: preset.burden_sol, bi: preset.burden_ib, tk: preset.tox_k };
    }
    return {
      ib: num(rec.ib_frac, DEFAULT_IB_FRAC),
      bs: num(rec.burden_sol, DEFAULT_BURDEN_SOL),
      bi: num(rec.burden_ib, DEFAULT_BURDEN_IB),
      tk: num(rec.tox_k, DEFAULT_TOX_K)
    };
  }

  // -------------------------------------------------------------------------
  // E. coli 高密度补料分批（24h）：Monod + 乙酸抑制 + 溢流 + T7 诱导表达
  // -------------------------------------------------------------------------
  function simulateEcoli(rec) {
    rec = rec || {};
    // 培养基预设配方覆盖（仅当用户未显式给定的字段；滑块/显式字段优先）
    var _mr = mediaRecipe("ecoli", rec.medium);
    for (var _mk in _mr) { if (rec[_mk] === undefined) rec[_mk] = _mr[_mk]; }
    var p = {
      mu_max: 0.55, ks: 0.05, yxs: 0.45, q_s_max: 1.2, q_p_max: 0.012,
      ki_acetate: 4.0, qo2_growth: 12.0, mo2: 1.5, c_star_mm: C_STAR_AIR,
      ki_do: 0.0002  // 与 Python mechanistic.py 对齐（≈0.08% 空气饱和度，氧限制半饱和常数）
    };
    // 培养基预设机理参数覆盖（复合培养基 ki_acetate↑ 缓解溢流乙酸）
    var _mp = mediaParams("ecoli", rec.medium);
    for (var _pk in _mp) { p[_pk] = _mp[_pk]; }
    // 阶段 FBA 耦合：用预计算查表灌入内禀参数（关掉两层各算各的）；默认不耦合保持回归安全
    var FBA = (rec.fba_coupled && typeof window !== "undefined" && window.FBA_TABLE)
      ? lookupFba(rec.carbon_source || "glucose", num(rec.fba_temp, 37.0), num(rec.strain_factor, 1.0),
                  rec.expression_form || "none", num(rec.expression_level, 0.15))
      : null;
    var FBA_DROP = 1.0, FBA_YAC = 0.7;
    if (FBA) {
      p.mu_max = MU_REF_37;
      p.q_s_max = effectiveUptake(rec.carbon_source || "glucose", num(rec.fba_temp, 37.0), num(rec.strain_factor, 1.0))
                 * MM_CARBON[rec.carbon_source || "glucose"] / 1000.0;
      FBA_DROP = (FBA.drop != null) ? FBA.drop : 1.0;
      FBA_YAC = (FBA.yac != null) ? FBA.yac : 0.7;
      if (typeof console !== "undefined") console.log("[FBA-coupled] carbon=" + (rec.carbon_source || "glucose")
        + " T=" + num(rec.fba_temp, 37.0) + " strain=" + num(rec.strain_factor, 1.0)
        + ": mu_max=" + p.mu_max.toFixed(3) + ", q_s_max=" + p.q_s_max.toFixed(3) + " gS/gX/h, drop_frac="
        + FBA_DROP.toFixed(3) + ", y_ac=" + FBA_YAC.toFixed(4) + " g/g");
    }
    var rc = {
      batch_volume_l: num(rec.batch_volume_l, 10), t_end_h: num(rec.t_end_h, 24),
      s0_g_l: num(rec.s0_g_l, 20), x0_g_l: num(rec.x0_g_l, 0.3),
      induction_h: num(rec.induction_h, 16), temp_pre: num(rec.temp_pre, 37),
      temp_post: num(rec.temp_post, 30), feed_start_h: num(rec.feed_start_h, 7),
      mu_set: num(rec.mu_set, 0.11), feed_s_g_l: num(rec.feed_s_g_l, 500),
      feed_s_g_l_post: num(rec.feed_s_g_l_post, 0),
      ph_set: num(rec.ph_set, 7.0), rpm: num(rec.rpm, 800),
      airflow_lmin: num(rec.airflow_lmin, 10), press_bar: num(rec.press_bar, 1.0),
      kla_scale: num(rec.kla_scale, 1.0), di_m: num(rec.di_m, 0.08),
      n_imp: num(rec.n_imp, 2), impeller_type: rec.impeller_type || "rushton",
      rho: num(rec.rho, 1030),
      ph_curve: rec.ph_curve, temp_curve: rec.temp_curve,
      ib_frac: num(rec.ib_frac, DEFAULT_IB_FRAC), burden_sol: num(rec.burden_sol, DEFAULT_BURDEN_SOL),
      burden_ib: num(rec.burden_ib, DEFAULT_BURDEN_IB), tox_k: num(rec.tox_k, DEFAULT_TOX_K),
      expression_form: rec.expression_form || null,
      // 阶段 FBA 耦合字段（默认不耦合：手调默认值，不继承 FBA μ_het）
      carbon_source: rec.carbon_source || "glucose",
      strain_factor: num(rec.strain_factor, 1.0),
      fba_coupled: !!rec.fba_coupled,
      fba_temp: num(rec.fba_temp, 37.0),
      expression_level: num(rec.expression_level, 0.15)
    };
    // 阶段 H · 表达负担耦合：给定表达形式用其预设（主开关），否则用滑块/字段（默认 0 ⇒ 回归安全）
    var _eb = resolveBurden(rec);
    var IB_FRAC = _eb.ib, BURDEN_SOL = _eb.bs, BURDEN_IB = _eb.bi, TOX_K = _eb.tk;
    var dt = num(rec.dt_h, 0.1), noise = num(rec.noise, 0), n = Math.floor(rc.t_end_h / dt) + 1;
    var rng = mulberry(Math.floor((rec.seed || 42) * 1) || 42);

    var X = rc.x0_g_l, S = rc.s0_g_l, P = 0.0, A = 0.0, Psol = 0.0, Pib = 0.0, vol = rc.batch_volume_l;
    var arrs = emptyArrays(n);
    var vvm = rc.airflow_lmin / rc.batch_volume_l;

    for (var i = 0; i < n; i++) {
      var t = i * dt;
      var induced = t >= rc.induction_h;
      var temp = curveAt(rc, "temp_curve", t, induced ? rc.temp_post : rc.temp_pre);
      var ph = curveAt(rc, "ph_curve", t, rc.ph_set);
      // 对称惩罚：高温（热应激）与低温（代谢减慢）两侧都应降低 μ；37 ℃ ±1 ℃ 视为最优区
      var fTemp = FBA ? growthTemperatureFactor(temp)
                      : (Math.abs(temp - 37.0) < 1.0 ? 1.0 : Math.exp(-0.06 * Math.abs(37.0 - temp)));
      var fPh = Math.abs(ph - rc.ph_set) < 0.3 ? 1.0 : Math.exp(-0.25 * Math.abs(ph - rc.ph_set));

      // 诱导期可切换补料液浓度（实验设计 Induction 阶段碳源浓度）：质量流率由 mu_set 决定，
      // 浓度只改变体积流率 → 影响稀释与最终装液量（浓度越高、稀释越小）
      var feedC = (induced && rc.feed_s_g_l_post > 0) ? rc.feed_s_g_l_post : rc.feed_s_g_l;
      var F_lh = 0.0;
      if (t >= rc.feed_start_h && X >= 0.5) F_lh = (rc.mu_set * X * vol / p.yxs) / feedC;
      var feed_g = F_lh * feedC / vol;

      var inhibition = 1.0 / (1.0 + A / p.ki_acetate);
      var monod = S > 0 ? S / (p.ks + S) : 0.0;
      var muPotRaw = p.mu_max * fTemp * fPh * monod * inhibition;

      // 碳平衡（先用未限 DO 的潜在 μ 估 realized OUR，供 DO 级联）
      var avail = S / dt + feed_g;
      var uptakePot = p.q_s_max * fTemp * monod * X;
      var uptake = Math.min(uptakePot, avail);
      var growthRaw = Math.min(muPotRaw * X, uptake * p.yxs);

      // DO 限制（X 生长耦合 DO）：用 realized OUR 估 CL，再 Monod 限生长
      var ourTent = X * (p.qo2_growth * (growthRaw / X) + p.mo2);
      var rpm = rc.rpm, air = rc.airflow_lmin, o2frac = 0.21, doSat = 100.0;
      var klaLast = 0.0, cStarLast = 0.0, cLTent = 0.0;
      for (var _ = 0; _ < 16; _++) {
        var pv = powerPerVolumeWm3(rpm, rc.di_m, rc.impeller_type, vol, rc.n_imp, rc.rho);
        klaLast = rc.kla_scale * klaH(pv, vvm);
        cStarLast = p.c_star_mm * (o2frac / 0.21) * rc.press_bar;
        cLTent = Math.max(cStarLast - ourTent / klaLast, 0.0);
        doSat = clamp(cLTent / cStarLast * 100.0, 0.0, 100.0);
        if (doSat >= 25.0 || rpm >= rc.rpm * 1.6) break;
        rpm = Math.min(rpm * 1.15, rc.rpm * 1.6);
        o2frac = Math.min(o2frac + 0.12, 0.95);
      }
      var fDo = cLTent / (p.ki_do + cLTent);
      var muPot = muPotRaw * fDo;
      var growth = Math.min(muPot * X, uptake * p.yxs);

      var acetateGen = Math.max(0.0, uptake - growth / p.yxs) * (FBA ? FBA_YAC : 0.7);
      var acetateCons = 0.08 * A * X / (1.0 + S);
      var qP = p.q_p_max * (induced ? 1.0 : 0.0) * inhibition * fPh;

      // 阶段 H · 产物分叉（可溶/包涵体）+ 表达负担对生长的耦合
      var burdenEff = BURDEN_SOL * (1.0 - IB_FRAC) + BURDEN_IB * IB_FRAC;
      var dPsol = (1.0 - IB_FRAC) * qP * X;
      var dPib = IB_FRAC * qP * X;
      // 阶段 H · 可溶毒性 = 额外维持/死亡能耗（与可溶毒性产物浓度成正比；IB 因隔离而不贡献）
      // 阶段 FBA 耦合：表达负担用 FBA μ_het 降幅 drop_frac 直接压生长（关掉内生 burden_eff·q_p）；
      // 可溶毒性 tox_k 两种模式都保留（FBA 不建模毒性死亡）。
      var dX = FBA ? (growth * FBA_DROP - 0.01 * X - TOX_K * Psol * X)
                   : (growth - 0.01 * X - burdenEff * qP * X - TOX_K * Psol * X);
      var dS = feed_g - uptake;
      var dA = acetateGen - acetateCons;
      var dP = dPsol + dPib;

      X = Math.max(X + dX * dt, 0.0);
      S = Math.max(S + dS * dt, 0.0);
      A = Math.max(A + dA * dt, 0.0);
      P = Math.max(P + dP * dt, 0.0);
      Psol = Math.max(Psol + dPsol * dt, 0.0);
      Pib = Math.max(Pib + dPib * dt, 0.0);
      vol += F_lh * dt;

      var muEff = X > 1e-9 ? growth / X : 0.0;
      var our = X * (p.qo2_growth * muEff + p.mo2);
      var cer = Math.max(our * 1.05 - 0.5 * acetateCons, 0.0);
      var rq = our > 1e-6 ? cer / our : 1.0;

      // DO 最终值：用实生长 OUR 与已确定执行器重算（≥ 级联 tentative）
      var cLFinal = klaLast > 0 ? Math.max(cStarLast - our / klaLast, 0.0) : 0.0;
      doSat = cStarLast > 0 ? clamp(cLFinal / cStarLast * 100.0, 0.0, 100.0) : 0.0;

      arrs.t[i] = t; arrs.X[i] = X; arrs.S[i] = S; arrs.P[i] = P; arrs.A[i] = A;
      arrs.Psol[i] = Psol; arrs.Pib[i] = Pib;
      arrs.DO[i] = clamp(nz(doSat, noise, rng), 0, 100);
      arrs.OUR[i] = nz(our, noise, rng);
      arrs.CER[i] = nz(cer, noise * 0.5, rng);
      arrs.RQ[i] = nz(rq, noise * 0.5, rng);
      arrs.FEED[i] = F_lh; arrs.RPM[i] = rpm; arrs.AIR[i] = air;
      arrs.TEMP[i] = nz(temp, noise, rng); arrs.PH[i] = nz(ph, noise, rng);
      arrs.INDUCED[i] = induced ? 1.0 : 0.0;
    }
    return finalize(arrs, "ecoli", rc, dt);
  }

  // -------------------------------------------------------------------------
  // Pichia pastoris 动态发酵（默认 96h），支持两种表达模式：
  //   · constitutive（组成型）：碳源全程为甘油（甘油批 → 甘油补料），蛋白全程表达
  //   · inducible（诱导型）：甘油批 → 甘油补料 → 碳饥饿 → 甲醇诱导；
  //       诱导型补料模式 meoh_feed_mode = methanol（纯甲醇）或 mixed（甲醇+甘油混合共利用）
  // 阶段 FBA 耦合：读 window.FBA_TABLE.data_pichia[carbon]（碳源 glycerol/methanol），
  // 表缺失或 fba_coupled=false 时退化为手调默认（关掉两层耦合，回归安全）。
  // -------------------------------------------------------------------------
  function simulatePichia(rec) {
    rec = rec || {};
    // 培养基预设配方覆盖（仅当用户未显式给定的字段）
    var _mr = mediaRecipe("pichia", rec.medium);
    for (var _mk in _mr) { if (rec[_mk] === undefined) rec[_mk] = _mr[_mk]; }

    var exprMode = rec.expression_mode || "inducible";      // constitutive / inducible
    var meohMode = rec.meoh_feed_mode || "methanol";        // methanol / mixed
    var p = {
      mu_max_gly: 0.28, mu_max_meoh: 0.14, ks_gly: 0.08, ks_meoh: 0.25,
      yxs_gly: 0.50, yxs_meoh: 0.42, q_p_max: 0.004, q_meoh_maint: 0.03,
      ki_meoh: 4.0, qo2_growth: 10.0, mo2: 1.2, o2_per_g_meoh: 46.9, c_star_mm: C_STAR_AIR
    };
    // 培养基预设机理参数覆盖（BSM 预设 params 为空，行为不变）
    var _mp = mediaParams("pichia", rec.medium);
    for (var _pk in _mp) { p[_pk] = _mp[_pk]; }

    // 阶段 FBA 耦合（毕赤酵母）：表达发生在哪个碳源，就用哪个碳源的 FBA 校准
    // （组成型在甘油上表达 → glycerol；诱导型在甲醇上表达 → methanol）
    var fbaCarbon = (exprMode === "constitutive") ? "glycerol" : "methanol";
    var FBAp = (rec.fba_coupled && typeof window !== "undefined" && window.FBA_TABLE && window.FBA_TABLE.data_pichia)
      ? lookupFbaPichia(fbaCarbon, num(rec.fba_temp, 28.0), num(rec.strain_factor, 1.0),
                        rec.expression_form || "none", num(rec.expression_level, 0.15))
      : null;
    var FBA_DROP = 1.0, FBA_YAC = 0.0;
    if (FBAp) {
      // FBA 内禀上限（替代手调 μ_max），生物质得率也用 FBA 值
      p.mu_max_gly = MU_REF_P_GLY; p.mu_max_meoh = MU_REF_P_MEOM;
      p.yxs_gly = YXS_P.glycerol; p.yxs_meoh = YXS_P.methanol;
      FBA_DROP = (FBAp.drop != null) ? FBAp.drop : 1.0;
      FBA_YAC = (FBAp.yac != null) ? FBAp.yac : 0.0;        // 甲醇 DO 受限溢出 yield（甲醛/g 甲醇）
      if (typeof console !== "undefined") console.log("[FBA-coupled:Pichia] carbon=" + fbaCarbon
        + " mode=" + exprMode + " T=" + num(rec.fba_temp, 28.0) + " strain=" + num(rec.strain_factor, 1.0)
        + ": mu_max_gly=" + p.mu_max_gly.toFixed(3) + ", mu_max_meoh=" + p.mu_max_meoh.toFixed(3)
        + ", drop_frac=" + FBA_DROP.toFixed(3) + ", overflow_yac=" + FBA_YAC.toFixed(4));
    }

    var rc = {
      batch_volume_l: num(rec.batch_volume_l, 10), t_end_h: num(rec.t_end_h, 96),
      gly0_g_l: num(rec.s0_g_l, 20), x0_g_l: num(rec.x0_g_l, 0.5),
      gly_feed_start_h: num(rec.gly_feed_start_h, 18), gly_feed_end_h: num(rec.gly_feed_end_h, 32),
      mu_set_gly: num(rec.mu_set, 0.08), gly_feed_g_l: num(rec.feed_s_g_l, 500),
      meoh_start_h: num(rec.meoh_start_h, 34), meoh_feed_rate_lh: num(rec.meoh_feed_rate_lh, 0.012),
      meoh_density: 0.792, temp_c: num(rec.temp_pre, 28), ph_set: num(rec.ph_set, 5.5),
      rpm: num(rec.rpm, 900), airflow_lmin: num(rec.airflow_lmin, 10), press_bar: num(rec.press_bar, 1.0),
      kla_scale: num(rec.kla_scale, 1.0), di_m: num(rec.di_m, 0.08),
      n_imp: num(rec.n_imp, 2), impeller_type: rec.impeller_type || "rushton", rho: num(rec.rho, 1030),
      ph_curve: rec.ph_curve, temp_curve: rec.temp_curve,
      ib_frac: num(rec.ib_frac, DEFAULT_IB_FRAC), burden_sol: num(rec.burden_sol, DEFAULT_BURDEN_SOL),
      burden_ib: num(rec.burden_ib, DEFAULT_BURDEN_IB), tox_k: num(rec.tox_k, DEFAULT_TOX_K),
      expression_form: rec.expression_form || null,
      // 阶段 FBA 耦合字段
      fba_coupled: !!rec.fba_coupled, fba_temp: num(rec.fba_temp, 28.0),
      strain_factor: num(rec.strain_factor, 1.0), expression_level: num(rec.expression_level, 0.15),
      // 表达模式 / 诱导型补料模式
      expression_mode: exprMode, meoh_feed_mode: meohMode,
      mixed_gly_feed_g_l: num(rec.mixed_gly_feed_g_l, 200.0),
      mixed_gly_start_h: num(rec.mixed_gly_start_h, 0),  // 0 ⇒ 沿用 meoh_start_h
      mixed_gly_end_h: num(rec.mixed_gly_end_h, 0)        // 0 ⇒ 沿用 t_end_h
    };
    // 阶段 H · 表达负担耦合：给定表达形式用其预设（主开关），否则用滑块/字段（默认 0 ⇒ 回归安全）
    var _eb = resolveBurden(rec);
    var IB_FRAC = _eb.ib, BURDEN_SOL = _eb.bs, BURDEN_IB = _eb.bi, TOX_K = _eb.tk;
    var dt = num(rec.dt_h, 0.25), noise = num(rec.noise, 0), n = Math.floor(rc.t_end_h / dt) + 1;
    var rng = mulberry(Math.floor((rec.seed || 7) * 1) || 7);

    var X = rc.x0_g_l, S = rc.gly0_g_l, M = 0.0, P = 0.0, Psol = 0.0, Pib = 0.0, vol = rc.batch_volume_l, aox = 0.0;
    var arrs = emptyArrays(n);
    // 组成型无甲醇诱导阶段
    var hasMeohPhase = (exprMode === "inducible");
    var vvm = rc.airflow_lmin / rc.batch_volume_l;

    for (var i = 0; i < n; i++) {
      var t = i * dt;
      var onMeoh = hasMeohPhase && (t >= rc.meoh_start_h);
      var temp = curveAt(rc, "temp_curve", t, rc.temp_c);
      var ph = curveAt(rc, "ph_curve", t, rc.ph_set);
      // Pichia 专属 cardinal 温度因子（FBA 耦合与否都用同一物理基准，替代原 exp 近似）
      var fTemp = growthTemperatureFactorPichia(temp);
      var fPh = Math.abs(ph - rc.ph_set) < 0.3 ? 1.0 : Math.exp(-0.30 * Math.abs(ph - rc.ph_set));
      // AOX（醇氧化酶）诱导：仅诱导型甲醇相累积
      if (onMeoh && M > 0.1) aox = Math.min(aox + 0.12 * dt, 1.0);
      else if (!onMeoh) aox = Math.max(aox - 0.05 * dt, 0.0);

      // 甘油补料：组成型全程；诱导型仅 meoh 前
      var glyFeedActive = (t >= rc.gly_feed_start_h && t < rc.gly_feed_end_h && X > 0.5);
      var fGly = 0.0;
      if (glyFeedActive)
        fGly = (rc.mu_set_gly * X * vol / p.yxs_gly) / rc.gly_feed_g_l;
      var feedGlyG = fGly * rc.gly_feed_g_l / vol;

      // 甲醇补料（仅诱导型）
      var fMeoh = onMeoh ? rc.meoh_feed_rate_lh : 0.0;
      var feedMeohG = fMeoh * rc.meoh_density * 1000.0 / vol;

      // 诱导型 · 甲醇+甘油混合补料（共利用）：诱导期内额外甘油 co-feed
      if (exprMode === "inducible" && meohMode === "mixed" && onMeoh) {
        var mgStart = rc.mixed_gly_start_h > 0 ? rc.mixed_gly_start_h : rc.meoh_start_h;
        var mgEnd = rc.mixed_gly_end_h > 0 ? rc.mixed_gly_end_h : rc.t_end_h;
        if (t >= mgStart && t < mgEnd && X > 0.5) {
          var fGlyMix = (rc.mu_set_gly * 0.6 * X * vol / p.yxs_gly) / rc.mixed_gly_feed_g_l;
          feedGlyG += fGlyMix * rc.mixed_gly_feed_g_l / vol;
        }
      }

      var monodG = S > 0 ? S / (p.ks_gly + S) : 0.0;
      var availG = S / dt + feedGlyG;
      var growthG = Math.min(p.mu_max_gly * fTemp * fPh * monodG * X, availG * p.yxs_gly);
      var uptakeG = growthG / p.yxs_gly + 0.015 * X;

      var meohInhib = 1.0 / (1.0 + M / p.ki_meoh);
      var monodM = M > 0 ? M / (p.ks_meoh + M) : 0.0;
      var availM = M / dt + feedMeohG;
      var growthM = Math.min(p.mu_max_meoh * fTemp * fPh * aox * monodM * meohInhib * X, availM * p.yxs_meoh);
      var maintMeoh = p.q_meoh_maint * aox * X * (onMeoh ? 1.0 : 0.0);
      var consMeoh = Math.min(growthM / p.yxs_meoh + maintMeoh, availM);

      var growth = growthG + growthM;
      // 表达活跃？组成型全程；诱导型仅在甲醇相
      var exprActive = (exprMode === "constitutive") ? 1.0 : (onMeoh ? 1.0 : 0.0);
      var qP = (exprMode === "constitutive")
        ? p.q_p_max * fTemp * fPh
        : p.q_p_max * fTemp * fPh * aox * (onMeoh ? 1.0 : 0.0);

      // 阶段 H · 产物分叉（可溶/包涵体）+ 表达负担 + 可溶毒性
      var burdenEff = BURDEN_SOL * (1.0 - IB_FRAC) + BURDEN_IB * IB_FRAC;
      var dPsol = (1.0 - IB_FRAC) * qP * X;
      var dPib = IB_FRAC * qP * X;
      // 阶段 FBA 耦合：表达负担用 FBA μ_het 降幅 drop_frac 直接压表达相生长（关掉内生 burden_eff·q_p）；
      // 甲醇 DO 受限溢出（甲醇→甲醛毒性）用 FBA overflow yield 额外压生长。
      var growthEff = (FBAp && exprActive) ? growth * FBA_DROP : growth;
      var meohTox = FBAp ? FBA_YAC * Math.max(0.0, M - p.ki_meoh) : 0.0;  // g 甲醛/g 甲醇 × 超量甲醇
      var dX = growthEff - 0.008 * X
             - (FBAp ? 0.0 : burdenEff * qP * X)   // 非 FBA 耦合时用内生负担系数
             - TOX_K * Psol * X - meohTox * X;
      var dS = feedGlyG - uptakeG;
      var dM = feedMeohG - consMeoh;
      var dP = dPsol + dPib;

      X = Math.max(X + dX * dt, 0.0);
      S = Math.max(S + dS * dt, 0.0);
      M = Math.max(M + dM * dt, 0.0);
      P = Math.max(P + dP * dt, 0.0);
      Psol = Math.max(Psol + dPsol * dt, 0.0);
      Pib = Math.max(Pib + dPib * dt, 0.0);
      vol += (fGly + fMeoh) * dt;

      var muEff = X > 1e-9 ? growth / X : 0.0;
      var our = X * (p.qo2_growth * muEff + p.mo2);
      our += (consMeoh - growthM / p.yxs_meoh) * p.o2_per_g_meoh;
      var rq = onMeoh ? 0.75 - 0.05 * aox : 1.0;
      var cer = Math.max(our * rq, 0.0);

      var rpm = rc.rpm, air = rc.airflow_lmin, o2frac = 0.21, doSat = 100.0;
      for (var _ = 0; _ < 3; _++) {
        var pv = powerPerVolumeWm3(rpm, rc.di_m, rc.impeller_type, vol, rc.n_imp, rc.rho);
        var kla = rc.kla_scale * klaH(pv, vvm);
        var cStar = p.c_star_mm * (o2frac / 0.21) * rc.press_bar;
        var cL = Math.max(cStar - our / kla, 0.0);
        doSat = clamp(cL / cStar * 100.0, 0.0, 100.0);
        if (doSat >= 25.0 || rpm >= rc.rpm * 1.6) break;
        rpm = Math.min(rpm * 1.2, rc.rpm * 1.6);
        o2frac = Math.min(o2frac + 0.15, 0.95);
      }

      // 操作阶段标注（供渲染/诊断）
      var phase;
      if (exprMode === "constitutive") {
        phase = (t < rc.gly_feed_start_h) ? "gly_bat" : "gly_feed";
      } else if (t < rc.meoh_start_h) {
        phase = (t < rc.gly_feed_end_h) ? "gly_feed" : "starv";
      } else {
        phase = (meohMode === "mixed") ? "mix_feed" : "meoh_induce";
      }

      arrs.t[i] = t; arrs.X[i] = X; arrs.S[i] = S; arrs.P[i] = P; arrs.M[i] = M;
      arrs.Psol[i] = Psol; arrs.Pib[i] = Pib; arrs.MODE[i] = phase;
      arrs.DO[i] = clamp(nz(doSat, noise, rng), 0, 100);
      arrs.OUR[i] = nz(our, noise, rng);
      arrs.CER[i] = nz(cer, noise * 0.5, rng);
      arrs.RQ[i] = nz(rq, noise * 0.5, rng);
      arrs.FEED[i] = fMeoh; arrs.RPM[i] = rpm; arrs.AIR[i] = air;
      arrs.TEMP[i] = nz(temp, noise, rng); arrs.PH[i] = nz(ph, noise, rng);
      arrs.INDUCED[i] = (onMeoh ? 1.0 : 0.0); arrs.AOX[i] = aox;
    }
    return finalize(arrs, "pichia", rc, dt);
  }

  // -------------------------------------------------------------------------
  // 公用辅助
  // -------------------------------------------------------------------------
  function num(v, d) { return (v == null || isNaN(v)) ? d : v; }
  // 分段曲线插值（与 Python _interp_curve 逐位一致）：rec[key] = [[t_h, value], ...]（按 t 升序）
  function curveAt(rec, key, t, def) {
    var c = rec[key];
    if (!c || !c.length) return def;
    if (t <= c[0][0]) return c[0][1];
    for (var j = 0; j < c.length - 1; j++) {
      var t0 = c[j][0], v0 = c[j][1], t1 = c[j + 1][0], v1 = c[j + 1][1];
      if (t0 <= t && t <= t1) {
        if (t1 === t0) return v1;
        return v0 + (v1 - v0) * (t - t0) / (t1 - t0);
      }
    }
    return c[c.length - 1][1];
  }
  function emptyArrays(n) {
    var keys = ["t", "X", "S", "P", "M", "A", "DO", "OUR", "CER", "RQ",
      "FEED", "RPM", "AIR", "TEMP", "PH", "INDUCED", "AOX", "Psol", "Pib", "MODE"];
    var o = {}; keys.forEach(function (k) { o[k] = new Array(n); }); return o;
  }
  function nz(v, rel, rng) {
    if (!(rel > 0)) return v;
    return v * (1 + rng() * 0 - 0) * (1 + gauss(rng) * rel);
  }
  function gauss(rng) {
    var u = 0, v = 0;
    while (u === 0) u = rng();
    while (v === 0) v = rng();
    return Math.sqrt(-2.0 * Math.log(u)) * Math.cos(2.0 * Math.PI * v);
  }
  // 轻量确定性 PRNG（保证每次"相同配方 → 相同动画"，便于复现与对照）
  function mulberry(seed) {
    var a = seed >>> 0;
    return function () {
      a |= 0; a = (a + 0x6D2B79F5) | 0;
      var t = Math.imul(a ^ (a >>> 15), 1 | a);
      t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
      return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
    };
  }
  function finalize(arrs, host, rc, dt) {
    var t = arrs.t, P = arrs.P, X = arrs.X, A = arrs.A, DO = arrs.DO, OUR = arrs.OUR;
    var n = t.length;
    var acetatePeak = 0, doMin = 100, ourPeak = 0;
    for (var i = 0; i < n; i++) {
      if (A[i] > acetatePeak) acetatePeak = A[i];
      if (DO[i] < doMin) doMin = DO[i];
      if (OUR[i] > ourPeak) ourPeak = OUR[i];
    }
    arrs.host = host; arrs.dt = dt; arrs.recipe = rc; arrs.n = n;
    arrs.summary = {
      titer: P[n - 1], dcw: X[n - 1], acetatePeak: acetatePeak,
      doMin: doMin, ourPeak: ourPeak, endH: t[n - 1]
    };
    return arrs;
  }

  function simulate(rec) {
    rec = rec || {};
    var host = (rec.host || "ecoli").toLowerCase();
    if (host.indexOf("pichia") >= 0 || host.indexOf("pastoris") >= 0) return simulatePichia(rec);
    return simulateEcoli(rec);
  }

  // 默认配方（与 Python 生成器一致），供页面初始化与"重置"复用；medium 命中预设时取预设配方
  function defaultRecipe(host, medium) {
    host = (host || "ecoli").toLowerCase();
    var hk = (host.indexOf("pichia") >= 0 || host.indexOf("pastoris") >= 0) ? "pichia" : "ecoli";
    if (!medium) medium = DEFAULT_MEDIUM[hk];
    var base;
    if (hk === "pichia") {
      base = { batch_volume_l: 10, t_end_h: 96, s0_g_l: 20, x0_g_l: 0.5,
        gly_feed_start_h: 18, gly_feed_end_h: 32, mu_set: 0.08, feed_s_g_l: 500,
        meoh_start_h: 34, meoh_feed_rate_lh: 0.012, temp_pre: 28, ph_set: 5.5,
        rpm: 900, airflow_lmin: 10, press_bar: 1.0, kla_scale: 1.0, di_m: 0.08,
        n_imp: 2, impeller_type: "rushton", rho: 1030, dt_h: 0.25, noise: 0, seed: 7,
        ib_frac: 0.3, burden_sol: 0.0, burden_ib: 0.0, tox_k: 0.0, expression_form: null,
        // 表达模式 / 诱导型补料模式（默认诱导型·纯甲醇，与原行为一致）
        expression_mode: "inducible", meoh_feed_mode: "methanol", mixed_gly_feed_g_l: 200.0,
        // 阶段 FBA 耦合（默认关闭：手调默认值，不继承 FBA μ_het；表缺失时同样退化为默认）
        fba_coupled: false, fba_temp: 28.0, strain_factor: 1.0, expression_level: 0.15 };
    } else {
      base = { batch_volume_l: 10, t_end_h: 24, s0_g_l: 20, x0_g_l: 0.3,
        induction_h: 16, temp_pre: 37, temp_post: 30, feed_start_h: 7, mu_set: 0.11,
        feed_s_g_l: 500, ph_set: 7.0, rpm: 800, airflow_lmin: 10, press_bar: 1.0,
        kla_scale: 1.0, di_m: 0.08, n_imp: 2, impeller_type: "rushton", rho: 1030,
        dt_h: 0.1, noise: 0, seed: 42,
        ib_frac: 0.3, burden_sol: 0.0, burden_ib: 0.0, tox_k: 0.0, expression_form: null };
    }
    // 培养基预设配方覆盖（预设优先于硬编码默认）
    var mr = mediaRecipe(hk, medium);
    for (var k in mr) base[k] = mr[k];
    base.host = hk;
    base.medium = medium;
    return base;
  }

  // 参数范围（驱动滑块 UI；min/max/step/unit/cn）
  var PARAM_RANGES = {
    ecoli: [
      { key: "rpm", label: "搅拌转速", unit: "rpm", min: 300, max: 1400, step: 10, def: 800 },
      { key: "airflow_lmin", label: "空气通量", unit: "L/min", min: 2, max: 30, step: 0.5, def: 10 },
      { key: "induction_h", label: "诱导时刻", unit: "h", min: 4, max: 22, step: 0.5, def: 16 },
      { key: "mu_set", label: "补料比生长速率 μset", unit: "1/h", min: 0.04, max: 0.20, step: 0.005, def: 0.11 },
      { key: "temp_post", label: "诱导后温度", unit: "℃", min: 16, max: 37, step: 1, def: 30 },
      { key: "s0_g_l", label: "初始葡萄糖", unit: "g/L", min: 5, max: 40, step: 1, def: 20 },
      { key: "x0_g_l", label: "接种菌浓 X0", unit: "g/L", min: 0.05, max: 2, step: 0.05, def: 0.3 },
      { key: "kla_scale", label: "传质折扣 kLa×", unit: "", min: 0.5, max: 1.5, step: 0.05, def: 1.0 },
      { key: "ib_frac", label: "包涵体比例", unit: "", min: 0, max: 1, step: 0.05, def: 0.3 },
      { key: "burden_sol", label: "可溶表达负担 λsol", unit: "", min: 0, max: 2, step: 0.1, def: 0.0 },
      { key: "burden_ib", label: "包涵体负担 λib", unit: "", min: 0, max: 2, step: 0.1, def: 0.0 },
      { key: "tox_k", label: "可溶毒性系数 kT", unit: "1/(g/L·h)", min: 0, max: 0.02, step: 0.001, def: 0.0 }
    ],
    pichia: [
      { key: "rpm", label: "搅拌转速", unit: "rpm", min: 300, max: 1400, step: 10, def: 900 },
      { key: "airflow_lmin", label: "空气通量", unit: "L/min", min: 2, max: 30, step: 0.5, def: 10 },
      { key: "meoh_start_h", label: "甲醇诱导时刻", unit: "h", min: 20, max: 60, step: 1, def: 34 },
      { key: "meoh_feed_rate_lh", label: "甲醇补料速率", unit: "L/h", min: 0.002, max: 0.04, step: 0.001, def: 0.012 },
      { key: "mu_set", label: "甘油补料 μset", unit: "1/h", min: 0.03, max: 0.15, step: 0.005, def: 0.08 },
      { key: "temp_pre", label: "培养温度", unit: "℃", min: 20, max: 32, step: 1, def: 28 },
      { key: "x0_g_l", label: "接种菌浓 X0", unit: "g/L", min: 0.05, max: 2, step: 0.05, def: 0.5 },
      { key: "kla_scale", label: "传质折扣 kLa×", unit: "", min: 0.5, max: 1.5, step: 0.05, def: 1.0 },
      { key: "ib_frac", label: "包涵体比例", unit: "", min: 0, max: 1, step: 0.05, def: 0.3 },
      { key: "burden_sol", label: "可溶表达负担 λsol", unit: "", min: 0, max: 2, step: 0.1, def: 0.0 },
      { key: "burden_ib", label: "包涵体负担 λib", unit: "", min: 0, max: 2, step: 0.1, def: 0.0 },
      { key: "tox_k", label: "可溶毒性系数 kT", unit: "1/(g/L·h)", min: 0, max: 0.02, step: 0.001, def: 0.0 }
    ]
  };

  window.FermentSim = {
    simulate: simulate, defaultRecipe: defaultRecipe, paramRanges: PARAM_RANGES,
    MEDIA: MEDIA, DEFAULT_MEDIUM: DEFAULT_MEDIUM, EXPRESSION_FORMS: EXPRESSION_FORMS,
    resolveBurden: resolveBurden,
    lookupFba: lookupFba, growthTemperatureFactor: growthTemperatureFactor,
    effectiveUptake: effectiveUptake, muCap: muCap,
    // 毕赤酵母专属 FBA 两层耦合
    lookupFbaPichia: lookupFbaPichia, growthTemperatureFactorPichia: growthTemperatureFactorPichia,
    effectiveUptakePichia: effectiveUptakePichia, muCapPichia: muCapPichia,
    KLA_K: KLA_K, _internals: { klaH: klaH, powerPerVolumeWm3: powerPerVolumeWm3 }
  };
})();
