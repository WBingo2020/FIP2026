/* ============================================================================
 * FIP 发酵智能平台 · 自包含交互引擎 (dash.js)
 * 纯原生 JS + SVG，零外部依赖。数据由 window.DASH_DATA 注入。
 * 计算逻辑移植自 fip.twins (Expression 7层 / Scale-up 几何相似放大 / Process 软测量)。
 * ==========================================================================*/
(function () {
  "use strict";

  // 宿主键归一化（与 ferment_sim.js 中 hostKeyOf 同源）
  function hostKeyOf(host) {
    var h = (host || "ecoli").toLowerCase();
    // 实验设计模块 host_type 取值为 yeast/ecoli；仿真引擎为 pichia/ecoli。两者都归到 pichia 组。
    return (h.indexOf("pichia") >= 0 || h.indexOf("pastoris") >= 0 || h.indexOf("yeast") >= 0) ? "pichia" : "ecoli";
  }

  // 数值兜底（与 ferment_sim.js 的 num 同源语义）：NaN/非数 → 默认值
  function num(v, d) { v = parseFloat(v); return (isFinite(v) && v === v) ? v : d; }

  function batchById(bid) {
    for (var i = 0; i < DATA.batches.length; i++) if (DATA.batches[i].batch_id === bid) return DATA.batches[i];
    return null;
  }

  // 批次培养基类型 → FermentSim 培养基预设 key（命中则返回，否则由调用方回退宿主默认）
  function batchMedium(b) {
    var host = hostKeyOf(b.host_organism || "");
    var mt = b.medium_type;
    if (mt && FermentSim.MEDIA[host] && FermentSim.MEDIA[host][mt]) return mt;
    return null;
  }

  // -------------------------------------------------------------------------
  // 全局状态
  // -------------------------------------------------------------------------
  var DATA = null;
  var state = {
    page: "home",
    filter: { host: "", scale: "", success: "" },
    selBatch: null,
    expr: { host: "Pichia_pastoris_X33", promoter: "PAOX1", sp: "alpha_factor", seq: "", cai: "", gc: "" },
    ecoli: { seq: "", strain: "K12", tox: false, secret: false, sp: "PelB", model: "builtin", bias: 0, fbOutcome: 0, fbSol: 70, fbInc: 30, fbYield: 120 },
    scale: {
      src: "5L", tgt: "500L", criterion: "pv", do_set: 30, otr_target: 150, mu: 0.0015,
      chart: { metric: "pv", include: null },
      sweep: { tag: "5L", cmp: null, cmpMetric: "pv" },
      host: "ecoli",
      savedAt: null,
      proc: { X: 100, mu: 0.20, muCrit: 0.20, qO2: 8, overflowThr: 1.5, pco2Thr: 0.4,
        shearTipThr: 7, tmixThr: 30, coolU: 500, coolDT: 15, presMax: 1.5,
        corr: { kla: 1, tmix: 1, heat: 1 } },
      fermenters: buildDefaultFermenters()
    },
    copilot: "",
    navFoldYeast: true,
    navFoldET: true,
    process: { host: "all", ref: null, do: 30, temp: 30, feed: 1.0, dur: 72, _lastBid: null,
      params: { DO: true, pH: true, OUR: true, CER: true, DCW: true,
        TEMP: false, RPM: false, AIR: false, O2_FLOW: false, WCW: false, FEED: false },
      cmp: { param: "DCW", batches: [] },
      wipar: { DCW: true, TITER: true, DO: false, pH: false, TEMP: false, RPM: false, AIR: false, O2_FLOW: false, OUR: false, CER: false, WCW: false, FEED: false },
      exp: {
        host_type: "yeast",
        strain: "",
        vessel: "5L",
        fill_kg: 3.0,
        od0: 0.5, mu_batch: 0.18, mu_fed: 0.10, target_od600: 30, induction_trigger: "time",
        do_ctrl: { air_vvm: 0.5, do_set: 30, rpm_min: 200, rpm_max: 1200, o2_start: 25 },
        phases: { batch: true, fedbatch: false, induction: false },
        batch: { dur: 12, ph: 6.8, temp: 30, feed_mode: "none", feed_const: 0, feed_steps: [], c_depletion_rate: 1.5,
          c_source: "甘油", c_conc: 20, n_source: "蛋白胨/酵母浸出物", cn_ratio: "8:1", recipe: [] },
        fedbatch: { dur: 24, ph: 6.8, temp: 30, feed_mode: "constant", feed_const: 1.0, feed_steps: [{ t: 0, rate: 0.5 }, { t: 12, rate: 1.2 }],
          c_source: "甘油", c_conc: 500, n_source: "酵母浸出物", cn_ratio: "6:1", recipe: [] },
        induction: { dur: 20, ph: 6.8, temp: 25, feed_mode: "constant", feed_const: 1.0, feed_steps: [{ t: 0, rate: 1.0 }],
          c_source: "甲醇", c_conc: 792, n_source: "—", cn_ratio: "—", recipe: [] }
      }
    },
  };

  var PRIMARY = "#1f9e89", ACCENT = "#2b6cb0", BG = "#0c1216", PANEL = "#13202a";

  // -------------------------------------------------------------------------
  // 小工具
  // -------------------------------------------------------------------------
  function $(id) { return document.getElementById(id); }
  function fmt(x, d) { if (x == null || isNaN(x)) return "—"; d = (d == null) ? 1 : d; return Number(x).toFixed(d); }
  function hostColor(h) { return (h && h.indexOf("E.coli") >= 0) ? ACCENT : "#b06ab3"; }
  function esc(s) { return String(s == null ? "" : s).replace(/[&<>"]/g, function (c) { return { "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[c]; }); }

  // -------------------------------------------------------------------------
  // SVG 图表
  // -------------------------------------------------------------------------
  function svgNS(tag) { return document.createElementNS("http://www.w3.org/2000/svg", tag); }

  // 多系列折线图
  function lineChart(host, series, opts) {
    opts = opts || {};
    var W = opts.w || 680, H = opts.h || 320;
    var padL = 56, padR = 18, padT = 14, padB = 48;   // 左/右/上/下留白：容纳轴标题与刻度
    var allx = [], ally = [];
    series.forEach(function (s) { s.points.forEach(function (p) { allx.push(p[0]); ally.push(p[1]); }); });
    var xmin = Math.min.apply(null, allx), xmax = Math.max.apply(null, allx);
    var ymin = opts.ymin != null ? opts.ymin : Math.min.apply(null, ally);
    var ymax = opts.ymax != null ? opts.ymax : Math.max.apply(null, ally);
    if (ymin === ymax) { ymax = ymin + 1; }
    if (xmin === xmax) { xmax = xmin + 1; }
    // 0 基线：若含负值或 ymin>0 不让 0 混入导致压缩；仅当 ymin<0 时纳入
    var ylo = ymin, yhi = ymax;
    function X(v) { return padL + (v - xmin) / (xmax - xmin) * (W - padL - padR); }
    function Y(v) { return H - padB - (v - ylo) / (yhi - ylo) * (H - padT - padB); }
    var sb = [];
    sb.push('<svg viewBox="0 0 ' + W + ' ' + H + '" width="100%" height="' + H + '" preserveAspectRatio="xMidYMid meet" style="background:' + BG + ';border-radius:10px">');
    // 网格 + y 轴刻度标签
    for (var i = 0; i <= 4; i++) {
      var yv = ylo + (yhi - ylo) * i / 4;
      var yy = Y(yv);
      sb.push('<line x1="' + padL + '" y1="' + yy.toFixed(1) + '" x2="' + (W - padR) + '" y2="' + yy.toFixed(1) + '" stroke="#1c2b33" stroke-width="1"/>');
      sb.push('<text x="' + (padL - 8) + '" y="' + (yy + 4).toFixed(1) + '" fill="#8fa3b8" font-size="11" text-anchor="end">' + fmt(yv, 1) + '</text>');
    }
    // x 轴网格 + 刻度标签（5 个均匀刻度）
    for (var j = 0; j <= 4; j++) {
      var xv = xmin + (xmax - xmin) * j / 4;
      var xx = X(xv);
      sb.push('<line x1="' + xx.toFixed(1) + '" y1="' + (H - padB) + '" x2="' + xx.toFixed(1) + '" y2="' + (H - padB + 4) + '" stroke="#2a3a44" stroke-width="1"/>');
      sb.push('<text x="' + xx.toFixed(1) + '" y="' + (H - padB + 18) + '" fill="#8fa3b8" font-size="11" text-anchor="middle">' + fmt(xv, 0) + '</text>');
    }
    // 坐标轴主线（左 x、底 y）
    sb.push('<line x1="' + padL + '" y1="' + padT + '" x2="' + padL + '" y2="' + (H - padB) + '" stroke="#3a4d59" stroke-width="1.5"/>');
    sb.push('<line x1="' + padL + '" y1="' + (H - padB) + '" x2="' + (W - padR) + '" y2="' + (H - padB) + '" stroke="#3a4d59" stroke-width="1.5"/>');
    // 工艺分期背景带（batch / fed-batch / induction）
    if (opts.phases && opts.phases.length) {
      opts.phases.forEach(function (ph) {
        var x0 = X(ph.start), x1 = X(ph.end);
        if (x1 - x0 < 0.5) x1 = x0 + 0.5;
        sb.push('<rect x="' + x0.toFixed(1) + '" y="' + padT + '" width="' + (x1 - x0).toFixed(1) + '" height="' + (H - padT - padB).toFixed(1) + '" fill="' + ph.color + '" fill-opacity="0.10"/>');
        sb.push('<line x1="' + x0.toFixed(1) + '" y1="' + padT + '" x2="' + x0.toFixed(1) + '" y2="' + (H - padB) + '" stroke="' + ph.color + '" stroke-width="1" stroke-dasharray="3 3" stroke-opacity="0.5"/>');
        if (x1 - x0 > 34) sb.push('<text x="' + (x0 + 4).toFixed(1) + '" y="' + (padT + 12) + '" fill="' + ph.color + '" font-size="10.5" opacity="0.9">' + esc(ph.label) + '</text>');
      });
    }
    // 可选竖线标记（如 What-if 设定发酵时长）
    if (opts.vline && opts.vline.x != null) {
      var vx = X(opts.vline.x);
      sb.push('<line x1="' + vx.toFixed(1) + '" y1="' + padT + '" x2="' + vx.toFixed(1) + '" y2="' + (H - padB) + '" stroke="' + (opts.vline.color || "#ff6b6b") + '" stroke-width="1.5" stroke-dasharray="4 3"/>');
      if (opts.vline.label) sb.push('<text x="' + vx.toFixed(1) + '" y="' + (padT - 2) + '" fill="' + (opts.vline.color || "#ff6b6b") + '" font-size="10.5" text-anchor="middle">' + esc(opts.vline.label) + '</text>');
    }
    // 轴标题
    sb.push('<text x="' + (padL + (W - padL - padR) / 2) + '" y="' + (H - 6) + '" fill="#b7c7d6" font-size="12" text-anchor="middle">' + (opts.xtitle || "时间 (h)") + '</text>');
    sb.push('<text x="16" y="' + (padT + (H - padT - padB) / 2) + '" fill="#b7c7d6" font-size="12" text-anchor="middle" transform="rotate(-90 16 ' + (padT + (H - padT - padB) / 2) + ')">' + (opts.ytitle || "数值") + '</text>');
    // 数据线
    series.forEach(function (s) {
      var d = "";
      s.points.forEach(function (p, i) { d += (i === 0 ? "M" : "L") + X(p[0]).toFixed(1) + " " + Y(p[1]).toFixed(1) + " "; });
      sb.push('<path d="' + d + '" fill="none" stroke="' + (s.color || PRIMARY) + '" stroke-width="2" stroke-linejoin="round"/>');
    });
    sb.push("</svg>");
    // 图例（HTML 行，色块 + 名称）
    var legend = '<div style="display:flex;gap:16px;flex-wrap:wrap;margin:8px 2px 2px;font-size:12px;color:#cdd9e5">';
    series.forEach(function (s) {
      legend += '<span style="display:inline-flex;align-items:center;gap:6px"><span style="width:14px;height:3px;border-radius:2px;background:' + (s.color || PRIMARY) + ';display:inline-block"></span>' + esc(s.name || "") + (s.unit ? ' <span style="color:#7d93a6">(' + esc(s.unit) + ')</span>' : '') + '</span>';
    });
    legend += "</div>";
    host.innerHTML = legend + sb.join("");
  }

  // 多参数合并曲线：每个参数一张独立纵坐标（真实量纲），左右分列，曲线按各自量程映射
  function lineChartMulti(host, series, opts) {
    opts = opts || {};
    var n = series.length;
    if (!n) { host.innerHTML = '<div style="color:#7d93a6;font-size:12px">未选择参数</div>'; return; }
    var colW = 28, leftN = Math.ceil(n / 2), rightN = n - leftN;
    var padT = 16, padB = 46, padL = 30 + leftN * colW, padR = 26 + rightN * colW;
    var innerW = 600, W = padL + innerW + padR, H = 360;
    var x0 = padL, x1 = W - padR, y0 = padT, y1 = H - padB;
    var hrs = series[0].points.map(function (p) { return p[0]; });
    var hmin = Math.min.apply(null, hrs), hmax = Math.max.apply(null, hrs);
    if (hmin === hmax) hmax = hmin + 1;
    function X(v) { return x0 + (v - hmin) / (hmax - hmin) * (x1 - x0); }
    function Yy(t) { return y1 - t * (y1 - y0); }
    var sb = [];
    sb.push('<svg viewBox="0 0 ' + W + ' ' + H + '" width="100%" height="' + H + '" preserveAspectRatio="xMidYMid meet" style="background:' + BG + ';border-radius:10px">');
    for (var g = 0; g <= 4; g++) {
      var gy = Yy(g / 4);
      sb.push('<line x1="' + x0 + '" y1="' + gy.toFixed(1) + '" x2="' + x1 + '" y2="' + gy.toFixed(1) + '" stroke="#17242c" stroke-width="1"/>');
    }
    if (opts.phases && opts.phases.length) {
      opts.phases.forEach(function (ph) {
        var xa = X(ph.start), xb = X(ph.end); if (xb - xa < 0.5) xb = xa + 0.5;
        sb.push('<rect x="' + xa.toFixed(1) + '" y="' + y0 + '" width="' + (xb - xa).toFixed(1) + '" height="' + (y1 - y0).toFixed(1) + '" fill="' + ph.color + '" fill-opacity="0.10"/>');
        sb.push('<line x1="' + xa.toFixed(1) + '" y1="' + y0 + '" x2="' + xa.toFixed(1) + '" y2="' + y1 + '" stroke="' + ph.color + '" stroke-width="1" stroke-dasharray="3 3" stroke-opacity="0.5"/>');
        if (xb - xa > 34) sb.push('<text x="' + (xa + 4).toFixed(1) + '" y="' + (y0 + 12) + '" fill="' + ph.color + '" font-size="10.5" opacity="0.9">' + esc(ph.label) + '</text>');
      });
    }
    for (var j = 0; j <= 4; j++) {
      var xv = hmin + (hmax - hmin) * j / 4, xx = X(xv);
      sb.push('<line x1="' + xx.toFixed(1) + '" y1="' + y1 + '" x2="' + xx.toFixed(1) + '" y2="' + (y1 + 4) + '" stroke="#2a3a44" stroke-width="1"/>');
      sb.push('<text x="' + xx.toFixed(1) + '" y="' + (y1 + 18) + '" fill="#8fa3b8" font-size="11" text-anchor="middle">' + fmt(xv, 0) + '</text>');
    }
    // 每个参数独立纵坐标
    series.forEach(function (s, i) {
      var lo = s.ylo, hi = s.yhi; if (hi === lo) hi = lo + 1;
      var left = i < leftN, k = left ? i : (i - leftN);
      var axX = left ? (x0 - (k + 1) * colW) : (x1 + (k + 1) * colW);
      sb.push('<line x1="' + axX.toFixed(1) + '" y1="' + y0 + '" x2="' + axX.toFixed(1) + '" y2="' + y1 + '" stroke="' + (s.color || PRIMARY) + '" stroke-width="1.4" stroke-opacity="0.85"/>');
      for (var t = 0; t <= 4; t++) {
        var tv = t / 4, yv = lo + (hi - lo) * tv, ty = Yy(tv);
        sb.push('<line x1="' + (axX - 2.5) + '" y1="' + ty.toFixed(1) + '" x2="' + (axX + 2.5) + '" y2="' + ty.toFixed(1) + '" stroke="' + (s.color || PRIMARY) + '" stroke-width="1" stroke-opacity="0.8"/>');
        sb.push('<text x="' + (left ? (axX - 4) : (axX + 4)).toFixed(1) + '" y="' + (ty + 3.5).toFixed(1) + '" fill="' + (s.color || PRIMARY) + '" font-size="9" text-anchor="' + (left ? "end" : "start") + '">' + fmt(yv, 1) + '</text>');
      }
      var midY = (y0 + y1) / 2;
      sb.push('<text x="' + (axX + (left ? -3 : 3)).toFixed(1) + '" y="' + midY.toFixed(1) + '" fill="' + (s.color || PRIMARY) + '" font-size="10" text-anchor="middle" transform="rotate(' + (left ? -90 : 90) + ' ' + axX.toFixed(1) + ' ' + midY.toFixed(1) + ')">' + esc(s.name) + (s.unit ? '(' + esc(s.unit) + ')' : '') + '</text>');
    });
    sb.push('<line x1="' + x0 + '" y1="' + y1 + '" x2="' + x1 + '" y2="' + y1 + '" stroke="#3a4d59" stroke-width="1.5"/>');
    sb.push('<text x="' + ((x0 + x1) / 2) + '" y="' + (H - 6) + '" fill="#b7c7d6" font-size="12" text-anchor="middle">' + (opts.xtitle || "时间 (h)") + '</text>');
    series.forEach(function (s) {
      var lo = s.ylo, hi = s.yhi; if (hi === lo) hi = lo + 1;
      var d = "";
      s.points.forEach(function (p, i2) { var t = (p[1] - lo) / (hi - lo); t = Math.max(0, Math.min(1, t)); d += (i2 === 0 ? "M" : "L") + X(p[0]).toFixed(1) + " " + Yy(t).toFixed(1) + " "; });
      sb.push('<path d="' + d + '" fill="none" stroke="' + (s.color || PRIMARY) + '" stroke-width="2" stroke-linejoin="round"' + (s.dash ? ' stroke-dasharray="6 4" stroke-opacity="0.9"' : '') + '/>');
    });
    sb.push("</svg>");
    var legend = '<div style="display:flex;gap:14px;flex-wrap:wrap;margin:8px 2px 2px;font-size:12px;color:#cdd9e5">';
    series.forEach(function (s) {
      var lo = s.ylo, hi = s.yhi, last = s.points[s.points.length - 1][1];
      legend += '<span style="display:inline-flex;align-items:center;gap:6px"><span style="width:14px;height:3px;border-radius:2px;background:' + (s.color || PRIMARY) + ';display:inline-block"></span>' + esc(s.name || "") + (s.unit ? ' <span style="color:#7d93a6">(' + esc(s.unit) + ')</span>' : '') + ' <span style="color:#9fb3c8">量程 ' + fmt(lo, 1) + '–' + fmt(hi, 1) + ' · 当前 ' + fmt(last, 1) + '</span></span>';
    });
    legend += "</div>";
    host.innerHTML = legend + sb.join("");
  }

  // 柱状图
  function barChart(host, items, opts) {
    opts = opts || {};
    var W = opts.w || 360, H = opts.h || 260, pad = 36;
    var maxv = Math.max.apply(null, items.map(function (x) { return x.value; })) || 1;
    var n = items.length, bw = (W - pad * 2) / n * 0.6;
    var gap = (W - pad * 2) / n;
    var sb = ['<svg viewBox="0 0 ' + W + ' ' + H + '" width="100%" style="background:' + BG + ';border-radius:10px">'];
    items.forEach(function (it, i) {
      var h = (it.value / maxv) * (H - pad * 2);
      var x = pad + gap * i + (gap - bw) / 2;
      var y = H - pad - h;
      sb.push('<rect x="' + x.toFixed(1) + '" y="' + y.toFixed(1) + '" width="' + bw.toFixed(1) + '" height="' + h.toFixed(1) + '" rx="4" fill="' + (it.color || PRIMARY) + '"/>');
      sb.push('<text x="' + (x + bw / 2).toFixed(1) + '" y="' + (y - 6).toFixed(1) + '" fill="#e6edf3" font-size="11" text-anchor="middle">' + fmt(it.value, 0) + '</text>');
      sb.push('<text x="' + (x + bw / 2).toFixed(1) + '" y="' + (H - pad + 14).toFixed(1) + '" fill="#7d93a6" font-size="11" text-anchor="middle">' + esc(it.label) + '</text>');
    });
    sb.push("</svg>");
    host.innerHTML = sb.join("");
  }

  // 雷达图 (0..1)
  function radarChart(host, axes, values) {
    var W = 420, H = 320, cx = W / 2, cy = H / 2 + 6, R = 110, n = axes.length;
    var sb = ['<svg viewBox="0 0 ' + W + ' ' + H + '" width="100%" style="background:' + BG + ';border-radius:10px">'];
    // 网格环
    for (var g = 1; g <= 4; g++) {
      var pts = [];
      for (var i = 0; i < n; i++) {
        var a = -Math.PI / 2 + 2 * Math.PI * i / n;
        var r = R * g / 4;
        pts.push((cx + r * Math.cos(a)).toFixed(1) + "," + (cy + r * Math.sin(a)).toFixed(1));
      }
      sb.push('<polygon points="' + pts.join(" ") + '" fill="none" stroke="#1c2b33" stroke-width="1"/>');
    }
    // 轴 + 标签
    for (var j = 0; j < n; j++) {
      var aj = -Math.PI / 2 + 2 * Math.PI * j / n;
      var ex = cx + R * Math.cos(aj), ey = cy + R * Math.sin(aj);
      sb.push('<line x1="' + cx + '" y1="' + cy + '" x2="' + ex.toFixed(1) + '" y2="' + ey.toFixed(1) + '" stroke="#1c2b33" stroke-width="1"/>');
      var lx = cx + (R + 26) * Math.cos(aj), ly = cy + (R + 22) * Math.sin(aj);
      var anchor = Math.abs(Math.cos(aj)) < 0.3 ? "middle" : (Math.cos(aj) > 0 ? "start" : "end");
      sb.push('<text x="' + lx.toFixed(1) + '" y="' + ly.toFixed(1) + '" fill="#9fb3c8" font-size="11" text-anchor="' + anchor + '">' + esc(axes[j]) + '</text>');
    }
    // 数据多边形
    var dpts = [];
    for (var k = 0; k < n; k++) {
      var ak = -Math.PI / 2 + 2 * Math.PI * k / n;
      var v = Math.max(0, Math.min(1, values[k]));
      dpts.push((cx + R * v * Math.cos(ak)).toFixed(1) + "," + (cy + R * v * Math.sin(ak)).toFixed(1));
    }
    sb.push('<polygon points="' + dpts.join(" ") + '" fill="' + PRIMARY + '55" stroke="' + PRIMARY + '" stroke-width="2"/>');
    for (var m = 0; m < n; m++) {
      var am = -Math.PI / 2 + 2 * Math.PI * m / n;
      var vm = Math.max(0, Math.min(1, values[m]));
      sb.push('<circle cx="' + (cx + R * vm * Math.cos(am)).toFixed(1) + '" cy="' + (cy + R * vm * Math.sin(am)).toFixed(1) + '" r="3" fill="' + PRIMARY + '"/>');
    }
    sb.push("</svg>");
    host.innerHTML = sb.join("");
  }

  // 热力图
  function heatmap(host, rows, cols, matrix) {
    var W = 380, H = 80 + rows.length * 34, cellW = 90, padL = 90;
    var sb = ['<svg viewBox="0 0 ' + W + ' ' + H + '" width="100%" style="background:' + BG + ';border-radius:10px">'];
    sb.push('<text x="6" y="22" fill="#9fb3c8" font-size="12">宿主 \\ 尺度</text>');
    cols.forEach(function (c, j) {
      sb.push('<text x="' + (padL + j * cellW + cellW / 2) + '" y="22" fill="#9fb3c8" font-size="12" text-anchor="middle">' + esc(c) + '</text>');
    });
    rows.forEach(function (r, i) {
      sb.push('<text x="6" y="' + (44 + i * 34 + 18) + '" fill="#e6edf3" font-size="12">' + esc(r) + '</text>');
      cols.forEach(function (c, j) {
        var v = matrix[i][j] || 0;
        var maxv = Math.max.apply(null, matrix.flat()) || 1;
        var t = v / maxv;
        var col = t > 0 ? "rgb(" + Math.round(30 + t * 0) + "," + Math.round(120 + t * 100) + "," + Math.round(120 + t * 40) + ")" : "#16242c";
        sb.push('<rect x="' + (padL + j * cellW) + '" y="' + (44 + i * 34) + '" width="' + (cellW - 6) + '" height="28" rx="4" fill="' + (v ? "#1f9e89" : "#16242c") + '" opacity="' + (0.35 + 0.65 * t) + '"/>');
        sb.push('<text x="' + (padL + j * cellW + (cellW - 6) / 2) + '" y="' + (44 + i * 34 + 18) + '" fill="#e6edf3" font-size="12" text-anchor="middle">' + (v || "") + '</text>');
      });
    });
    sb.push("</svg>");
    host.innerHTML = sb.join("");
  }

  // -------------------------------------------------------------------------
  // Expression Twin · 7 层评分 (移植自 fip/twins/expression/yeast.py)
  // -------------------------------------------------------------------------
  var KD = { A: 1.8, R: -4.5, N: -3.5, D: -3.5, C: 2.5, Q: -3.5, E: -3.5, G: -0.4, H: -3.2, I: 4.5, L: 3.8, K: -3.9, M: 1.9, F: 2.8, P: -1.6, S: -0.8, T: -0.7, W: -0.9, Y: -1.3, V: 4.2 };
  var PROMOTER_STRENGTH = { PAOX1: 0.95, PGAP: 0.70, GAP: 0.70, TEF1: 0.65, GAL1: 0.60, T7lac: 0.90, tac: 0.85, trc: 0.85 };
  var SIGNAL_PEPTIDE_SCORE = { "alpha-MF": 0.95, "alpha-factor": 0.92, "Ost1": 0.88, "PHO5-SP": 0.85, "GlaA-SP": 0.82, "SUC2-SP": 0.80, "Native SP": 0.60, "none": 0 };
  var LAYER_WEIGHTS = { L1_transcription: 0.15, L2_translation: 0.15, L3_folding: 0.15, L4_er: 0.15, L5_glycosylation: 0.10, L6_secretion: 0.20, L7_full: 0.10 };
  var DISORDER_RES = "GSQPNKRTYE";

  function seqFeatures(seq) {
    var s = "", i;
    for (i = 0; i < seq.length; i++) { var c = seq[i].toUpperCase(); if (KD[c] !== undefined) s += c; }
    var n = Math.max(s.length, 1);
    var gravy = 0, aro = 0, cys = 0, pro = 0, chg = 0, dis = 0, k;
    for (k = 0; k < s.length; k++) { var a = s[k]; gravy += (KD[a] || 0); if ("FWY".indexOf(a) >= 0) aro++; if (a === "C") cys++; if (a === "P") pro++; if (a === "K" || a === "R") chg++; if (a === "D" || a === "E") chg--; if (DISORDER_RES.indexOf(a) >= 0) dis++; }
    return { length: s.length, mw_kda: Math.round(s.length * 110.0 / 10) / 10, gravy: Math.round(gravy / n * 1000) / 1000, aromatic_frac: Math.round(aro / n * 1000) / 1000, cys_count: cys, pro_frac: Math.round(pro / n * 1000) / 1000, net_charge: chg, n_glyco_sites: nGlyco(s).length, disorder_frac: Math.round(dis / n * 1000) / 1000 };
  }
  function nGlyco(seq) {
    var s = "", i; for (i = 0; i < seq.length; i++) { var c = seq[i].toUpperCase(); if (KD[c] !== undefined) s += c; }
    var out = [];
    for (i = 0; i < s.length - 2; i++) { if (s[i] === "N" && s[i + 1] !== "P" && "ST".indexOf(s[i + 2]) >= 0) out.push(s.substr(i, 3)); }
    return out;
  }
  function layerTranscription(promoter) { var s = PROMOTER_STRENGTH[promoter] != null ? PROMOTER_STRENGTH[promoter] : 0.55; return [s, "启动子 " + promoter + "：强度代理 " + s.toFixed(2) + "。PAOX1 强但需甲醇诱导；PGAP 组成型。真实评分建议用 ExpressYeaself。"]; }
  function layerTranslation(cai, gc) { cai = cai == null || cai === "" ? 0.5 : parseFloat(cai); gc = gc == null || gc === "" ? 0.5 : parseFloat(gc); var pen = Math.abs(gc - 0.5); var s = Math.max(0, Math.min(1, 0.7 * cai + 0.3 * (1 - pen))); return [s, "CAI=" + cai.toFixed(2) + "、GC=" + gc.toFixed(2) + "：翻译效率代理 " + s.toFixed(2) + "。"]; }
  function layerFolding(f) { var gravy = f.gravy || 0, length = f.length || 300, cys = f.cys_count || 0, pro = f.pro_frac || 0; var s = 1.0 - 0.5 * Math.max(0, gravy) - 0.0003 * Math.max(0, length - 300) - 0.03 * cys - 0.2 * pro; s = Math.max(0.05, Math.min(1, s)); return [s, "GRAVY=" + gravy.toFixed(2) + "、长度=" + length + "、Cys=" + cys + "：折叠难度代理 " + s.toFixed(2) + "。真实结构风险建议用 AlphaFold pLDDT。"]; }
  function layerER(host, sp, chap) { var sps = SIGNAL_PEPTIDE_SCORE[sp] != null ? SIGNAL_PEPTIDE_SCORE[sp] : 0.10; var isP = host && host.indexOf("Pichia") >= 0; var c = (chap && chap !== "none" && chap !== "None" && chap !== "") ? 0.15 : 0; var s = Math.max(0, Math.min(1, sps + c)); return [s, "信号肽=" + sp + "(" + sps.toFixed(2) + ") + 伴侣=" + (chap || "none") + "(" + c.toFixed(2) + ")：ER 转位代理 " + s.toFixed(2) + (isP ? "。毕赤酵母依赖 alpha-MF/SUC2 进入分泌途径。" : "")]; }
  function layerGlyco(f, host) { var nn = f.n_glyco_sites || 0; var isP = host && host.indexOf("Pichia") >= 0; var s; if (isP) { s = nn === 0 ? 0.85 : nn <= 3 ? 0.75 : nn <= 6 ? 0.45 : 0.20; } else { s = nn <= 3 ? 0.80 : nn <= 6 ? 0.55 : 0.35; } return [s, "N-糖基化位点=" + nn + "：糖基化风险代理 " + s.toFixed(2) + (isP ? "。Pichia 高甘露糖化风险随位点数上升。" : "")]; }
  function layerSecretion(f, sp, sol) { var sps = SIGNAL_PEPTIDE_SCORE[sp] != null ? SIGNAL_PEPTIDE_SCORE[sp] : 0.10; sol = sol != null ? sol : (1.0 - 0.5 * Math.max(0, f.gravy || 0)); var s = Math.max(0, Math.min(1, 0.6 * sps + 0.4 * sol)); return [s, "信号肽分泌(" + sps.toFixed(2) + ") × 溶解度代理(" + sol.toFixed(2) + ")：分泌效率代理 " + s.toFixed(2) + "。"]; }
  function mapOutcome(overall, host, f) {
    var isP = host && host.indexOf("Pichia") >= 0;
    var ymax = isP ? 3000.0 : 2500.0;
    var y = Math.round(ymax * Math.pow(overall, 2.2) * 10) / 10;
    var lvl = overall >= 0.65 ? "high" : overall >= 0.40 ? "mid" : overall >= 0.20 ? "low" : "none";
    return [lvl, y];
  }
  function sevenLayer(host, promoter, sp, aaSeq, cai, gc, chaperone) {
    var f;
    if (aaSeq && aaSeq.trim()) { f = seqFeatures(aaSeq); } else {
      f = { length: 300, mw_kda: 33.0, gravy: 0.0, aromatic_frac: 0.05, cys_count: 2, pro_frac: 0.05, net_charge: 0, n_glyco_sites: 2, disorder_frac: 0.4 };
    }
    var sol = Math.max(0.1, Math.min(0.95, 1.0 - 0.5 * Math.max(0, f.gravy || 0) - 0.0002 * Math.max(0, (f.length || 300) - 300)));
    var s1 = layerTranscription(promoter), s2 = layerTranslation(cai, gc), s3 = layerFolding(f), s4 = layerER(host, sp, chaperone), s5 = layerGlyco(f, host), s6 = layerSecretion(f, sp, sol);
    var comp = { L1_transcription: s1[0], L2_translation: s2[0], L3_folding: s3[0], L4_er: s4[0], L5_glycosylation: s5[0], L6_secretion: s6[0] };
    var overall = 0; for (var k in comp) overall += LAYER_WEIGHTS[k] * comp[k];
    overall = Math.max(0, Math.min(1, overall));
    var o = mapOutcome(overall, host, f);
    return {
      overall: Math.round(overall * 1000) / 1000, features: f, solubility: Math.round(sol * 1000) / 1000,
      levels: [
        { key: "L1_transcription", name: "转录", score: s1[0], note: s1[1] },
        { key: "L2_translation", name: "翻译", score: s2[0], note: s2[1] },
        { key: "L3_folding", name: "蛋白折叠", score: s3[0], note: s3[1] },
        { key: "L4_er", name: "ER加工", score: s4[0], note: s4[1] },
        { key: "L5_glycosylation", name: "糖基化", score: s5[0], note: s5[1] },
        { key: "L6_secretion", name: "分泌", score: s6[0], note: s6[1] }
      ],
      expression_level: o[0], yield_mg_l: o[1],
      evidence: aaSeq && aaSeq.trim() ? "序列驱动 7 层评分（启发式代理）；接入 ExpressYeaself/SignalP/NetNGlyc/AlphaFold 后可升级为证据级预测。" : "无序列，基于构架特征的平均代理（结构不确定性高）。"
    };
  }

  // -------------------------------------------------------------------------
  // Scale-up Twin · 几何相似放大 (工程代理)
  // -------------------------------------------------------------------------
  function volOfScale(tag) {
    // 典型工作体积 (L)
    if (tag === "S") return 2;
    if (tag === "M") return 30;
    if (tag === "L") return 200;
    return parseFloat(tag) || 2;
  }
  function diamOfVol(V) { // H/D = 3 时 V = (3π/4) D^3
    return Math.pow(4 * V / (3 * Math.PI), 1 / 3);
  }
  function scaleupCalc(srcTag, tgtTag, criterion, srcRPM) {
    srcRPM = srcRPM || 800;
    var Vs = volOfScale(srcTag), Vt = volOfScale(tgtTag);
    var Ds = diamOfVol(Vs), Dt = diamOfVol(Vt);
    var f = Math.pow(Vt / Vs, 1 / 3); // 线性放大因子
    var Nt;
    if (criterion === "pv" || criterion === "kla") { Nt = srcRPM * Math.pow(f, -2 / 3); } // 等 P/V ≡ 等 kLa（同 vvm，二者同源）
    else if (criterion === "tip") { Nt = srcRPM * Math.pow(f, -1); }    // 等桨尖速度
    else { Nt = srcRPM; }                                               // 等混合时间 ≡ 等转速
    var tipS = Math.PI * Dt * Nt / 60; // m/s (Dt 单位 m)
    var pv = Math.pow(Nt, 3) * Math.pow(Dt, 2) / 1000; // 相对 P/V (kW/m3 量级代理)
    var kla = criterion === "kla" ? 1.0 : Math.pow(Nt / srcRPM, 2); // 相对 kLa
    var tmix = 4.5 * (Math.pow(Dt, 2)) / (Nt * Math.pow(Dt, 3) / 60 * 60); // 经验混合时间代理(s)
    if (!isFinite(tmix) || tmix <= 0) tmix = 30;
    // 风险
    var risks = [];
    if (tipS > 7) risks.push({ lvl: "red", t: "桨尖速度 " + tipS.toFixed(1) + " m/s > 7，剪切损伤风险" });
    else if (tipS > 5.5) risks.push({ lvl: "yellow", t: "桨尖速度 " + tipS.toFixed(1) + " m/s，偏高需关注" });
    if (tmix > 120) risks.push({ lvl: "red", t: "混合时间 " + tmix.toFixed(0) + " s 过长，混合不足" });
    else if (tmix > 60) risks.push({ lvl: "yellow", t: "混合时间 " + tmix.toFixed(0) + " s，偏长" });
    if (kla < 0.6) risks.push({ lvl: "red", t: "传氧能力相对下降明显，可能氧限制" });
    // 综合
    var lvl = "green";
    risks.forEach(function (r) { if (r.lvl === "red") lvl = "red"; else if (r.lvl === "yellow" && lvl !== "red") lvl = "yellow"; });
    if (!risks.length) risks.push({ lvl: "green", t: "各工程判据在推荐区间内" });
    return {
      Vs: Vs, Vt: Vt, Ds: Ds, Dt: Dt, factor: f, Nt: Nt, tipS: tipS, pv: pv, kla: kla, tmix: tmix,
      risks: risks, level: lvl
    };
  }

  // -------------------------------------------------------------------------
  // 发酵罐结构数据 / Reactor Engineering · 参数计算引擎
  // 由罐体几何 + 操作条件，计算 功率 P、P/V、桨尖速度、kLa、OTR、混合时间
  // -------------------------------------------------------------------------
  // 桨型 → 功率数 Np 与局部剪切系数 k_s（多桨时 Np 合计=ΣNp_i；剪切取各层最大值）
  function impNp(t) { return t === "hydrofoil" ? 0.5 : (t === "pitched" ? 1.3 : 5.0); }
  function impKShear(t) { return t === "rushton" ? 10 : (t === "pitched" ? 7 : 4); }
  var IMP_TYPES = [["rushton", "Rushton(平直叶)"], ["pitched", "斜叶桨"], ["hydrofoil", "翼型桨"]];
  // 各桨型 Np（功率数）与剪切特性说明：用于分层桨叶选择界面与“Np 值说明”
  var IMP_LEGEND = [
    { k: "rushton", name: "Rushton 6 平直叶", np: 5.0, ks: 10, note: "高功率、强气液分散，局部剪切最大，适合好氧/通气" },
    { k: "pitched", name: "斜叶桨 45°", np: 1.3, ks: 7, note: "折中混合与分散，剪切中等，适用范围广" },
    { k: "hydrofoil", name: "翼型桨", np: 0.5, ks: 4, note: "低功耗、低剪切，适合剪切敏感（菌丝/融合蛋白）体系" }
  ];

  // 由工作体积生成标准几何（H/D=3，罐径 D，液高 H=3D，桨径 Di=D/3，桨数随规模递增）
  function genFermenter(V_L, rpm) {
    var V = V_L / 1000;                                   // m³
    var D = Math.pow(4 * V / (3 * Math.PI), 1 / 3);        // m
    var H = 3 * D, Di = D / 3;
    var n_imp = V_L <= 10 ? 1 : (V_L <= 200 ? 2 : (V_L <= 2000 ? 3 : 4));
    var imps = []; for (var i = 0; i < n_imp; i++) imps.push("rushton");
    var n_min = Math.max(10, Math.round(rpm * 0.25)), n_max = Math.round(rpm * 3);
    return { tag: "", name: "", V_L: V_L, D: +D.toFixed(4), H: +H.toFixed(4), Di: +Di.toFixed(4),
      n_imp: n_imp, imps: imps, n_min: n_min, n_max: n_max, vvm: 1.0, o2_vvm: 0, N: rpm, air_vvm: 0.5 };
  }
  // 标准发酵罐工作体积预设（最大工作体积 L, 默认转速 rpm）；内联以避免初始化顺序依赖
  function buildDefaultFermenters() {
    var PRESETS = [
      [2, 700], [3, 650], [5, 550], [10, 450], [50, 350],
      [100, 300], [200, 240], [500, 180], [2000, 130], [10000, 90]
    ];
    var o = {};
    PRESETS.forEach(function (p) {
      var f = genFermenter(p[0], p[1]);
      f.tag = p[0] + "L"; f.name = p[0] + " L 标准罐";
      o[f.tag] = f;
    });
    return o;
  }


  // 单个发酵罐的派生工程参数
  // f: {V_L,D,H,Di,n_imp,imps[],n_min,n_max,vvm,N,air_vvm}; doSet: 设定 DO(%); mu: 发酵液表观黏度 Pa·s
  function fermenterDerived(f, doSet, mu) {
    var rho = 1100;                       // 发酵液密度 kg/m³
    var V = (f.V_L || 0) / 1000;          // 工作体积 m³
    var D = f.D || 0, Di = f.Di || 0, H = f.H || 0;
    var N = (f.N || 0) / 60;              // rev/s
    var imps = (f.imps && f.imps.length) ? f.imps : [f.imp || "rushton"];   // 分层桨型
    var Np = 0; imps.forEach(function (t) { Np += impNp(t); });              // 多桨 Np 合计
    var air_vvm = (f.vvm != null && f.vvm > 0) ? f.vvm : 1.0;    // 空气底通（默认 1 VVM）
    var o2_vvm = (f.o2_vvm != null && f.o2_vvm > 0) ? f.o2_vvm : 0;  // 纯氧共通 (vvm)
    var vvm_tot = air_vvm + o2_vvm;        // 总通气量（驱动 kLa 的气液分散）
    if (V <= 0 || Di <= 0 || N <= 0) {
      return { V: V, D: D, Di: Di, H: H, N: N, Np: Np, vvm: air_vvm, vvm_tot: vvm_tot, yO2: 0.21,
        air_vvm: air_vvm, o2_vvm: o2_vvm,
        P_W: 0, Pv_kW: 0, tip: 0, kla_s: 0, kla_h: 0, OTR: 0, tmix: 0, gdot: 0, shear: 0, ok: false };
    }
    // 功率（基于功率数）：P = Σ Np_i·ρ·N³·Di⁵ = Np·ρ·N³·Di⁵  (W)，Np 为多桨合计
    var P = Np * rho * Math.pow(N, 3) * Math.pow(Di, 5);
    var Pv_kW = P / V / 1000;             // kW/m³
    // 桨尖速度：v_tip = π·Di·N  (m/s)
    var tip = Math.PI * Di * N;
    // kLa — van't Riet 关联：kLa = K·(P/V)^a·(vvm)^b  （P/V 单位 kW/m³，vvm 单位 vvm 总通气）
    //   非凝聚(水相) K=0.032, a=0.4, b=0.5；凝聚体系 K=0.026, a=0.7, b=0.2（本平台取非凝聚默认）
    var kl = 0.032 * Math.pow(Pv_kW, 0.4) * Math.pow(vvm_tot, 0.5);   // 1/s（总通气 vvm_tot）
    var kla_h = kl * 3600;                // 1/h
    // 气体氧分压分数 y_O2 = (0.21·air + 1.0·o2)/(air+o2)：纯氧共通升高 y_O2、抬升溶氧饱和 C*
    // C*≈0.008 kg/m³ @ 21% O₂ / 1 atm 空气；故 C* = 0.008·(y_O2/0.21)
    var y_O2 = vvm_tot > 0 ? (0.21 * air_vvm + 1.0 * o2_vvm) / vvm_tot : 0.21;
    var Cstar = 0.008 * y_O2 / 0.21;      // kg/m³
    var dC = Cstar * (1 - (doSet || 30) / 100);
    var OTR = kl * dC * 112500;           // mmol O₂/L/h  (1 kg/m³=31.25 mmol/L, ×3600 s/h)
    // 混合时间 — Bates 关联（等径 Rushton）：N·t_mix = 4.5·(T/Di)²·(H/T)；多桨近似除以 cbrt(n_imp)
    var Nimp = Math.max(1, f.n_imp || 1);
    var tmix = 4.5 * Math.pow(D / Di, 2) * (H > 0 ? H / D : 1) / (N * Math.pow(Nimp, 1 / 3)); // s
    if (!isFinite(tmix) || tmix <= 0) tmix = 0;
    // 剪切：最大局部剪切速率 γ̇ ≈ k_s·N（rps），取各层桨叶最大值；k_s：Rushton 10 / 斜叶 7 / 翼型 4
    var kShearMax = 0; imps.forEach(function (t) { kShearMax = Math.max(kShearMax, impKShear(t)); });
    var gdot = kShearMax * N;             // 1/s
    var tau = (mu || 0.0015) * gdot;      // Pa（剪切应力）
    return { V: V, D: D, Di: Di, H: H, N: N, Np: Np, vvm: air_vvm, vvm_tot: vvm_tot, yO2: y_O2,
      air_vvm: air_vvm, o2_vvm: o2_vvm,
      P_W: P, Pv_kW: Pv_kW, tip: tip, kla_s: kl, kla_h: kla_h, OTR: OTR, tmix: tmix, gdot: gdot, shear: tau, ok: true };
  }

  // 发酵罐结构数据模块渲染（可编辑不同规模几何 + 计算列）
  function renderFermenters(host) {
    if (!host) return;
    var sc = state.scale, fers = sc.fermenters;
    var tags = Object.keys(fers).sort(function (a, b) { return fers[a].V_L - fers[b].V_L; });
    var doSet = sc.do_set, otrT = sc.otr_target, mu = sc.mu;
    function impSel(tag, idx, val) {
      return '<select title="Np=' + impNp(val).toFixed(1) + '、k_s=' + impKShear(val) + '" style="font-size:11px;padding:1px 3px" onchange="FIP.scaleSetImp(\'' + tag + '\',' + idx + ',this.value)">' +
        IMP_TYPES.map(function (t) { return '<option value="' + t[0] + '"' + (val === t[0] ? " selected" : "") + '>' + t[1] + '</option>'; }).join("") + '</select>';
    }
    function impCell(tag, f) {
      var arr = (f.imps && f.imps.length) ? f.imps : [f.imp || "rushton"];
      if (arr.length <= 1) return impSel(tag, 0, arr[0]);
      return '<div style="display:flex;flex-direction:column;gap:3px">' + arr.map(function (t, i) {
        return '<div style="display:flex;gap:3px;align-items:center"><span class="ksub" style="font-size:10px;width:13px">' + (i + 1) + '</span>' + impSel(tag, i, t) + '</div>';
      }).join("") + '</div>';
    }
    function row(tag) {
      var f = fers[tag], d = fermenterDerived(f, doSet, mu);
      var inp = function (field, val, step, w) {
        return '<input type="number" step="' + (step || "any") + '" value="' + val + '" style="width:' + (w || 56) + 'px" onchange="FIP.scaleSetFermenter(\'' + tag + '.' + field + '\', this.value)">';
      };
      return '<tr>' +
        '<td style="white-space:nowrap"><b>' + tag + '</b></td>' +
        '<td>' + inp("name", f.name || "", "", 92) + '</td>' +
        '<td>' + inp("V_L", f.V_L, "any", 52) + '</td>' +
        '<td>' + inp("D", f.D, "any", 52) + '</td>' +
        '<td>' + inp("H", f.H, "any", 52) + '</td>' +
        '<td>' + inp("Di", f.Di, "any", 52) + '</td>' +
        '<td>' + inp("n_imp", f.n_imp, "1", 40) + '</td>' +
        '<td>' + impCell(tag, f) + '</td>' +
        '<td style="text-align:right;color:#ffd27f">' + fmt(d.Np, 1) + '</td>' +
        '<td>' + inp("N", f.N, "any", 52) +
        '<div style="margin-top:2px" class="ksub">范围 ' + inp("n_min", f.n_min, "any", 42) + '~' + inp("n_max", f.n_max, "any", 42) + ' rpm</div></td>' +
        '<td>' + inp("vvm", f.vvm, "any", 48) + '<div class="ksub" style="margin-top:2px">空气底通</div></td>' +
        '<td>' + inp("o2_vvm", f.o2_vvm || 0, "any", 44) + '<div class="ksub" style="margin-top:2px">纯氧共通</div></td>' +
        '<td style="text-align:right;color:#7bdff2">' + fmt(d.yO2 * 100, 0) + '</td>' +
        '<td style="text-align:right;color:#ffce4d">' + fmt(d.Pv_kW, 2) + '</td>' +
        '<td style="text-align:right;color:#4cc9f0">' + fmt(d.tip, 2) + '</td>' +
        '<td style="text-align:right;color:#80ed99">' + fmt(d.kla_h, 1) + '</td>' +
        '<td style="text-align:right;color:#5b8def">' + fmt(d.tmix, 0) + '</td>' +
        '<td style="text-align:right;color:#ff8c8c">' + fmt(d.shear, 3) + '</td>' +
        '<td style="text-align:right;color:#b06ab3">' + fmt(d.OTR, 0) + '</td>' +
        '<td><button class="btn2" onclick="FIP.scaleDelFermenter(\'' + tag + '\')">✕</button></td>' +
        '</tr>';
    }
    var legend = '<div style="margin-top:6px">桨型 Np（功率数）与剪切系数 k_s：' +
      IMP_LEGEND.map(function (x) { return '<b>' + esc(x.name) + '</b> Np=' + x.np.toFixed(1) + '、k_s=' + x.ks + '（' + x.note + '）'; }).join('；') +
      '。多桨时 Np 合计 = ΣNp_i，局部剪切取各层桨叶最大值。</div>';
    var modeTxt = fipStorageMode() === "local" ? "浏览器本地 (localStorage，刷新仍保留)" : "本次会话内存（沙箱预览，刷新会丢失）";
    var savedTxt = state.scale.savedAt ? " · 已保存 " + _fmtTime(state.scale.savedAt) : "";
    var warnTxt = fipStorageMode() === "memory" ? " · 建议「导出」存为文件，下次「导入」恢复" : "";
    // —— 液位高度估算器：给定罐体直径 D 与工作体积 V，圆柱几何 H = 4V/(π·D²) ——
    var hc = state.scale.hcalc || (state.scale.hcalc = { D: 0.12, V: 2 });
    var hcD = +hc.D, hcV = +hc.V;
    var hcH = (hcD > 0) ? 4 * (hcV / 1000) / (Math.PI * hcD * hcD) : 0;
    var hcRatio = (hcD > 0) ? hcH / hcD : 0;
    var hcNote = !isFinite(hcRatio) || hcRatio <= 0 ? "请输入有效的罐体直径 D 与工作体积 V" :
      (hcRatio < 1.5 ? "长径比偏低（浅液）：混合 / 传氧受限风险，需更高 vvm 或转速" :
      (hcRatio <= 4 ? "长径比处于常规发酵罐范围（经验 1.5–4，立式搅拌罐常取 2–3.5）" :
      "长径比偏高：液柱压头 / 功率与混合难度增大，注意底部缺氧与剪切"));
    var applyOpts = tags.map(function (t) { return '<option value="' + t + '">' + t + '</option>'; }).join("");
    var calcHtml = '<div class="card" style="background:#0d1620;border:1px solid #1c2b34;border-radius:8px;padding:10px 12px;margin:8px 0">' +
      '<div class="ksub" style="margin-bottom:6px">📏 液位高度估算 · 给定罐体直径 D + 工作体积 V（圆柱几何 H = 4V / (π·D²)，V 由 L 换算为 m³）</div>' +
      '<div style="display:flex;gap:14px;flex-wrap:wrap;align-items:center">' +
      '<label>罐体直径 D (m) <input type="number" step="any" min="0" value="' + fmt(hcD, 3) + '" style="width:84px" oninput="FIP.scaleSetHcalc(\'D\', this.value)"></label>' +
      '<label>工作体积 V (L) <input type="number" step="any" min="0" value="' + fmt(hcV, 1) + '" style="width:84px" oninput="FIP.scaleSetHcalc(\'V\', this.value)"></label>' +
      '<span>→ 液位高度 <b style="color:#4cc9f0;font-size:15px">' + fmt(hcH, 3) + ' m</b> · 长径比 H/D = <b style="color:#ffce4d">' + fmt(hcRatio, 2) + '</b></span>' +
      '</div>' +
      '<div style="display:flex;gap:10px;flex-wrap:wrap;align-items:center;margin-top:8px">' +
      '<span class="hint" style="margin:0">' + hcNote + '</span>' +
      '<span style="margin-left:auto;display:inline-flex;gap:6px;align-items:center">应用到 ' +
      '<select id="hc_apply_tag">' + applyOpts + '</select>' +
      '<button class="btn2" onclick="FIP.scaleApplyHcalc(document.getElementById(\'hc_apply_tag\').value)">应用液高</button></span>' +
      '</div></div>';
    var html = '<div class="sub">🛢 发酵罐结构数据 · 不同规模（几何 + 操作 → 工程参数）</div>' +
      '<div style="display:flex;gap:14px;flex-wrap:wrap;align-items:center;margin-bottom:8px">' +
      '<label>设定 DO (%) <input type="number" step="any" value="' + doSet + '" style="width:60px" onchange="FIP.set(\'scale.do_set\', parseFloat(this.value))"></label>' +
      '<label>目标需氧 OTR (mmol/L/h) <input type="number" step="any" value="' + otrT + '" style="width:72px" onchange="FIP.set(\'scale.otr_target\', parseFloat(this.value))"></label>' +
      '<label>发酵液黏度 μ (Pa·s) <input type="number" step="any" value="' + mu + '" style="width:66px" onchange="FIP.scaleSetMu(parseFloat(this.value))"></label>' +
      '<button class="btn2" onclick="FIP.scaleAddFermenter()">＋ 添加尺度</button>' +
      '<button class="btn2" onclick="FIP.scaleSaveFermenters()">💾 保存结构参数</button>' +
      '<button class="btn2" onclick="FIP.scaleResetFermenters()">↺ 恢复默认</button>' +
      '<button class="btn2" onclick="FIP.scaleExportFermenters()">⬇ 导出</button>' +
      '<button class="btn2" onclick="(function(){var e=document.getElementById(\'s_fermenters_file\'); if(e) e.click();})()">⬆ 导入</button>' +
      '<input type="file" id="s_fermenters_file" accept="application/json,.json" style="display:none" onchange="FIP.scaleImportFermenters(this.files[0]); this.value=\'\';">' +
      '<span class="hint" style="margin:0 0 0 4px">存储：' + modeTxt + savedTxt + warnTxt + '</span>' +
      '</div>' +
      calcHtml +
      '<div style="overflow:auto"><table class="tbl"><tr>' +
      '<th>规模</th><th>名称</th><th>体积 L</th><th>罐径 D m</th><th>液高 H m</th><th>桨径 Di m</th><th>桨数</th><th>桨型<br>(分层 L1..Ln)</th><th>Np<br>合计</th><th>转速 rpm<br>(操作/范围)</th><th>空气底通<br>vvm</th><th>纯氧<br>vvm</th><th>气体O₂<br>%</th>' +
      '<th>P/V<br>kW/m³</th><th>桨尖<br>m/s</th><th>kLa<br>1/h</th><th>t_mix<br>s</th><th>剪切力<br>Pa</th><th>OTR<br>mmol/L/h</th><th></th></tr>' +
      tags.map(row).join("") + '</table></div>' +
      '<div class="method-note">功率 P=ΣNp_i·ρ·N³·Di⁵（Np 多桨合计）；P/V=P/V；桨尖速度=π·Di·N；kLa 采用 van\'t Riet（K=0.032,α=0.4,β=0.5，非凝聚体系），由<b>总通气</b>（空气底通+纯氧）驱动；OTR=kLa·ΔC（ΔC=C*·(1−DO/100)）；<b>C* 随气体氧分压升高</b>：默认 1 VVM 空气底通（≈21% O₂，C*≈8 mg/L），纯氧共通 o2_vvm 抬高气体 O₂ 分数 y_O₂=(0.21·air+1.0·o2)/(air+o2) → C*↑ → OTR↑；总通气 >1.5 VVM 时建议降低空气底通、提高纯氧占比。混合时间 Bates 关联 N·t_mix=4.5·(T/Di)²·(H/T)，多桨除以 ∛n；剪切力 τ=μ·γ̇，γ̇≈k_s·N。改动几何即时重算，并同步驱动下方放大换算、跨反应器对比图与转速扫掠图。' + legend + '</div>';
    host.innerHTML = html;
  }

  // 跨反应器工程参数对比图（指标可切换，反应器可勾选；P/V、剪切力用对数轴）
  var FERMENTER_METRICS = [
    { k: "pv", name: "P/V（体积功率）", unit: "kW/m³", color: "#ffce4d", log: true, get: function (d) { return d.Pv_kW; } },
    { k: "kla", name: "kLa（体积传质系数）", unit: "1/h", color: "#80ed99", log: false, get: function (d) { return d.kla_h; } },
    { k: "tmix", name: "混合时间", unit: "s", color: "#5b8def", log: false, get: function (d) { return d.tmix; } },
    { k: "tip", name: "桨尖速度", unit: "m/s", color: "#4cc9f0", log: false, get: function (d) { return d.tip; } },
    { k: "shear", name: "剪切力（剪切应力）", unit: "Pa", color: "#ff8c8c", log: true, get: function (d) { return d.shear; } }
  ];
  function renderFermenterChart(host) {
    if (!host) return;
    var sc = state.scale, fers = sc.fermenters;
    var ch = sc.chart || (sc.chart = { metric: "pv", include: null });
    var M = FERMENTER_METRICS[0];
    FERMENTER_METRICS.forEach(function (m) { if (m.k === ch.metric) M = m; });
    var allTags = Object.keys(fers).sort(function (a, b) { return fers[a].V_L - fers[b].V_L; });
    var inc = (ch.include && ch.include.length) ? ch.include : allTags.slice();
    // 指标选择器
    var sel = '<div style="display:flex;gap:8px;flex-wrap:wrap;margin-bottom:8px">';
    FERMENTER_METRICS.forEach(function (m) {
      sel += '<button class="btn2' + (m.k === ch.metric ? " btn2-on" : "") + '" onclick="FIP.scaleSetChartMetric(\'' + m.k + '\')">' + m.name + '</button>';
    });
    sel += '</div>';
    // 反应器勾选
    var chips = '<div style="display:flex;gap:6px;flex-wrap:wrap;align-items:center;margin-bottom:8px"><span class="hint" style="margin:0">对比反应器：</span>';
    allTags.forEach(function (t) {
      var on = inc.indexOf(t) >= 0;
      chips += '<button class="btn2' + (on ? " btn2-on" : "") + '" onclick="FIP.scaleToggleFermenter(\'' + t + '\')">' + t + '</button>';
    });
    chips += '<button class="btn2" onclick="FIP.scaleChartAll()">全选</button></div>';

    var data = allTags.filter(function (t) { return inc.indexOf(t) >= 0; }).map(function (t) {
      var d = fermenterDerived(fers[t], sc.do_set, sc.mu);
      return { tag: t, v: M.get(d) };
    });
    if (!data.length) {
      host.innerHTML = '<div class="sub">📊 跨反应器工程参数对比</div>' + sel + chips +
        '<div class="hint">请至少勾选一个反应器进行对比。</div>';
      return;
    }
    // 绘制柱状图（支持对数/线性）
    var W = 720, Hh = 340, padL = 50, padR = 14, padT = 14, padB = 46;
    var innerW = W - padL - padR, innerH = Hh - padT - padB;
    var x0 = padL, y0 = padT, x1 = W - padR, y1 = Hh - padB;
    var vs = data.map(function (d) { return d.v; });
    var vmax = Math.max.apply(null, vs) || 1;
    var useLog = M.log;
    var lmin = 0, lmax = 1, vFloor = 1e-9;
    if (useLog) {
      var posVs = vs.filter(function (v) { return v > 0; });
      vFloor = posVs.length ? Math.min.apply(null, posVs) : 1e-9;
      lmin = Math.floor(Math.log10(vFloor) * 2) / 2;
      lmax = Math.ceil(Math.log10(vmax) * 2) / 2;
      if (lmax - lmin < 0.5) lmax = lmin + 0.5;
    }
    function yMap(v) {
      if (!useLog) { var f = vmax > 0 ? v / vmax : 0; return y1 - f * innerH; }
      var lv = Math.log10(Math.max(v, vFloor));
      var f = (lv - lmin) / (lmax - lmin); f = Math.max(0, Math.min(1, f));
      return y1 - f * innerH;
    }
    var sb = ['<svg viewBox="0 0 ' + W + ' ' + Hh + '" width="100%" style="background:' + BG + ';border-radius:10px">'];
    if (useLog) {
      for (var k = lmin; k <= lmax + 1e-9; k += 0.5) {
        var vv = Math.pow(10, k), yy = y1 - (k - lmin) / (lmax - lmin) * innerH;
        sb.push('<line x1="' + x0 + '" y1="' + yy.toFixed(1) + '" x2="' + x1 + '" y2="' + yy.toFixed(1) + '" stroke="#17242c" stroke-width="1"/>');
        sb.push('<text x="' + (x0 - 6) + '" y="' + (yy + 3.5).toFixed(1) + '" fill="#8fa3b8" font-size="10" text-anchor="end">' + fmt(vv, vv < 0.01 ? 0 : 2) + '</text>');
      }
    } else {
      for (var g = 0; g <= 4; g++) {
        var yy2 = y1 - g / 4 * innerH, tv = vmax * g / 4;
        sb.push('<line x1="' + x0 + '" y1="' + yy2.toFixed(1) + '" x2="' + x1 + '" y2="' + yy2.toFixed(1) + '" stroke="#17242c" stroke-width="1"/>');
        sb.push('<text x="' + (x0 - 6) + '" y="' + (yy2 + 3.5).toFixed(1) + '" fill="#8fa3b8" font-size="10" text-anchor="end">' + fmt(tv, tv < 1 ? 2 : 0) + '</text>');
      }
    }
    sb.push('<line x1="' + x0 + '" y1="' + y0 + '" x2="' + x0 + '" y2="' + y1 + '" stroke="#3a4d59" stroke-width="1.5"/>');
    sb.push('<line x1="' + x0 + '" y1="' + y1 + '" x2="' + x1 + '" y2="' + y1 + '" stroke="#3a4d59" stroke-width="1.5"/>');
    var n = data.length, gap = innerW / n, bw = Math.min(48, gap * 0.62);
    data.forEach(function (d, i) {
      var cx = x0 + gap * (i + 0.5), bx = cx - bw / 2, by = yMap(d.v), bh = y1 - by;
      sb.push('<rect x="' + bx.toFixed(1) + '" y="' + by.toFixed(1) + '" width="' + bw.toFixed(1) + '" height="' + Math.max(0, bh).toFixed(1) + '" rx="4" fill="' + M.color + '"/>');
      sb.push('<text x="' + cx.toFixed(1) + '" y="' + (by - 6).toFixed(1) + '" fill="#e6edf3" font-size="11" text-anchor="middle">' + fmt(d.v, d.v < 0.01 ? 3 : (d.v < 10 ? 2 : 1)) + '</text>');
      sb.push('<text x="' + cx.toFixed(1) + '" y="' + (y1 + 16).toFixed(1) + '" fill="#7d93a6" font-size="11" text-anchor="middle">' + esc(d.tag) + '</text>');
    });
    sb.push('<text x="14" y="' + (y0 + innerH / 2) + '" fill="#b7c7d6" font-size="12" text-anchor="middle" transform="rotate(-90 14 ' + (y0 + innerH / 2) + ')">' + esc(M.name) + ' (' + esc(M.unit) + ')' + (useLog ? ' · 对数' : '') + '</text>');
    sb.push('</svg>');
    var out = '<div class="sub">📊 跨反应器工程参数对比</div>' + sel + chips +
      '<div style="font-size:12px;color:#9fb3c8;margin-bottom:4px">当前指标：<b style="color:' + M.color + '">' + M.name + '</b> · 单位 ' + M.unit + (useLog ? ' · 纵轴对数刻度' : ' · 纵轴线性刻度') + '</div>' +
      sb.join("") +
      '<div class="method-note">按工作体积升序对比各规格反应器；勾选/取消可增减对比对象。P/V、剪切力跨数量级，纵轴采用对数刻度以便同图对比；kLa、混合时间、桨尖速度为线性刻度。剪切力 τ=μ·γ̇，随黏度 μ 与桨型变化。</div>';
    host.innerHTML = out;
  }

  // 搅拌转速扫掠：选定反应器在 [n_min, n_max] 转速范围内，各工程参数随转速变化的折线图
  var SWEEP_METRICS = [
    { k: "pv", name: "P/V（体积功率）", unit: "kW/m³", color: "#ffce4d", log: true, get: function (d) { return d.Pv_kW; } },
    { k: "tip", name: "桨尖速度", unit: "m/s", color: "#4cc9f0", log: false, get: function (d) { return d.tip; } },
    { k: "kla", name: "kLa（传质系数）", unit: "1/h", color: "#80ed99", log: false, get: function (d) { return d.kla_h; } },
    { k: "tmix", name: "混合时间", unit: "s", color: "#5b8def", log: false, get: function (d) { return d.tmix; } },
    { k: "shear", name: "剪切力（剪切应力）", unit: "Pa", color: "#ff8c8c", log: false, get: function (d) { return d.shear; } },
    { k: "otr", name: "OTR（传氧速率）", unit: "mmol/L/h", color: "#b06ab3", log: false, get: function (d) { return d.OTR; } }
  ];
  function sweepLineSVG(f, M, nmin, nmax, W, Hh) {
    var padL = 46, padR = 8, padT = 10, padB = 24;
    var iw = W - padL - padR, ih = Hh - padT - padB;
    var NP = 44;
    var pts = [];
    for (var i = 0; i < NP; i++) {
      var rpm = nmin + (nmax - nmin) * i / (NP - 1);
      var ff = Object.assign({}, f, { N: Math.round(rpm) });
      var d = fermenterDerived(ff, state.scale.do_set, state.scale.mu);
      pts.push([rpm, M.get(d)]);
    }
    var vmax = Math.max.apply(null, pts.map(function (p) { return p[1]; })) || 1;
    function xMap(rpm) { return padL + (nmax > nmin ? (rpm - nmin) / (nmax - nmin) : 0) * iw; }
    function yMap(v) { var f2 = vmax > 0 ? v / vmax : 0; return (Hh - padB) - f2 * ih; }
    var sb = ['<svg viewBox="0 0 ' + W + ' ' + Hh + '" width="100%" style="background:#0a1116;border-radius:8px">'];
    for (var g = 0; g <= 3; g++) {
      var yy = (Hh - padB) - g / 3 * ih, tv = vmax * g / 3;
      sb.push('<line x1="' + padL + '" y1="' + yy.toFixed(1) + '" x2="' + (W - padR) + '" y2="' + yy.toFixed(1) + '" stroke="#17242c" stroke-width="1"/>');
      sb.push('<text x="' + (padL - 5) + '" y="' + (yy + 3.5).toFixed(1) + '" fill="#8fa3b8" font-size="9" text-anchor="end">' + fmt(tv, tv < 1 ? 2 : 0) + '</text>');
    }
    sb.push('<text x="' + padL + '" y="' + (Hh - 8) + '" fill="#8fa3b8" font-size="9" text-anchor="middle">' + fmt(nmin, 0) + '</text>');
    sb.push('<text x="' + (W - padR) + '" y="' + (Hh - 8) + '" fill="#8fa3b8" font-size="9" text-anchor="middle">' + fmt(nmax, 0) + '</text>');
    var path = "";
    pts.forEach(function (p, i) { path += (i ? "L" : "M") + xMap(p[0]).toFixed(1) + " " + yMap(p[1]).toFixed(1) + " "; });
    sb.push('<path d="' + path + '" fill="none" stroke="' + M.color + '" stroke-width="2"/>');
    pts.forEach(function (p) { sb.push('<circle cx="' + xMap(p[0]).toFixed(1) + '" cy="' + yMap(p[1]).toFixed(1) + '" r="1.6" fill="' + M.color + '"/>'); });
    sb.push('</svg>');
    return sb.join("");
  }
  function renderSweep(host) {
    if (!host) return;
    var sc = state.scale, fers = sc.fermenters;
    var sw = sc.sweep || (sc.sweep = { tag: "5L" });
    if (!fers[sw.tag]) sw.tag = Object.keys(fers).sort(function (a, b) { return fers[a].V_L - fers[b].V_L; })[0];
    var f = fers[sw.tag];
    var nmin = (f.n_min != null && !isNaN(f.n_min)) ? f.n_min : Math.max(10, Math.round(f.N * 0.25));
    var nmax = (f.n_max != null && !isNaN(f.n_max)) ? f.n_max : Math.round(f.N * 3);
    if (nmax <= nmin) nmax = nmin + 10;
    var tagsAll = Object.keys(fers).sort(function (a, b) { return fers[a].V_L - fers[b].V_L; });
    var sel = '<select onchange="FIP.scaleSetSweepTag(this.value)">' + tagsAll.map(function (t) {
      return '<option value="' + t + '"' + (sw.tag === t ? " selected" : "") + '>' + t + '</option>';
    }).join("") + '</select>';
    var W = 300, Hh = 168;
    var grid = '<div style="display:flex;flex-wrap:wrap;gap:12px">';
    SWEEP_METRICS.forEach(function (M) {
      grid += '<div style="flex:1 1 280px;min-width:260px"><div style="font-size:12px;color:' + M.color + ';margin-bottom:2px">' + M.name + ' (' + M.unit + ')</div>' + sweepLineSVG(f, M, nmin, nmax, W, Hh) + '</div>';
    });
    grid += '</div>';
    var html = '<div class="sub">🌀 搅拌转速扫掠 · 各工程参数随转速变化（单罐）</div>' +
      '<div style="display:flex;gap:12px;flex-wrap:wrap;align-items:center;margin-bottom:8px">' +
      '反应器 ' + sel +
      '<span class="hint" style="margin:0">转速范围 <b style="color:#e6edf3">' + fmt(nmin, 0) + '</b> ~ <b style="color:#e6edf3">' + fmt(nmax, 0) + '</b> rpm（来自该反应器“范围”设定，可在上方表格修改）</span>' +
      '</div>' + grid +
      '<div class="method-note">横轴为转速 rpm（' + fmt(nmin, 0) + '→' + fmt(nmax, 0) + '），取自该反应器“转速范围”列；其余几何/操作参数固定。可见 P/V、kLa、OTR、剪切力随转速（≈N³/N²/N）快速上升，桨尖速度（∝N）线性、混合时间（∝1/N）下降。据此选取满足传氧与混合、又不超出剪切/功率上限的操作窗口。</div>';

    // ——— 跨规格转速扫掠对比：多规格同图 ———
    var sw2 = sc.sweep;
    if (!sw2.cmp) sw2.cmp = tagsAll.slice();
    var cmpTags = sw2.cmp.filter(function (t) { return fers[t]; });
    var cmpMetric = sw2.cmpMetric || "pv";
    var cmpM = SWEEP_METRICS[0]; SWEEP_METRICS.forEach(function (m) { if (m.k === cmpMetric) cmpM = m; });
    var cmpSel = '<div style="display:flex;gap:8px;flex-wrap:wrap;margin:6px 0">';
    SWEEP_METRICS.forEach(function (m) {
      cmpSel += '<button class="btn2' + (m.k === cmpMetric ? " btn2-on" : "") + '" onclick="FIP.scaleSweepCmpMetric(\'' + m.k + '\')">' + m.name + '</button>';
    });
    cmpSel += '</div>';
    var cmpChips = '<div style="display:flex;gap:6px;flex-wrap:wrap;align-items:center;margin-bottom:8px"><span class="hint" style="margin:0">对比规格（下拉多选）：</span>';
    tagsAll.forEach(function (t) {
      var on = cmpTags.indexOf(t) >= 0;
      cmpChips += '<button class="btn2' + (on ? " btn2-on" : "") + '" onclick="FIP.scaleSweepCmpToggle(\'' + t + '\')">' + t + '</button>';
    });
    cmpChips += '<button class="btn2" onclick="FIP.scaleSweepCmpAll()">全选</button></div>';
    var cmpChart = sweepCmpSVG(cmpTags, cmpMetric, 760, 360);
    html += '<div class="sub" style="margin-top:14px">🌀 跨规格转速扫掠对比 · 多规格同图（横轴 rpm 对数）</div>' +
      cmpSel + cmpChips +
      '<div style="font-size:12px;color:#9fb3c8;margin-bottom:4px">指标：<b style="color:' + cmpM.color + '">' + cmpM.name + '</b> · 单位 ' + cmpM.unit + '；每条曲线为该规格在其转速范围 [n_min, n_max] 内的扫掠，横轴为对数转速。小体积罐操作转速高（右侧）、大体积罐操作转速低（左侧），可直观比较同转速或同几何下的参数差异。</div>' +
      cmpChart +
      '<div class="method-note">跨规格对比采用共享对数转速轴：各规格曲线仅在其自身工作区间内展开。因体积/几何差异，同一指标在不同规格的“量级”不同；比较时应结合曲线所处转速区与纵轴量值，判断放大会否超出剪切/功率/混合上限。</div>';
    host.innerHTML = html;
  }

  // 跨规格对比调色板（按体积升序分配，保证颜色稳定）
  var CMP_COLORS = ["#ffce4d", "#4cc9f0", "#80ed99", "#ff8c8c", "#b06ab3", "#5b8def",
    "#f7a35c", "#9fe6c8", "#e0aaff", "#ffd166", "#7bdff2", "#ff9b9b"];

  // 跨规格转速扫掠对比图：多规格的扫掠曲线绘于同一张图，横轴为对数转速
  function sweepCmpSVG(tags, metricK, W, Hh) {
    var M = SWEEP_METRICS[0];
    SWEEP_METRICS.forEach(function (m) { if (m.k === metricK) M = m; });
    var fers = state.scale.fermenters, doSet = state.scale.do_set, mu = state.scale.mu;
    var lines = [], rmin = Infinity, rmax = -Infinity, vmax = 0;
    tags.forEach(function (t) {
      var f = fers[t]; if (!f) return;
      var nmin = (f.n_min != null && !isNaN(f.n_min)) ? f.n_min : Math.max(10, Math.round(f.N * 0.25));
      var nmax = (f.n_max != null && !isNaN(f.n_max)) ? f.n_max : Math.round(f.N * 3);
      if (nmax <= nmin) nmax = nmin + 10;
      rmin = Math.min(rmin, nmin); rmax = Math.max(rmax, nmax);
      var pts = [], NP = 40;
      for (var i = 0; i < NP; i++) {
        var rpm = nmin + (nmax - nmin) * i / (NP - 1);
        var ff = Object.assign({}, f, { N: Math.round(rpm) });
        var d = fermenterDerived(ff, doSet, mu);
        var v = M.get(d);
        pts.push([rpm, v]);
        if (v > vmax) vmax = v;
      }
      lines.push({ tag: t, pts: pts });
    });
    if (!lines.length) return '<div class="hint">请至少选择一个对比规格。</div>';
    if (!isFinite(rmin) || !isFinite(rmax) || rmax <= rmin) { rmin = 10; rmax = 2000; }
    vmax = vmax || 1;
    var useLog = M.log;
    var lmin = 0, lmax = 1, vFloor = 1e-9;
    if (useLog) {
      var posVs = [];
      lines.forEach(function (ln) { ln.pts.forEach(function (p) { if (p[1] > 0) posVs.push(p[1]); }); });
      vmax = posVs.length ? Math.max.apply(null, posVs) : 1;
      vFloor = posVs.length ? Math.min.apply(null, posVs) : 1e-9;
      lmin = Math.floor(Math.log10(vFloor) * 2) / 2;
      lmax = Math.ceil(Math.log10(vmax) * 2) / 2;
      if (lmax - lmin < 0.5) lmax = lmin + 0.5;
    }
    var padL = 56, padR = 12, padT = 12, padB = 34;
    var iw = W - padL - padR, ih = Hh - padT - padB;
    function xMap(rpm) { return padL + (rmax > rmin ? (rpm - rmin) / (rmax - rmin) : 0) * iw; }
    function yMap(v) {
      if (!useLog) { var f0 = vmax > 0 ? v / vmax : 0; return (Hh - padB) - f0 * ih; }
      var lv = Math.log10(Math.max(v, vFloor)); var f1 = (lv - lmin) / (lmax - lmin); f1 = Math.max(0, Math.min(1, f1)); return (Hh - padB) - f1 * ih;
    }
    var sb = ['<svg viewBox="0 0 ' + W + ' ' + Hh + '" width="100%" style="background:#0a1116;border-radius:8px">'];
    if (useLog) {
      for (var kk = lmin; kk <= lmax + 1e-9; kk += 0.5) {
        var vv = Math.pow(10, kk), yyL = (Hh - padB) - (kk - lmin) / (lmax - lmin) * ih;
        sb.push('<line x1="' + padL + '" y1="' + yyL.toFixed(1) + '" x2="' + (W - padR) + '" y2="' + yyL.toFixed(1) + '" stroke="#17242c" stroke-width="1"/>');
        sb.push('<text x="' + (padL - 6) + '" y="' + (yyL + 3.5).toFixed(1) + '" fill="#8fa3b8" font-size="10" text-anchor="end">' + fmt(vv, vv < 0.01 ? 0 : (vv < 1 ? 2 : 1)) + '</text>');
      }
    } else {
      for (var g = 0; g <= 4; g++) {
        var yy = (Hh - padB) - g / 4 * ih, tv = vmax * g / 4;
        sb.push('<line x1="' + padL + '" y1="' + yy.toFixed(1) + '" x2="' + (W - padR) + '" y2="' + yy.toFixed(1) + '" stroke="#17242c" stroke-width="1"/>');
        sb.push('<text x="' + (padL - 6) + '" y="' + (yy + 3.5).toFixed(1) + '" fill="#8fa3b8" font-size="10" text-anchor="end">' + fmt(tv, tv < 1 ? 2 : (tv < 100 ? 1 : 0)) + '</text>');
      }
    }
    for (var gx = 0; gx <= 5; gx++) {
      var rv = rmin + (rmax - rmin) * gx / 5;
      var xx = xMap(rv);
      sb.push('<line x1="' + xx.toFixed(1) + '" y1="' + padT + '" x2="' + xx.toFixed(1) + '" y2="' + (Hh - padB) + '" stroke="#10202a" stroke-width="1"/>');
      sb.push('<text x="' + xx.toFixed(1) + '" y="' + (Hh - padB + 14) + '" fill="#8fa3b8" font-size="10" text-anchor="middle">' + fmt(rv, 0) + '</text>');
    }
    sb.push('<line x1="' + padL + '" y1="' + padT + '" x2="' + padL + '" y2="' + (Hh - padB) + '" stroke="#3a4d59" stroke-width="1.5"/>');
    sb.push('<line x1="' + padL + '" y1="' + (Hh - padB) + '" x2="' + (W - padR) + '" y2="' + (Hh - padB) + '" stroke="#3a4d59" stroke-width="1.5"/>');
    lines.forEach(function (ln, idx) {
      var col = CMP_COLORS[idx % CMP_COLORS.length];
      var path = "";
      ln.pts.forEach(function (p, i) { path += (i ? "L" : "M") + xMap(p[0]).toFixed(1) + " " + yMap(p[1]).toFixed(1) + " "; });
      sb.push('<path d="' + path + '" fill="none" stroke="' + col + '" stroke-width="2"/>');
      ln.pts.forEach(function (p) { sb.push('<circle cx="' + xMap(p[0]).toFixed(1) + '" cy="' + yMap(p[1]).toFixed(1) + '" r="1.5" fill="' + col + '"/>'); });
    });
    sb.push('<text x="14" y="' + (padT + ih / 2) + '" fill="#b7c7d6" font-size="11" text-anchor="middle" transform="rotate(-90 14 ' + (padT + ih / 2) + ')">' + esc(M.name) + ' (' + esc(M.unit) + ')' + (useLog ? ' · 对数' : '') + '</text>');
    sb.push('<text x="' + (W - padR) + '" y="' + (padT + 2) + '" fill="#8fa3b8" font-size="10" text-anchor="end">rpm</text>');
    sb.push('</svg>');
    var legend = '<div style="display:flex;gap:12px;flex-wrap:wrap;margin-top:6px">';
    lines.forEach(function (ln, idx) {
      var col = CMP_COLORS[idx % CMP_COLORS.length];
      legend += '<span style="display:inline-flex;align-items:center;gap:5px;font-size:11px;color:#cdd9e5"><span style="width:12px;height:3px;background:' + col + ';border-radius:2px;display:inline-block"></span>' + esc(ln.tag) + '</span>';
    });
    legend += '</div>';
    return sb.join("") + legend;
  }

  // -------------------------------------------------------------------------
  // 发酵罐结构参数持久化（保存 / 读取 / 恢复默认）—— 防止刷新或数据更新后丢失
  // -------------------------------------------------------------------------
  function _fmtTime(iso) {
    try { var d = new Date(iso); var p = function (x) { return (x < 10 ? "0" : "") + x; }; return p(d.getHours()) + ":" + p(d.getMinutes()) + ":" + p(d.getSeconds()); }
    catch (e) { return ""; }
  }
  var FERMENTER_STORE_KEY = "fip_fermenters_v1";
  function scaleFermenterValid(f) {
    if (!f || typeof f !== "object") return false;
    var req = ["V_L", "D", "H", "Di", "n_imp", "N", "n_min", "n_max", "vvm", "air_vvm", "o2_vvm"];
    for (var i = 0; i < req.length; i++) { if (typeof f[req[i]] !== "number" || isNaN(f[req[i]])) return false; }
    if (!Array.isArray(f.imps) || !f.imps.length) return false;
    return true;
  }
  function scaleSanitizeFermenter(f) {
    return {
      tag: f.tag, name: (f.name != null ? String(f.name) : ""),
      V_L: +f.V_L, D: +f.D, H: +f.H, Di: +f.Di, n_imp: +f.n_imp,
      imps: f.imps.slice(), N: +f.N, n_min: +f.n_min, n_max: +f.n_max, vvm: +f.vvm, air_vvm: +f.air_vvm, o2_vvm: +(f.o2_vvm || 0)
    };
  }
  function scaleLoadFermenters() {
    try {
      var raw = fipGet(FERMENTER_STORE_KEY);
      if (!raw) return null;
      var obj = JSON.parse(raw);
      if (!obj || typeof obj !== "object" || !Object.keys(obj).length) return null;
      var out = {}, ok = false;
      Object.keys(obj).forEach(function (t) {
        var f = obj[t];
        if (scaleFermenterValid(f)) { out[t] = scaleSanitizeFermenter(f); ok = true; }
      });
      return ok ? out : null;
    } catch (e) { return null; }
  }
  function scaleSaveFermenters() {
    try {
      var out = {};
      Object.keys(state.scale.fermenters).forEach(function (t) {
        out[t] = scaleSanitizeFermenter(state.scale.fermenters[t]);
      });
      fipSet(FERMENTER_STORE_KEY, JSON.stringify(out));
      state.scale.savedAt = new Date().toISOString();
      return true;
    } catch (e) { return false; }
  }
  function scaleResetFermenters() {
    state.scale.fermenters = buildDefaultFermenters();
    state.scale.savedAt = null;
    fipDel(FERMENTER_STORE_KEY);
  }
  // 导出：把当前各规格发酵罐结构参数下载为 JSON 文件（沙箱预览环境的可靠持久化方案）
  function scaleExportFermenters() {
    try {
      var out = {};
      Object.keys(state.scale.fermenters).forEach(function (t) { out[t] = scaleSanitizeFermenter(state.scale.fermenters[t]); });
      var blob = new Blob([JSON.stringify(out, null, 2)], { type: "application/json" });
      var url = (window.URL || window.webkitURL).createObjectURL(blob);
      var a = document.createElement("a");
      a.href = url; a.download = "fip_fermenters.json";
      document.body.appendChild(a); a.click(); document.body.removeChild(a);
      setTimeout(function () { try { (window.URL || window.webkitURL).revokeObjectURL(url); } catch (e) {} }, 1500);
      return true;
    } catch (e) { return false; }
  }
  // 导入：读取 JSON 文件，校验后恢复发酵罐结构参数，并同步写回存储
  function scaleImportFermenters(file) {
    if (!file) return false;
    var reader = new FileReader();
    reader.onload = function () {
      try {
        var obj = JSON.parse(reader.result);
        if (!obj || typeof obj !== "object" || !Object.keys(obj).length) { alert("导入失败：文件为空或不是有效的发酵罐参数 JSON。"); return; }
        var out = {}, ok = false;
        Object.keys(obj).forEach(function (t) { var f = obj[t]; if (scaleFermenterValid(f)) { out[t] = scaleSanitizeFermenter(f); ok = true; } });
        if (!ok) { alert("导入失败：未找到任何有效发酵罐参数（需含 V_L / D / H / Di / n_imp / 转速范围 / vvm 等数值字段）。"); return; }
        state.scale.fermenters = out;
        scaleSaveFermenters();
        if ($("s_fermenters")) renderFermenters($("s_fermenters"));
        if ($("s_fchart")) renderFermenterChart($("s_fchart"));
        if ($("s_sweep")) renderSweep($("s_sweep"));
        calcScale();
      } catch (e) { alert("导入失败：解析 JSON 出错 — " + (e && e.message ? e.message : e)); }
    };
    reader.readAsText(file);
    return true;
  }

  // -------------------------------------------------------------------------
  // Copilot · 知识库关键词匹配
  // -------------------------------------------------------------------------
  function copilotAnswer(q) {
    q = (q || "").toLowerCase();
    if (!q) return "请输入关于批次、放大、软测量或蛋白表达的问题，例如：“如何提高 Pichia 在 M 到 L 放大的成功率？”。";
    var kb = (DATA && DATA.knowledge) || [];
    var best = [], score = 0;
    kb.forEach(function (it) {
      var text = (it.title + " " + it.content + " " + (it.tags || []).join(" ") + " " + it.category).toLowerCase();
      var sc = 0;
      q.split(/\s+/).forEach(function (w) { if (w.length >= 2 && text.indexOf(w) >= 0) sc += 1; });
      // 中文粗匹配
      ["放大", "成功率", "软测量", "biomass", "表达", "糖基化", "剪切", "混合", "kla", "桨尖", "kpi", "传氧", "批次", "风险", "预测"].forEach(function (kw) { if (q.indexOf(kw) >= 0 && text.indexOf(kw) >= 0) sc += 2; });
      if (sc > score) { score = sc; best = [it]; } else if (sc === score && sc > 0) { best.push(it); }
    });
    if (!best.length || score === 0) {
      return "未在知识库中找到直接匹配。可尝试关键词：放大、成功率、软测量、表达、糖基化、剪切、传氧、KPI。平台当前有 " + (DATA ? DATA.stats.n_batches : 0) + " 个历史批次、" + (DATA ? DATA.knowledge.length : 0) + " 条知识库条目支撑决策。";
    }
    return best.map(function (it) {
      return "【" + esc(it.category) + (it.host ? " · " + esc(it.host) : "") + "】" + esc(it.title) + "\n" + esc(it.content) + (it.source ? "\n来源：" + esc(it.source) : "");
    }).join("\n\n");
  }

  // -------------------------------------------------------------------------
  // 页面渲染
  // -------------------------------------------------------------------------
  function navHTML() {
    var onET = state.page === "yeast" || state.page === "ecoli" || (state.page && state.page.indexOf("yL") === 0);
    var etFolded = !!state.navFoldET && !onET;
    var yeastFolded = !!state.navFoldYeast && !(state.page && state.page.indexOf("yL") === 0);
    var groups = [
      { title: null, items: [["home", "🏠 平台总览"]] },
      { title: "🧬 Expression Twin", nested: [
        { title: "🍶 酵母 Yeast", items: [["yeast", "🍶 酵母总览"]],
          subTitle: "↳ 评估层级 L1–L7", subItems: [
            ["yL1", "↳ L1 转录层"], ["yL2", "↳ L2 翻译层"], ["yL3", "↳ L3 蛋白折叠"],
            ["yL4", "↳ L4 ER加工"], ["yL5", "↳ L5 糖基化"], ["yL6", "↳ L6 分泌"], ["yL7", "↳ L7 全流程"]
          ] },
        { title: "🦠 大肠杆菌 E.coli", items: [["ecoli", "🦠 E.coli 异源表达顾问"]] }
      ] },
      { title: null, items: [["batches", "🗄 批次数据库"], ["process", "📈 Process Twin"], ["scaleup", "⚖️ Scale-up Twin"], ["coupling", "🔗 耦合 (E↔P)"], ["simanim", "🎬 动态仿真 (What-if)"], ["copilot", "💬 Copilot"]] }
    ];
    var html = '<nav style="display:flex;flex-direction:column;gap:4px;padding:10px">';
    function renderGroup(g, isSub) {
      if (g.title) html += '<div class="navsection' + (isSub ? " navsection-sub" : "") + '">' + g.title + '</div>';
      if (g.items) html += g.items.map(function (it) {
        return '<div class="navitem ' + (state.page === it[0] ? "active" : "") + '" onclick="FIP.nav(\'' + it[0] + '\')">' + it[1] + '</div>';
      }).join("");
      if (g.subTitle) {
        var arrow = yeastFolded ? "▸" : "▾";
        html += '<div class="navsubsection" onclick="FIP.navToggle(\'yeast\')">' + arrow + " " + g.subTitle + '</div>';
        html += '<div' + (yeastFolded ? ' style="display:none"' : "") + '>' + g.subItems.map(function (it) {
          return '<div class="navitem navitem-sub ' + (state.page === it[0] ? "active" : "") + '" onclick="FIP.nav(\'' + it[0] + '\')">' + it[1] + '</div>';
        }).join("") + '</div>';
      }
    }
    groups.forEach(function (g) {
      if (g.nested) {
        var arrow = etFolded ? "▸" : "▾";
        html += '<div class="navsection navsection-toggle" onclick="FIP.navToggle(\'et\')">' + arrow + " " + g.title + '</div>';
        if (!etFolded) g.nested.forEach(function (sg) { renderGroup(sg, true); });
      } else {
        renderGroup(g, false);
      }
    });
    return html + "</nav>";
  }

  function renderShell() {
    $("nav").innerHTML = navHTML();
    var titles = { home: "平台总览", yeast: "YeastExpress Pro · 酵母表达评分", ecoli: "E. coli 异源蛋白表达顾问",
      yL1: "L1 转录层 · 酵母表达", yL2: "L2 翻译层 · 酵母表达", yL3: "L3 蛋白折叠 · 酵母表达", yL4: "L4 ER加工 · 酵母表达",
      yL5: "L5 糖基化 · 酵母表达", yL6: "L6 分泌 · 酵母表达", yL7: "L7 全流程预测 · 酵母表达",
      batches: "批次数据库", process: "Process Twin · 工艺孪生", scaleup: "Scale-up Twin · 放大顾问",
      coupling: "Expression↔Process 耦合", simanim: "🎬 动态仿真 (What-if)", copilot: "Fermentation Copilot" };
    $("pagetitle").textContent = titles[state.page] || "";
    var body = $("body-content");
    if (state.page === "home") renderHome(body);
    else if (state.page === "batches") renderBatches(body);
    else if (state.page === "process") renderProcess(body);
    else if (state.page === "yeast") renderExpressionYeast(body);
    else if (state.page === "ecoli") renderExpressionEcoli(body);
    else if (state.page && state.page.indexOf("yL") === 0) renderYeastLevel(parseInt(state.page.slice(2), 10) - 1, body);
    else if (state.page === "scaleup") renderScaleup(body);
    else if (state.page === "coupling") renderCoupling(body);
    else if (state.page === "simanim") renderSimAnim(body);
    else if (state.page === "copilot") renderCopilot(body);
  }

  function card(html, extra) {
    return '<div class="card" ' + (extra || "") + '>' + html + '</div>';
  }

  // ---- 实验设计信息量反馈 (DOE) 面板（阶段 F：交互版 HTML 端口，复刻 Streamlit render_doe_panel）----
  // 数据由 build_interactive.py 预计算进 DATA.doe[host] = { construct, repo, repo_by_scale }，
  // 与 Streamlit 共用同一 DoeInfoAnalyzer；此处仅做只读渲染。
  var DOE_DIM_LABEL = { batch: "批次", construct: "构型", process: "工艺" };
  var DOE_RATING_CLS = { "较充分": "ec-good", "一般": "ec-warn", "偏弱": "ec-bad", "不足": "ec-bad" };

  function doeFactorRows(factors) {
    if (!factors || !factors.length) return '<div class="ksub">无因子数据</div>';
    return '<table class="dev-table"><tr><th>维度</th><th>因子</th><th>占用度</th><th>说明</th></tr>' +
      factors.map(function (f) {
        var occ = f.occupancy == null ? 0 : f.occupancy;
        var cls = occ >= 0.6 ? "dev-ok" : (occ >= 0.35 ? "dev-warn" : "dev-bad");
        return '<tr><td>' + (DOE_DIM_LABEL[f.dimension] || f.dimension) + '</td>' +
          '<td>' + esc(f.label) + '</td>' +
          '<td class="' + cls + '">' + Math.round(occ * 100) + '%</td>' +
          '<td style="color:#9fb3c8">' + esc(f.detail || "") + '</td></tr>';
      }).join("") + '</table>';
  }

  function doeInteractionHtml(interactions) {
    if (!interactions || !interactions.length) return '';
    var rows = interactions.map(function (i) {
      var miss = (i.missing_combos || []).slice(0, 4).join("、") || "—";
      return '<tr><td>' + esc(i.label) + '</td><td>' + i.n_observed_combos + "/" + i.n_possible_combos +
        ' (' + Math.round((i.coverage || 0) * 100) + '%)</td><td style="color:#9fb3c8">' + esc(miss) + '</td></tr>';
    }).join("");
    return '<details style="margin-top:8px"><summary style="cursor:pointer;color:#9fb3c8;font-size:13px">交互覆盖（格子填充率）</summary>' +
      '<table class="dev-table"><tr><th>组合</th><th>覆盖</th><th>缺失组合</th></tr>' + rows + '</table></details>';
  }

  function doeRecHtml(recs) {
    if (!recs || !recs.length) return '';
    return '<details style="margin-top:8px"><summary style="cursor:pointer;color:#9fb3c8;font-size:13px">下一步实验建议</summary>' +
      '<ul class="notes" style="margin:4px 0 4px">' +
      recs.map(function (r) { return '<li>' + esc(r) + '</li>'; }).join("") + '</ul></details>';
  }

  // 单份 DOE 报告卡片（复刻 Streamlit _doe_report_block）
  function doeReportCard(rep, title) {
    if (!rep) return '';
    var ratingCls = DOE_RATING_CLS[rep.rating] || "ec-warn";
    var html = '<div class="sub" style="margin-bottom:6px">' + esc(title) + '</div>';
    html += '<div style="margin:2px 0 8px"><span class="ec-badge ' + ratingCls + '">' + esc(rep.rating || "—") + '</span> ' +
      '信息量 <b style="color:#1f9e89">' + (rep.info_score != null ? Math.round(rep.info_score) : 0) + '/100</b>' +
      ' · 样本 ' + (rep.n_batches || 0) + '</div>';
    if (rep.satellite_warning) {
      html += '<div class="hint" style="color:#ffce4d">⚠ 检测到 satellite 式聚集：批量但设计空间未铺开，单纯加批增益有限——' +
        '优先铺开下方低覆盖因子/组合（论文 Set A→D：参数空间覆盖比样本数更重要）。</div>';
    }
    html += doeFactorRows(rep.factors);
    html += doeInteractionHtml(rep.interactions);
    html += doeRecHtml(rep.recommendations);
    if (rep.scope === "batch" && rep.transfer_readiness) {
      var tr = rep.transfer_readiness;
      var ready = tr.ready ? "✅ 可拟合迁移残差" : ("🟡 不足（" + tr.n_target_batches + "/" + tr.min_required + "）");
      html += '<div class="ksub" style="margin-top:8px">迁移就绪：' + ready + '</div>';
    }
    return html;
  }

  // 组合一个宿主的全部 DOE 区块（construct + repo / repo_by_scale）
  function doeHostBlock(host, showConstruct, scale) {
    var blk = (DATA.doe || {})[host];
    if (!blk) return '<div class="ksub">该宿主暂无 DOE 评估数据（构型库与批次库均无记录）。</div>';
    var html = '';
    if (showConstruct && blk.construct) {
      html += doeReportCard(blk.construct, "构型库设计空间覆盖（设计期）");
    }
    if (blk.repo_by_scale && scale && blk.repo_by_scale[scale]) {
      html += doeReportCard(blk.repo_by_scale[scale], "批次集 + 工艺设定点覆盖（" + scale + " 尺度）");
    } else if (blk.repo) {
      html += doeReportCard(blk.repo, "批次集 + 工艺设定点覆盖");
    }
    return html;
  }

  // DOE 区块（含可选宿主选择器）——Process 页 show_construct=false，Expression 页为 true
  function doeSectionHtml(opts) {
    opts = opts || {};
    var hosts = Object.keys(DATA.doe || {});
    var inner = '';
    if (opts.withHostSelect && hosts.length) {
      inner += '<div style="display:flex;gap:8px;align-items:center;margin-bottom:8px">' +
        '<span class="ksub">分析宿主</span>' +
        '<select id="doe_host_sel" onchange="FIP.doeSetHost(this.value)">' +
        hosts.map(function (h) {
          return '<option value="' + esc(h) + '"' + (h === opts.host ? " selected" : "") + '>' + esc(h) + '</option>';
        }).join("") + '</select></div>';
    }
    inner += '<div id="doe_inner">' + doeHostBlock(opts.host, opts.showConstruct, opts.scale) + '</div>';
    return '<div class="sub">🧪 实验设计信息量反馈 (DOE)</div>' + inner;
  }

  // ---- 首页 ----
  function renderHome(body) {
    var st = DATA.stats;
    var sr = st.n_batches ? Math.round(st.n_success / st.n_batches * 100) : 0;
    var kpis = [
      [st.n_batches, "历史批次", "成功 " + st.n_success + " / 失败 " + st.n_fail],
      [sr + "%", "放大成功率(代理)", "成功批次占比"],
      [st.n_time_points.toLocaleString(), "过程时序点", "TimescaleDB"],
      [st.n_constructs, "蛋白构架(PEFM)", "Expression Twin"],
      [st.n_knowledge, "知识库条目", "Scale-up DB"]
    ];
    var html = '<div class="kpirow">' + kpis.map(function (k) {
      return card('<div class="knum">' + k[0] + '</div><div class="klbl">' + k[1] + '</div><div class="ksub">' + k[2] + '</div>');
    }).join("") + "</div>";
    html += '<div style="display:flex;gap:14px;flex-wrap:wrap">';
    html += card('<div class="sub">宿主分布</div><div id="c_host"></div>', 'style="flex:1;min-width:320px"');
    html += card('<div class="sub">宿主 × 尺度</div><div id="c_heat"></div>', 'style="flex:1;min-width:320px"');
    html += "</div>";
    // 模型就绪度
    var thr = DATA.thresholds;
    var rows = [["Soft Sensor", "soft_sensor"], ["终点预测", "end_point"], ["Scale-up Twin", "scaleup_twin"], ["Expression Twin", "expression_twin"]].map(function (r) {
      var t = thr[r[1]];
      var act = { soft_sensor: Math.round(st.n_time_points / 1000), end_point: st.n_batches, scaleup_twin: Math.round(st.n_batches / 3), expression_twin: st.n_constructs }[r[1]];
      var status = act >= t.reliable ? "✅可靠" : (act >= t.start ? "🟡可用" : "🔴不足");
      return "<tr><td>" + r[0] + "</td><td>" + t.start + "</td><td>" + t.reliable + "</td><td>" + act + "</td><td>" + status + "</td></tr>";
    }).join("");
    html += card('<div class="sub">🧪 模型就绪度（对照《4.3 数据量门槛》）</div><table class="tbl"><tr><th>模型</th><th>启动门槛</th><th>较可靠门槛</th><th>当前样本</th><th>状态</th></tr>' + rows + "</table>");
    html += '<div class="hint">提示：左侧切换到 🧬 Expression Twin / ⚖️ Scale-up Twin / 📈 Process Twin 体验实时交互计算。</div>';
    body.innerHTML = html;
    var hostItems = Object.keys(st.hosts).map(function (h) { return { label: h.replace("_", " "), value: st.hosts[h], color: hostColor(h) }; });
    barChart($("c_host"), hostItems);
    var hosts = Object.keys(st.hosts);
    var scales = ["S", "M", "L"];
    var mat = hosts.map(function (h) {
      return scales.map(function (sc) {
        var c = 0; DATA.batches.forEach(function (b) { if (b.host_organism === h && b.scale_tag === sc) c++; }); return c;
      });
    });
    heatmap($("c_heat"), hosts.map(function (h) { return h.replace("_", " "); }), scales, mat);
  }

  // ---- 批次数据库 ----
  // 数据溯源 + 模型校准状态条：明确告知当前批次库的数据来源，以及哪些参数
  // 会随新数据自动更新、哪些必须在真实数据接入后人工重标定。
  function ptProvenanceBar() {
    var P = DATA.provenance;
    if (!P) return "";
    var synth = P.source !== "real";
    var bg = synth ? "#2a2410" : "#10241a", bd = synth ? "#aa3" : "#3a6";
    var auto = (P.auto_updated || []).map(function (s) { return "<li>" + esc(s) + "</li>"; }).join("");
    var need = (P.needs_recalibration || []).map(function (s) { return "<li>" + esc(s) + "</li>"; }).join("");
    var mc = P.softsensor_mode_counts || {};
    return '<div style="background:' + bg + ';border:1px solid ' + bd + ';border-radius:10px;padding:10px 12px;margin-bottom:10px">' +
      '<div style="font-weight:600;margin-bottom:4px">🔬 数据溯源与模型校准状态' +
      (synth ? ' · <span style="color:#ffce4d">当前为合成演示数据</span>' : ' · <span style="color:#3a6">真实生产数据</span>') + '</div>' +
      '<div class="hint" style="margin-bottom:6px">' + esc(P.generator || "") + '</div>' +
      '<div class="hint" style="margin-bottom:6px">' + esc(P.note || "") + '</div>' +
      '<div style="display:flex;gap:16px;flex-wrap:wrap">' +
      '<div style="flex:1;min-width:280px"><div class="ksub">✅ 随数据自动更新</div><ul style="margin:4px 0 0 18px;padding:0" class="hint">' + auto + '</ul></div>' +
      '<div style="flex:1;min-width:280px"><div class="ksub">⚠ 需人工/实验重标定</div><ul style="margin:4px 0 0 18px;padding:0" class="hint">' + need + '</ul></div>' +
      '</div>' +
      '<div class="hint" style="margin-top:6px">软测量模型 <b>' + esc(P.softsensor_model || "") + '</b> · 标定模式分布：' +
      'yield(得率标定) ' + (mc.yield || 0) + ' · scale(形状缩放) ' + (mc.scale || 0) + ' · prior(未标定) ' + (mc.prior || 0) + ' 批</div>' +
      '</div>';
  }
  // ---- 在线位号完备度面板 ----
  // 回答「批次库需要哪些在线参数原始数据」：按 Tier 列出模型实际消费的位号，
  // 并标出当前库里每个位号有几批真正提供了原始时序（缺失即意味着模型降级或用假设值）。
  function ptTagCoveragePanel() {
    var T = DATA.tag_tiers, S = DATA.tag_coverage_summary, C = DATA.tag_coverage;
    if (!T || !S) return "";
    var N = S.n_batch || 1;
    var tiers = [
      { key: "tier1_required", cn: "Tier-1 · 必需（缺失即降级）", color: "#e06c6c" },
      { key: "tier2_recommended", cn: "Tier-2 · 强烈推荐（缺失则用假设值）", color: "#e0a94d" },
      { key: "tier3_optional", cn: "Tier-3 · 可选（诊断增强）", color: "#5a9" }
    ];
    function row(t, spec, tierColor) {
      var n = (S.present || {})[t] || 0;
      var miss = (S.tier1_missing || {})[t] || 0;
      var ratio = n / N;
      var col = ratio >= 0.999 ? "#3a6" : (ratio > 0 ? "#e0a94d" : "#e06c6c");
      var bar = '<span style="display:inline-block;width:70px;height:6px;background:#333;border-radius:3px;vertical-align:middle;margin-right:6px">' +
        '<span style="display:inline-block;width:' + Math.round(ratio * 70) + 'px;height:6px;background:' + col + ';border-radius:3px"></span></span>';
      return '<tr>' +
        '<td style="color:' + tierColor + ';font-weight:600">' + esc(t) + '</td>' +
        '<td>' + esc(spec.cn) + '</td>' +
        '<td class="hint">' + esc(spec.unit) + '</td>' +
        '<td>' + bar + '<span class="hint">' + n + '/' + N + (miss ? ' <span style="color:#e06c6c">缺' + miss + '</span>' : '') + '</span></td>' +
        '<td class="hint">' + esc(spec.why) + '</td></tr>';
    }
    var blocks = tiers.map(function (tr) {
      var spec = T[tr.key] || {};
      var rows = Object.keys(spec).map(function (t) { return row(t, spec[t], tr.color); }).join("");
      return '<div style="margin-bottom:8px"><div class="ksub" style="color:' + tr.color + '">' + esc(tr.cn) + '</div>' +
        '<table class="tbl"><tr><th>位号</th><th>含义</th><th>单位</th><th>当前库覆盖</th><th>为什么需要</th></tr>' + rows + '</table></div>';
    }).join("");
    var o2m = S.o2_measured || 0;
    var warn = '<div class="hint" style="background:#2a2410;border:1px solid #aa3;border-radius:8px;padding:8px 10px;margin-top:6px">' +
      '⚠ <b>O₂ 流量实测批次 ' + o2m + '/' + N + '</b>：' +
      (o2m === 0
        ? '当前全库无独立纯氧位号，氧分压按空气 y<sub>O₂</sub>=20.95% 派生。<b>富氧/纯氧共通工艺无法核算</b>，Scale-up 的「降空气底通 + 提氧占比」方案缺少数据支撑。接入真实数据后请务必上传 <code>O2</code> 位号。'
        : '已有 ' + o2m + ' 批提供纯氧流量实测；其余批次氧分压仍按空气派生。') +
      '　罐压（PRESS）缺失时 C* 按名义压力取，OTR 放大精度下降；罐重（WT）缺失时工作体积按名义值。' +
      '</div>';
    return card('<div class="sub">🔌 在线参数位号完备度（模型实际消费的原始时序）</div>' +
      '<div class="hint" style="margin-bottom:6px">下表位号即为「批次数据库需要哪些在线原始数据」的答案：' +
      'Tier-1 缺任一会让对应模块直接降级或无法计算，Tier-2 缺失则用假设值（已在界面标出），Tier-3 仅增强诊断。' +
      '导入真实数据时列名应对齐这些位号（见「下载模板」）。</div>' + blocks + warn);
  }
  function renderBatches(body) {
    var hosts = ["", "E.coli_BL21(DE3)", "Pichia_pastoris_X33"];
    var scales = ["", "S", "M", "L"];
    var succ = ["", "true", "false"];
    var rows = DATA.batches.filter(function (b) {
      if (state.filter.host && b.host_organism !== state.filter.host) return false;
      if (state.filter.scale && b.scale_tag !== state.filter.scale) return false;
      if (state.filter.success && String(b.success) !== state.filter.success) return false;
      return true;
    });
    var html = card(
      '<div class="sub">数据接口</div><div style="display:flex;gap:10px;flex-wrap:wrap;align-items:center">' +
      '<button class="btn2" onclick="FIP.exportCSV()">⬇ 导出 CSV</button>' +
      '<button class="btn2" onclick="FIP.exportJSON()">⬇ 导出 JSON</button>' +
      '<button class="btn2" onclick="FIP.downloadTemplate()">⬇ 下载模板</button>' +
      '<button class="btn2" onclick="FIP.downloadTsTemplate()">⬇ 下载时序模板（在线位号）</button>' +
      '<label class="btn2" style="cursor:pointer;margin:0">⬆ 导入文件<input id="f_import" type="file" accept=".csv,.json" style="display:none" onchange="FIP.importFile(this.files[0])"></label>' +
      '<span id="io_msg" class="ksub"></span></div>'
    );
    html += card(
      '<div class="sub">筛选</div><div style="display:flex;gap:10px;flex-wrap:wrap">' +
      '<select id="f_host" onchange="FIP.set(\'filter.host\',this.value)">' + hosts.map(function (h) { return '<option value="' + h + '"' + (state.filter.host === h ? " selected" : "") + '>' + (h || "全部宿主") + '</option>'; }).join("") + '</select>' +
      '<select id="f_scale" onchange="FIP.set(\'filter.scale\',this.value)">' + scales.map(function (s) { return '<option value="' + s + '"' + (state.filter.scale === s ? " selected" : "") + '>' + (s || "全部尺度") + '</option>'; }).join("") + '</select>' +
      '<select id="f_succ" onchange="FIP.set(\'filter.success\',this.value)">' + succ.map(function (s) { return '<option value="' + s + '"' + (state.filter.success === s ? " selected" : "") + '>' + ({ "": "全部结果", "true": "成功", "false": "失败" }[s]) + '</option>'; }).join("") + '</select>' +
      '<span class="ksub">共 ' + rows.length + ' 批</span></div>'
    );
    html += card('<div class="sub">批次列表（点击查看工艺曲线）</div><div style="max-height:420px;overflow:auto"><table class="tbl" id="btab">' +
      "<tr><th>批次</th><th>宿主</th><th>尺度</th><th>体积(L)</th><th>时长(h)</th><th>OD</th><th>效价(g/L)</th><th>结果</th></tr>" +
      rows.map(function (b) {
        return '<tr class="brow" onclick="FIP.viewBatch(\'' + b.batch_id + '\')"><td>' + b.batch_id + '</td><td style="color:' + hostColor(b.host_organism) + '">' + b.host_organism.replace("_", " ") + '</td><td>' + b.scale_tag + '</td><td>' + fmt(b.working_volume_l, 1) + '</td><td>' + fmt(b.duration_h, 1) + '</td><td>' + fmt(b.harvest_od, 1) + '</td><td>' + fmt(b.titer_g_l, 2) + '</td><td>' + (b.success ? "✅" : "❌") + '</td></tr>';
      }).join("") + "</table></div>");
    html += card('<div class="sub">在线位号完备度</div>' +
      '<button class="btn2" onclick="FIP.toggleTags()">' + (state.showTags ? '▾ 收起' : '▸ 展开') + ' Tier-1/2/3 位号清单与当前库覆盖率</button>' +
      '<div id="tag_panel">' + (state.showTags ? ptTagCoveragePanel() : '') + '</div>');
    html += '<div id="batch_detail"></div>';
    body.innerHTML = ptProvenanceBar() + html;
  }
  function toggleTags() {
    state.showTags = !state.showTags;
    renderBatches($("body-content"));
  }

  function viewBatch(bid) {
    state.selBatch = bid;
    var b = DATA.batches.filter(function (x) { return x.batch_id === bid; })[0];
    var ts = DATA.timeseries[bid];
    var cv = (DATA.tag_coverage || {})[bid];
    var cvHtml = "";
    if (cv) {
      var miss = cv.tier1_missing || [];
      var col = miss.length ? "#e06c6c" : "#3a6";
      cvHtml = '<div class="hint" style="margin-top:4px">🔌 在线位号：已接入 ' + (cv.present || []).join(" / ") +
        '　·　<span style="color:' + col + '">Tier-1 缺失：' + (miss.length ? miss.join(" / ") : "无") + '</span>' +
        '　·　O₂ 来源：<b>' + (cv.o2_source === "measured" ? "实测纯氧位号" : "由空气×20.95% 派生") + '</b></div>';
    }
    var html = card('<div class="sub">📈 ' + bid + ' · ' + (b ? b.host_organism.replace("_", " ") : "") + ' · ' + (b ? b.scale_tag : "") + '</div>' +
      (b ? '<div class="ksub">体积 ' + fmt(b.working_volume_l, 1) + ' L · 时长 ' + fmt(b.duration_h, 1) + ' h · 效价 ' + fmt(b.titer_g_l, 2) + ' g/L · ' + (b.success ? "成功" : "失败" + (b.root_cause ? "（" + esc(b.root_cause) + "）" : "")) + '</div>' : "") +
      cvHtml +
      '<div id="d_lines"></div><div id="d_biomass"></div>');
    var d = $("batch_detail"); if (!d) { renderBatches($("body-content")); d = $("batch_detail"); }
    d.innerHTML = html;
    if (ts) {
      var hrs = ts.hours;
      lineChart($("d_lines"), [
        { name: "DO", color: PRIMARY, points: hrs.map(function (h, i) { return [h, ts.DO[i]]; }), unit: "%" },
        { name: "pH", color: "#ffce4d", points: hrs.map(function (h, i) { return [h, ts.pH[i]]; }), unit: "" },
        { name: "OUR", color: ACCENT, points: hrs.map(function (h, i) { return [h, ts.OUR[i]]; }), unit: "mmol/L·h" },
        { name: "CER", color: "#b06ab3", points: hrs.map(function (h, i) { return [h, ts.CER[i]]; }), unit: "mmol/L·h" }
      ], { ytitle: "过程变量" });
      lineChart($("d_biomass"), [
        { name: "DCW(软测量)", color: PRIMARY, points: hrs.map(function (h, i) { return [h, ts.DCW[i]]; }), unit: "g/L" }
      ], { ymin: 0, ytitle: "DCW (g/L)" });
    }
  }

  // -------------------------------------------------------------------------
  // Process Twin 深化计算（前端演示级数字孪生）
  // -------------------------------------------------------------------------
  // ---- 软测量：氧/碳衡算 ODE 反演生物量（与后端 soft_sensor.dcw_from_our_balance 同构）----
  // 生理方程（≡ 机理生成器）：OUR(t) = X·(qO2_g·μ + mO2)，μ = (1/X)·dX/dt
  //   ⇒ dX/dt = -(mO2/qO2_g)·X + OUR(t)/qO2_g        ← 关于 X 的线性 ODE
  // 对分段线性 OUR 解析积分（指数积分），无需假设比摄氧率 OTA 恒定。
  // 量纲：Y_OX = 1/qO2_g [g DCW / mmol O₂]，标定后典型 0.04–0.085。
  var SS_DEFAULT_Y_OX = 1 / 12;   // g DCW/mmol O₂（E. coli 先验 ≡ qO2_g=12）
  var SS_DEFAULT_MO2 = 1.5;       // mmol O₂/gX/h 维持
  var SS_DEFAULT_X0 = 0.3;        // g/L 接种生物量
  var SS_RQ = 1.05;               // CER/OUR（葡萄糖呼吸）
  // 生物学可达上限（g DCW/L）：超过即为情景外产物，需钳制并提示
  var HOST_XMAX = { ecoli: 200, yeast: 400, other: 250 };

  // 宿主生理常数（先验；真实数据接入后应由跨批次回归重标定，见 SS_PROVENANCE）
  var SS_HOST_CONST = {
    ecoli: { qo2_g: 12, mo2: 1.5, x0: 0.3, dryFrac: 0.22 },
    yeast: { qo2_g: 10, mo2: 1.0, x0: 0.5, dryFrac: 0.26 },
    other: { qo2_g: 12, mo2: 1.5, x0: 0.3, dryFrac: 0.24 }
  };
  function ptHostConst(host) { return SS_HOST_CONST[host] || SS_HOST_CONST.other; }
  function ptSoftMeta(bid) {
    var m = (DATA.softsensor && DATA.softsensor[bid]) || null;
    if (!m) {
      var h = ptHostOf(bid); if (h === "other") h = "ecoli";
      var c = ptHostConst(h);
      return { host: h, qo2_g: c.qo2_g, y_ox: 1 / c.qo2_g, mo2: c.mo2, x0: c.x0,
        calibrated: false, calib: 1, raw_end: 0, anchor: null, dry_frac: c.dryFrac, mode: "prior" };
    }
    return m;
  }
  // 终点锚定：二分标定 qo2_g 使 X(t_end) = 离线 harvest DCW（与后端 dcw_from_our_balance 同构）
  function ptCalibrateYield(hours, our, host, x0, mo2, anchor) {
    var C = ptHostConst(host), qg0 = C.qo2_g, lo = 0.5, hi = 400;
    if (!(anchor > 0)) return { qo2_g: qg0, y_ox: 1 / qg0, mode: "prior", diag: null };
    for (var i = 0; i < 60; i++) {
      var mid = (lo + hi) / 2;
      var xe = ptDcwFromGas(hours, our, 1 / mid, mo2, x0);
      if (xe[xe.length - 1] > anchor) lo = mid; else hi = mid;
    }
    var qg = (lo + hi) / 2, xe2 = ptDcwFromGas(hours, our, 1 / qg, mo2, x0), xeE = xe2[xe2.length - 1];
    if (Math.abs(xeE - anchor) <= 0.02 * anchor && qg >= 3 && qg <= 120)
      return { qo2_g: qg, y_ox: 1 / qg, mode: "yield", diag: null };
    // 生理窗口外 / 氧耗能量不平衡：保持 ODE 形状，按整体缩放至实测终点
    var raw = ptDcwFromGas(hours, our, 1 / qg0, mo2, x0), rawE = raw[raw.length - 1];
    var den = Math.max(rawE - x0, 1e-6), k = (anchor - x0) / den;
    return { qo2_g: qg0, y_ox: 1 / qg0, mode: "scale", scaleK: k,
      diag: Math.abs(xeE - anchor) > 0.02 * anchor
        ? "累计耗氧无法解释该生物量（可达上界 " + fmt(xeE, 1) + " g/L < 实测 " + fmt(anchor, 1) + " g/L），疑尾气分析仪漂移或采样丢峰"
        : "标定常数 qo2_g=" + fmt(qg, 1) + " mmol/gX 超出 3–120 生理窗口，疑维持/非生长耗氧占比异常" };
  }
  // 由尾气序列重建软测量列（真实数据通常不含 DCW/WCW，导入后需即时补算）
  function ptRebuildSoftSensor(bid, ts) {
    if (!ts || !ts.hours || !ts.OUR || !ts.hours.length) return null;
    var b = null; DATA.batches.forEach(function (x) { if (x.batch_id === bid) b = x; });
    var host = ptHostType(b || {}); if (host === "other") host = "ecoli";
    var C = ptHostConst(host);
    var anchor = (b && b.harvest_dcw_g_l > 0) ? b.harvest_dcw_g_l : null;
    var cal = ptCalibrateYield(ts.hours, ts.OUR, host, C.x0, C.mo2, anchor);
    var X = ptDcwFromGas(ts.hours, ts.OUR, cal.y_ox, C.mo2, C.x0);
    if (cal.mode === "scale") X = X.map(function (v) { return Math.max(0, C.x0 + (v - C.x0) * cal.scaleK); });
    ts.DCW = X.map(function (v) { return Math.round(v * 1000) / 1000; });
    ts.WCW = X.map(function (v) { return Math.round(v / C.dryFrac * 1000) / 1000; });
    if (!ts.O2_FLOW && ts.AIR) ts.O2_FLOW = ts.AIR.map(function (v) { return Math.round(v * 0.2095 * 1e4) / 1e4; });
    if (!ts.phases && ts.FEED) ts.phases = [{ label: "发酵全程", start: ts.hours[0], end: ts.hours[ts.hours.length - 1], color: "#3b6ea5" }];
    if (!DATA.softsensor) DATA.softsensor = {};
    DATA.softsensor[bid] = {
      host: host, qo2_g: Math.round(cal.qo2_g * 100) / 100, y_ox: Math.round(cal.y_ox * 1e4) / 1e4,
      mo2: C.mo2, x0: C.x0, raw_end: Math.round(X[X.length - 1] * 100) / 100,
      anchor: anchor, calibrated: !!anchor, mode: cal.mode,
      calib: anchor ? Math.round(anchor / Math.max(X[X.length - 1], 1e-9) * 1000) / 1000 : 1,
      diag: cal.diag, dry_frac: C.dryFrac, source: "frontend-rebuilt"
    };
    return DATA.softsensor[bid];
  }
  function ptDcwFromGas(hours, arr, yOx, mo2, x0) {
    var qg = 1 / Math.max(yOx, 1e-6), k = mo2 / qg, out = [], x = x0;
    if (!hours.length) return out;
    out.push(x0);
    for (var i = 0; i < hours.length - 1; i++) {
      var dt = hours[i + 1] - hours[i];
      if (!(dt > 0)) { out.push(x); continue; }
      var s = (arr[i + 1] - arr[i]) / dt, kd = k * dt, e = Math.exp(-kd);
      var phi1 = kd < 1e-6 ? dt : (1 - e) / k;
      var phi2 = kd < 1e-6 ? dt * dt / 2 : (dt - (1 - e) / k) / k;
      x = e * x + (arr[i] * phi1 + s * phi2) / qg;
      out.push(x < 0 ? 0 : x);
    }
    return out;
  }
  function ptTrapz(xs, ys) {
    var s = 0;
    for (var i = 1; i < Math.min(xs.length, ys.length); i++) {
      var dx = xs[i] - xs[i - 1];
      if (dx > 0) s += (ys[i] + ys[i - 1]) / 2 * dx;
    }
    return s;   // ∫y dx
  }
  // 【第三法】逻辑斯蒂（Verhulst）生长律拟合 —— 不依赖尾气信号，仅检验
  // 生物量曲线是否符合典型 S 型生长：X(t)=K/(1+A·e^{-μt})，A=(K−X₀)/X₀。
  // 替代旧的「全周期指数外推」：后者把补料期的 μ 衰减段强行指数化，会系统性高估终点 1.5× 以上。
  function ptDcwLogistic(hours, Xobs) {
    var n = hours.length;
    if (n < 3) return Xobs.slice();
    var x0 = Math.max(1e-3, Xobs[0]), xe = Xobs[n - 1];
    if (!(xe > 0)) return Xobs.slice();
    var bestK = xe * 1.2, bestMu = 0.15, bestErr = Infinity;
    for (var ki = 0; ki < 60; ki++) {
      var K = xe * (0.95 + 0.06 * ki);                       // 0.95 – 4.5 ×终点
      var A = (K - x0) / Math.max(x0, 1e-6);
      for (var mi = 0; mi < 40; mi++) {
        var mu = 0.02 + 0.012 * mi;                          // 0.02 – 0.50 1/h
        var err = 0;
        for (var i = 0; i < n; i++) {
          var pred = K / (1 + A * Math.exp(-mu * hours[i]));
          var d = Math.log(Math.max(pred, 1e-6)) - Math.log(Math.max(Xobs[i], 1e-6));
          err += d * d;
        }
        if (err < bestErr) { bestErr = err; bestK = K; bestMu = mu; }
      }
    }
    var A2 = (bestK - x0) / Math.max(x0, 1e-6);
    return hours.map(function (h) { return bestK / (1 + A2 * Math.exp(-bestMu * h)); });
  }
  function ptKPI(ts) {
    var hrs = ts.hours, OUR = ts.OUR, CER = ts.CER, DCW = ts.DCW, rq = [], mu = [], ota = [];
    for (var i = 0; i < hrs.length; i++) {
      rq.push(CER[i] ? OUR[i] / CER[i] : 0);
      if (i === 0) mu.push(0);
      else { var dt = hrs[i] - hrs[i - 1]; mu.push(dt > 0 ? (Math.log(Math.max(DCW[i], 0.01)) - Math.log(Math.max(DCW[i - 1], 0.01))) / dt : 0); }
      ota.push(DCW[i] ? OUR[i] / DCW[i] : 0);
    }
    var rqAvg = rq.reduce(function (a, b) { return a + b; }, 0) / Math.max(1, rq.length);
    var otaEnd = ota.length ? ota[ota.length - 1] : 0;
    var muEnd = mu.length ? mu[mu.length - 1] : 0;
    return { rq: rq, mu: mu, ota: ota, rqAvg: rqAvg, otaEnd: otaEnd, muEnd: muEnd };
  }
  function ptGolden(curB, curTs, curKpi, refBid) {
    var ref = null, rts = null;
    DATA.batches.forEach(function (x) { if (x.batch_id === refBid) ref = x; });
    if (DATA.timeseries[refBid]) rts = DATA.timeseries[refBid];
    if (!ref || !rts) return null;
    var rkpi = ptKPI(rts);
    var curDCWe = curTs.DCW[curTs.DCW.length - 1], refDCWe = rts.DCW[rts.DCW.length - 1];
    var curRQe = curTs.OUR[curTs.OUR.length - 1] / curTs.CER[curTs.CER.length - 1];
    var refRQe = rts.OUR[rts.OUR.length - 1] / rts.CER[rts.CER.length - 1];
    var rows = [
      ["终点 OD", curB.harvest_od, ref.harvest_od],
      ["终点 DCW (g/L)", curDCWe, refDCWe],
      ["效价 (g/L)", curB.titer_g_l, ref.titer_g_l],
      ["RQ(末)", curRQe, refRQe],
      ["μ(末) 1/h", curKpi.muEnd, rkpi.muEnd]
    ];
    return rows.map(function (r) {
      var dev = r[2] ? (r[1] - r[2]) / r[2] * 100 : 0;
      var lvl = Math.abs(dev) > 20 ? "bad" : (Math.abs(dev) > 10 ? "warn" : "ok");
      return { name: r[0], cur: r[1], ref: r[2], dev: dev, lvl: lvl };
    });
  }
  // What-if 情景因子：把 DO/温度/补料的影响拆成「生物量响应 modX」与「比产物生成速率响应 modQ」。
  // 效价 P = ∫ q_P·X dt ⇒ 相对变化 ≈ modX·modQ（**一次**，旧实现把 DCW 与 titer 同时乘同一因子 => mod² 错误放大）。
  function ptScenarioFactors(p, base, bid) {
    // base 可能为空（未选批次 / 首次进入页面时拖动滑块）：旧实现直接读 base.doRef，
    // 此时 base 为 null 会抛 TypeError 并使整个情景面板白屏，故先兜底为空对象。
    if (!base) base = {};
    var doRef = Math.max(1, base.doRef || 30);          // 批次实际工况 DO（中位）
    var doSet = Math.max(1, p.do || doRef);
    var T = Math.max(15, Math.min(45, p.temp || 30));
    var F = Math.max(0.1, Math.min(3, p.feed || 1));
    // DO：以**批次实际工况 DO** 为参考，不足→氧限制溢流；>60% 属过高溶氧（ROS/氧化应激）才惩罚
    var doX = 0.75 + 0.25 * Math.min(1, doSet / doRef);
    if (doSet > 60) doX *= 1 - Math.min(0.18, 0.05 * (doSet - 60) / 10);
    var doQ = 0.62 + 0.38 * Math.min(1.35, doSet / doRef);
    if (doSet > 60) doQ *= 1 - Math.min(0.12, 0.035 * (doSet - 60) / 10);
    // 温度：**相对批次实际温度**的高斯响应（T==参考温度 ⇒ 因子 1，滑块回基线时情景=基线）。
    //   生长最适偏 34–37 ℃，重组表达/可溶性偏 30 ℃ —— 两条不同的响应曲线。
    var tRef = Math.max(15, Math.min(45, (base && base.tempRef) || 30));
    function relGauss(t, opt, sg) { return Math.exp(-(Math.pow(t - opt, 2) - Math.pow(tRef - opt, 2)) / (2 * sg * sg)); }
    var tX = relGauss(T, 34, Math.sqrt(20));
    var tQ = relGauss(T, 30, Math.sqrt(16));
    // 补料：X 随补料上升但收益递减；补料过量→溢流/乙酸，比产率略降
    var fX = Math.pow(F, 0.65), fQ = Math.pow(F, -0.15);
    function cl(v) { return Math.max(0.25, Math.min(2.2, v)); }
    return {
      modX: cl(doX * tX * fX), modQ: cl(doQ * tQ * fQ),
      do: doSet, temp: T, feed: F, doRef: doRef,
      xMax: (bid && DATA.batches) ? ptXmax(ptHostOf(bid)) : 200
    };
  }
  function ptHostOf(bid) {
    var b = null; DATA.batches.forEach(function (x) { if (x.batch_id === bid) b = x; });
    return ptHostType(b || {});
  }
  function ptXmax(host) { return HOST_XMAX[host] || HOST_XMAX.other; }
  // 批次实际工况参考值：DO/TEMP 取**诱导期后的中位数**。
  // 不能用末值：DO 末值常因氧耗崩溃趋于 0（如取全局末值会令 referencias 因子失真）。
  function ptMedian(arr, floor) {
    if (!arr || !arr.length) return 0;
    var v = arr.slice().sort(function (a, b) { return a - b; });
    if (floor != null) { v = v.filter(function (x) { return x > floor; }); }
    if (!v.length) return arr[arr.length - 1];
    var m = Math.floor(v.length / 2);
    return v.length % 2 ? v[m] : (v[m - 1] + v[m]) / 2;
  }
  function ptDoRef(ts) { return ts && ts.DO && ts.DO.length ? ptMedian(ts.DO, 1) : 30; }
  function ptTempRef(ts) { return ts && ts.TEMP && ts.TEMP.length ? ptMedian(ts.TEMP) : 30; }
  function ptWhatIf(p, base, bid) {
    var fa = ptScenarioFactors(p, base, bid);
    var titer = (base.titer || 0) * fa.modX * fa.modQ;
    var dur = (base.dur || 0) * (1 / Math.max(0.35, Math.pow(fa.feed, 0.5))) *
      (fa.temp < 30 ? 1 + (30 - fa.temp) / 10 * 0.3 : 1);
    return { titer: titer, dur: dur, fa: fa, rel: base.titer ? (titer / base.titer - 1) * 100 : 0 };
  }
  function ptSliderRow(label, path, val, min, max, step) {
    var idv = "sl_" + path.replace(/\./g, "_");
    return '<div class="slider-row"><label>' + label + ' <span class="sv" id="' + idv + '">' + fmt(val, 2) + '</span></label>' +
      '<input type="range" min="' + min + '" max="' + max + '" step="' + step + '" value="' + val + '" onchange="FIP.set(\'' + path + '\',parseFloat(this.value))"></div>';
  }
  // What-if 发酵时长：进度条 + 手填数字，二者绑定到 process.dur
  function ptDurRow() {
    var v = state.process.dur;
    return '<div class="slider-row"><label>发酵时长 (h) <span class="sv" id="sl_process_dur">' + fmt(v, 0) + '</span></label>' +
      '<input type="range" min="12" max="160" step="1" value="' + v + '" onchange="FIP.set(\'process.dur\',parseFloat(this.value))">' +
      '<input type="number" min="12" max="160" step="1" value="' + v + '" style="width:72px" onchange="FIP.set(\'process.dur\',parseFloat(this.value))"></div>';
  }
  function ptKpiCards(kpi, b, ts, bid) {
    // 累计耗氧必须是 ∫OUR dt（对时间的梯形积分，mmol O₂/L）；旧实现为简单求和（未乘 Δt），
    // 既非正确量纲也与采样密度耦合，导致碳源效率失真。
    var hours = ts ? ts.hours : [];
    var cumO = (ts && hours.length) ? ptTrapz(hours, ts.OUR) : 0;                 // mmol O₂/L
    var cumCO2 = (ts && hours.length) ? ptTrapz(hours, ts.CER) : 0;               // mmol CO₂/L
    var AUC_X = (ts && hours.length && ts.DCW) ? ptTrapz(hours, ts.DCW) : 0;      // g DCW·h/L
    // 体积产率优先取实测 outcome 字段（= titer / 实际培养时长）
    var volRate = (b && b.productivity_g_l_h != null && !isNaN(b.productivity_g_l_h))
      ? b.productivity_g_l_h
      : (b && ts ? b.titer_g_l / Math.max(0.1, b.duration_h) : 0);
    var carbonEff = (b && cumO > 0) ? b.titer_g_l / cumO : 0;                     // g 产物 / mmol O₂
    // 得率交叉校核：实测 harvest DCW 与累计耗氧的一致性（净得率 = ΔX / ∫OUR dt）
    var SMk = ptSoftMeta(bid);
    var endX = (ts && ts.DCW && ts.DCW.length) ? ts.DCW[ts.DCW.length - 1] : 0;
    var yOxObs = (cumO > 0 && endX > 0) ? (endX - (SMk.x0 || 0)) / cumO : 0;      // 表观净得率
    // —— 累计耗氧是"总量"指标，读法必须带场景：批次 + 时长 + 时间均值 + 峰值，否则会被误读为"越大越好"。
    // 本 KPI 组服务于"看实测批次"：场景标签直接带上 batch_id，避免与 What-if 预测混淆。
    var spanH = (hours.length > 1) ? (hours[hours.length - 1] - hours[0]) : (b && b.duration_h ? b.duration_h : 0);
    var ourMean = (spanH > 0) ? cumO / spanH : 0;                                 // mmol O₂/L/h（时间均值）
    var carbonEffMol = carbonEff * 1000;                                          // g 产物 / mol O₂
    // 峰值 OUR：用 3 点滑动中位数，抑制"孤立单点抖动"（某一帧尾气/DO 换算异常）。
    // 注意它**不抑制趋势**：实测批次的末段上升（如 …214,240,263,281,309,337）是
    // 单调趋势而非离群，中位数会原样保留，此时的峰值就是末点，属真实读数，
    // 因此下面对"末段趋势"与"孤立尖峰"分别标注，不再一律写成"去抖"。
    var ourPk = { val: 0, t: 0, rawMax: 0, rawMaxT: 0, n: 0 };
    if (ts && ts.OUR && ts.OUR.length) {
      var _arr = ts.OUR;
      for (var _k = 0; _k < _arr.length; _k++) {
        var _vk = _arr[_k];
        if (isFinite(_vk) && _vk > ourPk.rawMax) { ourPk.rawMax = _vk; ourPk.rawMaxT = hours[_k] || 0; }
      }
      for (var _i = 0; _i < _arr.length; _i++) {
        var _v = _arr[_i];
        if (!isFinite(_v)) continue;
        ourPk.n++;
        var _a = isFinite(_arr[_i - 1]) ? _arr[_i - 1] : _v;
        var _c = isFinite(_arr[_i + 1]) ? _arr[_i + 1] : _v;
        var _med = Math.max(Math.min(_a, _v), Math.min(Math.max(_a, _v), _c));  // 三点中位数
        if (_med > ourPk.val) { ourPk.val = _med; ourPk.t = hours[_i] || 0; }
      }
    }
    // 峰/均 &gt; 2 视为存在瞬时供氧尖峰；样本过少（&lt;5）则统计不可信，不做尖峰判定
    var peakRatio = (ourMean > 0) ? ourPk.val / ourMean : 0;
    var spikeFlag = (ourPk.n >= 5 && peakRatio >= 2.0);
    var spikeTxt = spikeFlag ? "（峰/均 " + fmt(peakRatio, 2) + "，存在瞬时供氧尖峰）" : "";
    // 区分「末段单调趋势」与「孤立尖峰」：前者是真实读数（中位数抑制不了），
    // 后者才是需要排除的异常帧。
    var lastIdx = (ts && ts.OUR) ? ts.OUR.length - 1 : 0;
    var _tail = [];
    if (ts && ts.OUR && ts.OUR.length >= 3) {
      for (var _q = Math.max(0, ts.OUR.length - 3); _q < ts.OUR.length; _q++) {
        if (isFinite(ts.OUR[_q])) _tail.push(ts.OUR[_q]);
      }
    }
    var tailRising = (_tail.length >= 3) &&
      (_tail[_tail.length - 1] > _tail[_tail.length - 2] && _tail[_tail.length - 2] > _tail[_tail.length - 3]);
    var tailOutlier = (ts && ts.OUR && ts.OUR.length >= 5 && ourPk.val > 0 &&
      !tailRising && ourPk.rawMax > ourPk.val * 1.5);
    var scenario = (b && b.batch_id ? b.batch_id : (bid || "—")) + " · " + (b && b.host_organism ? String(b.host_organism).replace("_", " ") : "—") +
      (b && b.scale_tag ? " · " + b.scale_tag : "") + (spanH > 0 ? " · " + fmt(spanH, 0) + " h" : "");
    var cards = [
      { l: "RQ（OUR/CER）均值", v: fmt(kpi.rqAvg, 2), s: "呼吸商，≈1 为混合代谢" },
      { l: "比生长速率 μ(末) 1/h", v: fmt(kpi.muEnd, 3), s: "对数期外推" },
      { l: "OTA 比摄氧率(末)", v: fmt(kpi.otaEnd, 2), s: "OUR/DCW，mmol O₂/g DCW/h" },
      { l: "峰值 OUR", v: fmt(ourPk.val, 0), s: "mmol O₂/L/h · t=" + fmt(ourPk.t, 1) + " h（3 点滑动中位数，仅抑制孤立单点抖动）" + spikeTxt +
          (tailRising && ourPk.rawMaxT > 0 ? "<br>ℹ 峰值落在末段（末 3 点单调上升、属趋势末值而非离群，中位数不会抑制趋势）" : "") +
          (tailOutlier ? "<br>⚠ 原始最大 OUR=" + fmt(ourPk.rawMax, 0) + " 疑为孤立异常帧，已排除出峰值" : "") },
      { l: "OUR 时间均值", v: fmt(ourMean, 1), s: "∫OUR dt ÷ 时长，mmol O₂/L/h" },
      { l: "累计耗氧 ∫OUR dt", v: fmt(cumO, 0), s: "mmol O₂/L · " + (spanH > 0 ? "全程 " + fmt(spanH, 0) + " h" : "时长未知") +
          " 积分<br><b>总量指标：随培养时长线性增长，不可跨批次直接比大小</b>", hl: true },
      { l: "碳效率 g 产物/mol O₂", v: fmt(carbonEffMol, 2), s: "titer ÷ 累计耗氧<br>可跨批次比较的氧利用经济性" },
      { l: "单位时间产率 g/L/h", v: fmt(volRate, 3), s: "titer ÷ 实际培养时长" },
      { l: "表观得率 Y_X/O₂", v: fmt(yOxObs, 4), s: "（终点DCW−接种）/∫OUR dt，g/mmol" }
    ];
    return '<div class="kpi-grid">' + cards.map(function (c) {
      return '<div class="kpi-card"' + (c.hl ? ' style="border-color:#1f9e8955"' : "") + '><div class="kpi-val">' + c.v + '</div><div class="kpi-lab">' + c.l + '</div><div class="kpi-sub">' + c.s + '</div></div>';
    }).join("") + '</div>' +
      '<div class="method-note" style="margin-top:8px">🧪 <b>实测批次</b>：' + esc(scenario) + '（本组为历史批次实测指标，非 What-if 预测）</div>' +
      '<div class="method-note" style="margin-top:6px">📖 <b>怎么读累计耗氧</b>：它是"这一批总共烧了多少氧"的总量指标，同工艺下培养时间越长数值越大，<b>不是越大越好</b>。' +
      '判断供氧是否成为瓶颈看<b>峰值 OUR</b>（3 点滑动中位数，只抑制孤立单点抖动；若末段是单调上升趋势，峰值即末点，属真实读数、不会被抑制；峰/均 ≥ 2.0 提示存在瞬时供氧尖峰，越高越接近设备传质上限）；' +
      '判断氧花得值不值看<b>碳效率</b>（每消耗 1 mol O₂ 得到多少 g 产物，可跨批次比较）；' +
      '两者结合<b>累计耗氧</b>才能区分"长周期高产"与"长周期低效"。</div>';
  }
  function ptDevTable(dev) {
    var body = dev.map(function (r) {
      var cls = { ok: "dev-ok", warn: "dev-warn", bad: "dev-bad" }[r.lvl];
      return '<tr><td>' + r.name + '</td><td>' + fmt(r.cur, 2) + '</td><td>' + fmt(r.ref, 2) + '</td><td class="' + cls + '">' + (r.dev >= 0 ? "+" : "") + fmt(r.dev, 0) + '%</td></tr>';
    }).join("");
    return '<table class="dev-table"><tr><th>指标</th><th>当前批</th><th>标杆批</th><th>偏离</th></tr>' + body + '</table>' +
      '<div class="method-note">偏离 &gt;10% 提示、&gt;20% 预警（红色）；用于快速定位异常批次。</div>';
  }
  function ptPredBox(wi) {
    var risk = wi.rel < -15;
    return '<div class="predbox' + (risk ? " risk" : "") + '"><div class="pv">' + fmt(wi.titer, 2) + ' <small>g/L</small></div>' +
      '<div class="pl">预测终点效价（基线 ' + (wi.rel >= 0 ? "+" : "") + fmt(wi.rel, 0) + '%）</div>' +
      '<div class="pv" style="font-size:18px;margin-top:10px">' + fmt(wi.dur, 0) + ' <small>h</small></div><div class="pl">预测培养时长</div></div>';
  }

  // ---- Process Twin ----
  // 在线曲线可选参数（顺序对应需求：DO/pH/温度/转速/空气通量/氧气通量/OUR/CER/WCW/DCW/补料速率）
  var PT_PARAMS = [
    { key: "DO", label: "DO", unit: "%", color: PRIMARY },
    { key: "pH", label: "pH", unit: "", color: "#ffce4d" },
    { key: "TEMP", label: "温度", unit: "℃", color: "#ff8c42" },
    { key: "RPM", label: "转速", unit: "rpm", color: "#4cc9f0" },
    { key: "AIR", label: "空气通量", unit: "L/min", color: "#7bdff2" },
    { key: "O2_FLOW", label: "氧气通量", unit: "L/min", color: "#5b8def" },
    { key: "OUR", label: "OUR", unit: "mmol/L·h", color: ACCENT },
    { key: "CER", label: "CER", unit: "mmol/L·h", color: "#b06ab3" },
    { key: "WCW", label: "WCW", unit: "g/L", color: "#80ed99" },
    { key: "DCW", label: "DCW", unit: "g/L", color: "#a0e8af" },
    { key: "FEED", label: "补料速率", unit: "相对", color: "#ff6b6b" }
  ];
  function ptHostType(b) {
    if (b.host_type) return b.host_type;
    var ho = b.host_organism || "";
    if (/Pichia|pastoris|GS115|X33|KM71|Mut[sS]?/.test(ho)) return "yeast";   // 仅毕赤酵母归为酵母模型（不含 S. cerevisiae）
    if (/E\.coli|coli/.test(ho)) return "ecoli";
    return "other";
  }
  function ptParamChips() {
    return '<div class="sub" style="margin-top:6px">参数可选模块（点击显隐，各参数独立纵坐标）</div><div style="display:flex;gap:8px;flex-wrap:wrap">' +
      PT_PARAMS.map(function (p) {
        var on = state.process.params[p.key];
        return '<span class="chip' + (on ? " chip-on" : "") + '" onclick="FIP.set(\'process.params.' + p.key + '\', ' + (!on) + ')"><span class="chip-dot" style="background:' + p.color + '"></span>' + p.label + (p.unit ? ' <span class="chip-u">' + p.unit + '</span>' : '') + '</span>';
      }).join("") + '</div>';
  }
  // 将任意序列按小时基线 resample 到 [0,dur]（区间内插值、超区间按末段斜率线性外推）
  function ptResampleTo(arr, hours, dur, npts) {
    if (!hours.length) return [];
    var lastH = hours[hours.length - 1], lastV = arr[arr.length - 1];
    var inc = arr.length > 1 ? (arr[arr.length - 1] - arr[arr.length - 2]) : 0;
    function interp(xv) {
      if (xv <= hours[0]) return arr[0];
      for (var i = 1; i < hours.length; i++) {
        if (hours[i] >= xv) { var r = (xv - hours[i - 1]) / ((hours[i] - hours[i - 1]) || 1); return arr[i - 1] + (arr[i] - arr[i - 1]) * r; }
      }
      return lastV + inc * (xv - lastH);
    }
    var pts = [], step = dur / Math.max(1, npts - 1);
    for (var k = 0; k < npts; k++) { var x = k * step; pts.push([x, interp(x)]); }
    return pts;
  }
  // 实测基线轨迹（原小时轴）
  function ptBaselineSeries(key, base, ts) {
    var hrs = ts.hours;
    if (key === "TITER") {
      var s = ts.DCW.reduce(function (a, c) { return a + c; }, 0); var acc = 0;
      var gf = expInductionFrac(), gate = gf * (hrs.length ? hrs[hrs.length - 1] : 0);
      return hrs.map(function (h, i) { if (h >= gate) acc += ts.DCW[i]; return (base.titer || 0) * (s ? acc / s : 0); });
    }
    if (key === "DCW") return ts.DCW.slice();
    return ts[key].slice();
  }
  // What-if 情景轨迹（baseline 经 modX/modQ 工艺偏移后投影到 [0,dur]）
  function ptScenarioSeries(key, base, ts, fa, dur, doRef) {
    var hrs = ts.hours;
    var gf = expInductionFrac(), gate = gf * (hrs.length ? hrs[hrs.length - 1] : 0);
    var xMax = fa && fa.xMax ? fa.xMax : HOST_XMAX.other;
    function xScen(i) { return Math.min(ts.DCW[i] * (fa ? fa.modX : 1), xMax); }
    if (key === "TITER") {
      var accS = 0, accB = 0, mQ = fa ? fa.modQ : 1;
      var tFull = hrs.map(function (h, i) {
        if (h >= gate) { accS += xScen(i); accB += ts.DCW[i]; }
        return (base.titer || 0) * mQ * (accB > 0 ? accS / accB : 0);
      });
      return ptResampleTo(tFull, hrs, dur, 60);
    }
    if (key === "DCW") return ptResampleTo(ts.DCW.map(function (v, i) { return xScen(i); }), hrs, dur, 60);
    if (key === "OUR") return ptResampleTo(ts.OUR.map(function (v) { return v * (fa ? fa.modX : 1); }), hrs, dur, 60);
    if (key === "CER") return ptResampleTo(ts.CER.map(function (v) { return v * (fa ? fa.modX : 1); }), hrs, dur, 60);
    if (key === "WCW") return ptResampleTo(ts.WCW.map(function (v, i) { return v * (xScen(i) / Math.max(ts.DCW[i], 1e-9)); }), hrs, dur, 60);
    if (key === "DO") { var sh = (state.process.do || doRef) - doRef; return ptResampleTo(ts.DO.map(function (v) { return v + sh; }), hrs, dur, 60); }
    if (key === "TEMP") { var mT = ts.TEMP.reduce(function (a, c) { return a + c; }, 0) / ts.TEMP.length; var shT = (state.process.temp || 30) - mT; return ptResampleTo(ts.TEMP.map(function (v) { return v + shT; }), hrs, dur, 60); }
    if (key === "FEED") return ptResampleTo(ts.FEED.map(function (v) { return v * (state.process.feed || 1); }), hrs, dur, 60);
    return ptResampleTo(ts[key], hrs, dur, 60); // pH/RPM/AIR/O2_FLOW/OUR/CER/WCW：简化模型未覆盖，仅按时长截取/外推
  }
  function ptWhatIfTraj(p, base, ts, bid) {
    var fa = ptScenarioFactors(p, base, bid);
    var hrs = ts ? ts.hours : [], dcwBase = ts ? ts.DCW : [];
    var xMax = fa.xMax, capped = false;
    var dcwScenFull = dcwBase.map(function (v) {
      var u = v * fa.modX;
      if (u > xMax) { capped = true; return xMax; }   // 生物学上限钳制（可达湿菌体/氧传递限制）
      return u;
    });
    var titerBase = base.titer || 0;
    var gf = expInductionFrac(), gate = gf * (hrs.length ? hrs[hrs.length - 1] : 0);
    // P(t) = titer_base · [∫_gate^t q_P·X dτ] / [∫_gate^T q_P0·X0 dτ]
    //      = titer_base · modQ · (∫_gate^t X_scen dτ) / (∫_gate^T X_base dτ)   —— modX、modQ 各计一次
    var accS = 0, accB = 0;
    var titerFull = hrs.map(function (h, i) {
      if (h >= gate) { accS += dcwScenFull[i]; accB += dcwBase[i]; }
      return titerBase * fa.modQ * (accB > 0 ? accS / accB : 0);
    });
    var dur = (p.dur != null && p.dur > 0) ? p.dur : (base.dur || 0);
    // 情景曲线严格落在 [0,dur]：缩短则裁剪、延长则按末段速率外推（修复“时长调整图谱无变化”）
    var dcwScen = ptResampleTo(dcwScenFull, hrs, dur, 60);
    var titerScen = ptResampleTo(titerFull, hrs, dur, 60);
    var endDCW = dcwScen.length ? dcwScen[dcwScen.length - 1][1] : 0;
    var titer = titerScen.length ? titerScen[titerScen.length - 1][1] : 0;
    var lastH = hrs.length ? hrs[hrs.length - 1] : 0;
    return {
      hours: hrs, dcwBase: dcwBase, dcwScen: dcwScen, titerScen: titerScen,
      titer: titer, endDCW: endDCW, dur: dur, lastH: lastH,
      rel: titerBase ? (titer / titerBase - 1) * 100 : 0,
      fa: fa, xMax: xMax, capped: capped, xMaxHost: ptXmax(ptHostOf(bid))
    };
  }
  var BATCH_COLORS = ["#1f9e89", "#2b6cb0", "#ffce4d", "#ff6b6b", "#b06ab3", "#4cc9f0", "#80ed99", "#ff8c42", "#7bdff2", "#a06cd5", "#ef476f", "#06d6a0"];
  function ptWiparMeta(k) {
    if (k === "TITER") return { label: "效价(产物)", color: "#80ed99", unit: "g/L" };
    var m = null; PT_PARAMS.forEach(function (p) { if (p.key === k) m = p; });
    return { label: m ? m.label : k, color: m ? m.color : "#9fb3c8", unit: m ? m.unit : "" };
  }
  // What-if 情景参数轨迹（可选）
  function ptWhatIfParamChips() {
    var keys = PT_PARAMS.map(function (p) { return p.key; }).concat(["TITER"]);
    return '<div class="sub" style="margin-top:10px">情景参数轨迹（可选 · 显示 基线 vs 情景 至设定发酵时长）</div><div style="display:flex;gap:8px;flex-wrap:wrap">' +
      keys.map(function (k) {
        var m = ptWiparMeta(k), on = state.process.wipar[k];
        return '<span class="chip' + (on ? " chip-on" : "") + '" onclick="FIP.set(\'process.wipar.' + k + '\', ' + (!on) + ')"><span class="chip-dot" style="background:' + m.color + '"></span>' + m.label + '</span>';
      }).join("") + '</div><div id="p_whatif_params" style="margin-top:10px"></div>';
  }
  // What-if 情景参数轨迹：合并为一张多轴大图（参考「在线过程曲线」），每个选中参数独立纵坐标，
  // 每条参数双线：基线（灰虚）/ 情景（实线，按 DO/温度/补料 偏移投影），统一时间轴至设定发酵时长
  function renderWhatIfScenarioChart(host, keys, base, ts, traj) {
    if (!host) return;
    if (!keys || !keys.length) { host.innerHTML = '<div class="hint">请在上方勾选情景参数（基线 vs 情景），将在此合并为一张多轴图（参考「在线过程曲线」）</div>'; return; }
    var hrs = ts.hours, dur = (traj && traj.dur) || (base && base.dur) || 0;
    var params = keys.map(function (k) {
      var meta = ptWiparMeta(k);
      var baseVals = ptBaselineSeries(k, base, ts);
      var basePts = ptResampleTo(baseVals, hrs, dur, 60);          // 基线重采样到 [0,dur]
      var scenPts = ptScenarioSeries(k, base, ts, traj.fa, dur, base.doRef); // 情景投影到 [0,dur]
      return { key: k, label: meta.label, unit: meta.unit, color: meta.color, basePts: basePts, scenPts: scenPts };
    });
    params.forEach(function (p) {
      var all = p.basePts.concat(p.scenPts).map(function (q) { return q[1]; });
      var lo = Math.min.apply(null, all), hi = Math.max.apply(null, all);
      if (lo === hi) hi = lo + 1;
      var pad = (hi - lo) * 0.12 || 1;
      p.ylo = lo - pad; p.yhi = hi + pad;
    });
    var n = params.length, colW = 30, leftN = Math.ceil(n / 2), rightN = n - leftN;
    var padT = 16, padB = 46, padL = 30 + leftN * colW, padR = 26 + rightN * colW;
    var innerW = 620, W = padL + innerW + padR, H = 380;
    var x0 = padL, x1 = W - padR, y0 = padT, y1 = H - padB;
    function X(v) { return x0 + (v / (dur || 1)) * (x1 - x0); }
    function Yy(t) { return y1 - t * (y1 - y0); }
    function buildPath(pts, lo, hi) {
      var d = "";
      pts.forEach(function (p, i2) { var t = (p[1] - lo) / (hi - lo); t = Math.max(0, Math.min(1, t)); d += (i2 === 0 ? "M" : "L") + X(p[0]).toFixed(1) + " " + Yy(t).toFixed(1) + " "; });
      return d;
    }
    var sb = [];
    sb.push('<svg viewBox="0 0 ' + W + ' ' + H + '" width="100%" height="' + H + '" preserveAspectRatio="xMidYMid meet" style="background:' + BG + ';border-radius:10px">');
    for (var g = 0; g <= 4; g++) {
      var gy = Yy(g / 4);
      sb.push('<line x1="' + x0 + '" y1="' + gy.toFixed(1) + '" x2="' + x1 + '" y2="' + gy.toFixed(1) + '" stroke="#17242c" stroke-width="1"/>');
    }
    if (ts.phases && ts.phases.length) {
      ts.phases.forEach(function (ph) {
        var xa = X(ph.start), xb = X(ph.end); if (xb - xa < 0.5) xb = xa + 0.5;
        sb.push('<rect x="' + xa.toFixed(1) + '" y="' + y0 + '" width="' + (xb - xa).toFixed(1) + '" height="' + (y1 - y0).toFixed(1) + '" fill="' + ph.color + '" fill-opacity="0.10"/>');
        sb.push('<line x1="' + xa.toFixed(1) + '" y1="' + y0 + '" x2="' + xa.toFixed(1) + '" y2="' + y1 + '" stroke="' + ph.color + '" stroke-width="1" stroke-dasharray="3 3" stroke-opacity="0.5"/>');
        if (xb - xa > 34) sb.push('<text x="' + (xa + 4).toFixed(1) + '" y="' + (y0 + 12) + '" fill="' + ph.color + '" font-size="10.5" opacity="0.9">' + esc(ph.label) + '</text>');
      });
    }
    // 设定发酵时长竖线
    var vx = X(dur);
    sb.push('<line x1="' + vx.toFixed(1) + '" y1="' + y0 + '" x2="' + vx.toFixed(1) + '" y2="' + y1 + '" stroke="#ff6b6b" stroke-width="1.5" stroke-dasharray="4 3"/>');
    sb.push('<text x="' + vx.toFixed(1) + '" y="' + (y0 - 2).toFixed(1) + '" fill="#ff6b6b" font-size="10.5" text-anchor="middle">设定 ' + fmt(dur, 0) + 'h</text>');
    for (var j = 0; j <= 4; j++) {
      var xv = dur * j / 4, xx = X(xv);
      sb.push('<line x1="' + xx.toFixed(1) + '" y1="' + y1 + '" x2="' + xx.toFixed(1) + '" y2="' + (y1 + 4) + '" stroke="#2a3a44" stroke-width="1"/>');
      sb.push('<text x="' + xx.toFixed(1) + '" y="' + (y1 + 18) + '" fill="#8fa3b8" font-size="11" text-anchor="middle">' + fmt(xv, 0) + '</text>');
    }
    params.forEach(function (s, i) {
      var lo = s.ylo, hi = s.yhi, left = i < leftN, k = left ? i : (i - leftN);
      var axX = left ? (x0 - (k + 1) * colW) : (x1 + (k + 1) * colW);
      sb.push('<line x1="' + axX.toFixed(1) + '" y1="' + y0 + '" x2="' + axX.toFixed(1) + '" y2="' + y1 + '" stroke="' + (s.color || PRIMARY) + '" stroke-width="1.4" stroke-opacity="0.85"/>');
      for (var t = 0; t <= 4; t++) {
        var tv = t / 4, yv = lo + (hi - lo) * tv, ty = Yy(tv);
        sb.push('<line x1="' + (axX - 2.5) + '" y1="' + ty.toFixed(1) + '" x2="' + (axX + 2.5) + '" y2="' + ty.toFixed(1) + '" stroke="' + (s.color || PRIMARY) + '" stroke-width="1" stroke-opacity="0.8"/>');
        sb.push('<text x="' + (left ? (axX - 4) : (axX + 4)).toFixed(1) + '" y="' + (ty + 3.5).toFixed(1) + '" fill="' + (s.color || PRIMARY) + '" font-size="9" text-anchor="' + (left ? "end" : "start") + '">' + fmt(yv, 1) + '</text>');
      }
      var midY = (y0 + y1) / 2;
      sb.push('<text x="' + (axX + (left ? -3 : 3)).toFixed(1) + '" y="' + midY.toFixed(1) + '" fill="' + (s.color || PRIMARY) + '" font-size="10" text-anchor="middle" transform="rotate(' + (left ? -90 : 90) + ' ' + axX.toFixed(1) + ' ' + midY.toFixed(1) + ')">' + esc(s.label) + (s.unit ? '(' + esc(s.unit) + ')' : '') + '</text>');
      sb.push('<path d="' + buildPath(s.basePts, lo, hi) + '" fill="none" stroke="#6b7d8f" stroke-width="1.6" stroke-dasharray="5 4" stroke-linejoin="round"/>');
      sb.push('<path d="' + buildPath(s.scenPts, lo, hi) + '" fill="none" stroke="' + (s.color || PRIMARY) + '" stroke-width="2" stroke-linejoin="round"/>');
    });
    sb.push('<line x1="' + x0 + '" y1="' + y1 + '" x2="' + x1 + '" y2="' + y1 + '" stroke="#3a4d59" stroke-width="1.5"/>');
    sb.push('<text x="' + ((x0 + x1) / 2) + '" y="' + (H - 6) + '" fill="#b7c7d6" font-size="12" text-anchor="middle">时间 (h) · 灰虚=基线 / 实线=情景（按设定偏移投影）至设定发酵时长</text>');
    sb.push("</svg>");
    var legend = '<div style="display:flex;gap:18px;flex-wrap:wrap;margin:8px 2px 2px;font-size:12px;color:#cdd9e5">' +
      '<span style="display:inline-flex;align-items:center;gap:6px"><span style="width:18px;height:0;border-top:2px dashed #6b7d8f;display:inline-block"></span>基线 (实测)</span>' +
      '<span style="display:inline-flex;align-items:center;gap:6px"><span style="width:18px;height:3px;border-radius:2px;background:' + PRIMARY + ';display:inline-block"></span>情景 (DO/温度/补料 偏移)</span></div>';
    params.forEach(function (s) {
      var lastB = s.basePts[s.basePts.length - 1][1], lastS = s.scenPts[s.scenPts.length - 1][1];
      legend += '<span style="display:inline-flex;align-items:center;gap:6px"><span style="width:14px;height:3px;border-radius:2px;background:' + (s.color || PRIMARY) + ';display:inline-block"></span>' + esc(s.label) + (s.unit ? ' <span style="color:#7d93a6">(' + esc(s.unit) + ')</span>' : '') + ' <span style="color:#9fb3c8">基线 ' + fmt(lastB, 1) + ' → 情景 ' + fmt(lastS, 1) + '</span></span>';
    });
    legend += "</div>";
    host.innerHTML = legend + sb.join("");
  }
  // 历史批次横向对比：参数下拉框（单选） + 批次下拉框（多选弹层） + 叠加图
  // 数据源直接来自 DATA.batches / PT_PARAMS，与批次数据库同步
  function ptCmpParamSelect() {
    return '<div class="sub" style="margin-top:6px">对比参数（下拉选择）</div>' +
      '<select id="cmp_param" class="fip-select" onchange="FIP.cmpParam(this.value)">' +
      PT_PARAMS.map(function (p) {
        return '<option value="' + p.key + '"' + (state.process.cmp.param === p.key ? " selected" : "") + '>' + p.label + (p.unit ? ' (' + p.unit + ')' : '') + '</option>';
      }).join("") + '</select>';
  }
  function ptCmpBatchSelect(filt) {
    var sel = state.process.cmp.batches;
    var opts = filt.map(function (b) {
      var on = sel.indexOf(b.batch_id) >= 0;
      var label = b.batch_id + ' · ' + (b.host_organism || "").replace(/_/g, " ") + ' · ' + (b.scale_tag || "");
      return '<label class="cmp-opt"><input type="checkbox" data-bid="' + b.batch_id + '" ' + (on ? "checked" : "") + ' onchange="FIP.toggleCmpBatch(\'' + b.batch_id + '\')"> <span>' + label + '</span></label>';
    }).join("");
    return '<div class="sub" style="margin-top:8px">对比批次（下拉多选 · 当前宿主过滤 ' + filt.length + ' 批）</div>' +
      '<div class="cmp-drop-wrap"><button id="cmp_batch_btn" class="fip-select" type="button" onclick="FIP.cmpToggleDrop(event)">已选 ' + sel.length + ' 批 ▾</button>' +
      '<div id="cmp_drop" class="cmp-drop">' + opts + '</div></div>';
  }
  function renderCompare() {
    var cmp = state.process.cmp, host = $("p_cmp");
    if (!host) return;
    var meta = ptWiparMeta(cmp.param);
    var bids = cmp.batches.length ? cmp.batches : [];
    var series = bids.map(function (bid, idx) {
      var ts = DATA.timeseries[bid];
      if (!ts || !ts[cmp.param]) return null;
      return { name: bid, color: BATCH_COLORS[idx % BATCH_COLORS.length], points: ts.hours.map(function (h, i) { return [h, ts[cmp.param][i]]; }) };
    }).filter(Boolean);
    if (!series.length) { host.innerHTML = '<div class="hint">请至少选择一个对比批次（点击上方批次标签）</div>'; return; }
    lineChart(host, series, { ytitle: meta.label + (meta.unit ? " (" + meta.unit + ")" : ""), ymin: 0, phases: null });
  }
  // 本次实验参数输入（实验设计）模块
  var VESSELS = ["2L", "3L", "5L", "10L", "50L", "100L", "200L", "500L", "10000L"];
  var EXP_PHASE_PARAMS = [
    { k: "ph", label: "pH", unit: "" },
    { k: "temp", label: "温度", unit: "℃" }
  ];
  var DO_CTRL_PARAMS = [
    { k: "air_vvm", label: "空气固定通气量", unit: "VVM" },
    { k: "do_set", label: "目标 DO 设定", unit: "%" },
    { k: "rpm_min", label: "搅拌转速最低", unit: "rpm" },
    { k: "rpm_max", label: "搅拌转速最高", unit: "rpm" },
    { k: "o2_start", label: "补氧启动阈值 DO", unit: "%" }
  ];
  var EXP_PHASES = [
    { k: "batch", label: "Batch 期" },
    { k: "fedbatch", label: "Fed-batch 期" },
    { k: "induction", label: "Induction 期" }
  ];
  function expFeedModeChips(ph, e) {
    var modes = [["none", "无补料"], ["constant", "恒定流速"], ["stepwise", "阶梯补料"]];
    return modes.map(function (m) {
      var on = e[ph.k].feed_mode === m[0];
      return '<span class="chip' + (on ? " chip-on" : "") + '" onclick="FIP.expSet(\'process.exp.' + ph.k + '.feed_mode\', \'' + m[0] + '\')">' + m[1] + '</span>';
    }).join("");
  }
  function expFeedStepTable(ph, e) {
    var steps = e[ph.k].feed_steps || [];
    var rows = steps.map(function (s, idx) {
      return '<tr><td><input type="number" step="any" value="' + s.t + '" onchange="FIP.expSet(\'process.exp.' + ph.k + '.feed_steps.' + idx + '.t\', parseFloat(this.value))"></td>' +
        '<td><input type="number" step="any" value="' + s.rate + '" onchange="FIP.expSet(\'process.exp.' + ph.k + '.feed_steps.' + idx + '.rate\', parseFloat(this.value))"></td>' +
        '<td><button class="btn2" onclick="FIP.expDelStep(\'' + ph.k + '\',' + idx + ')">移除</button></td></tr>';
    }).join("");
    return '<div style="color:#7d93a6;font-size:12px;margin-bottom:4px">流速单位 mL/kg/h（kg = 起始装液量 ' + e.fill_kg + ' kg）</div><table class="exp-step-tbl"><tr><th>起始时间 (h)</th><th>流速 (mL/kg/h)</th><th></th></tr>' + rows + '</table>' +
      '<button class="btn2" style="margin-top:6px" onclick="FIP.expAddStep(\'' + ph.k + '\')">+ 添加阶梯</button>';
  }
  // 每阶段补料策略库控件：保存 / 选择 / 删除 / 上传
  function feedStratCtl(ph) {
    return '<div class="feed-strat" style="margin-top:8px;display:flex;gap:6px;align-items:center;flex-wrap:wrap">' +
      '<span class="hint" style="margin:0">补料策略库：</span>' +
      '<input id="feed_strat_name_' + ph.k + '" type="text" placeholder="策略名" style="flex:0 1 110px;background:#0c1922;border:1px solid #24323c;color:#e6edf3;border-radius:8px;padding:5px 7px;font-size:12px">' +
      '<button class="btn2" onclick="FIP.feedStratSave(\'' + ph.k + '\')">💾 存为策略</button>' +
      '<select id="feed_strat_sel_' + ph.k + '" class="fip-select" onchange="FIP.feedStratApply(\'' + ph.k + '\', this.value)"><option value="">— 选择已保存策略 —</option></select>' +
      '<button class="btn2" onclick="FIP.feedStratDel(\'' + ph.k + '\')">🗑</button>' +
      '<label class="btn2" style="cursor:pointer;margin:0">📤 上传<input type="file" accept=".json,application/json" onchange="FIP.feedStratUpload(\'' + ph.k + '\', this.files[0])" style="display:none"></label>' +
      '</div>';
  }
  // 每阶段培养基配方库控件（与补料策略库并列，直接在本阶段卡片上「调用」已保存配方）
  function recipeLibCtl(ph) {
    return '<div class="feed-strat" style="margin-top:8px;display:flex;gap:6px;align-items:center;flex-wrap:wrap">' +
      '<span class="hint" style="margin:0">培养基配方库：</span>' +
      '<select id="recipe_lib_sel_' + ph.k + '" class="fip-select" onchange="FIP.expRecipeLibApply2(\'' + ph.k + '\', this.value)"><option value="">— 选择已保存配方 —</option></select>' +
      '</div>';
  }
  function recipeLibRender2() {
    EXP_PHASES.forEach(function (ph) {
      var sel = $("recipe_lib_sel_" + ph.k); if (!sel) return;
      var lib = recipeLibGet();
      sel.innerHTML = '<option value="">— 选择已保存配方 —</option>' + lib.map(function (x) {
        return '<option value="' + esc(x.name) + '">' + esc(x.name) + ' · ' + ((x.rows || []).length) + ' 组分</option>';
      }).join("");
    });
  }
  function renderExpForm() {
    var e = state.process.exp, host = $("p_exp");
    if (!host) return;
    var mt = expModelTimes(e); // 各阶段时序 + 预期 OD600（起始/末值）
    var phaseHtml = EXP_PHASES.map(function (ph) {
      if (!e.phases[ph.k]) return "";
      var p = e[ph.k];
      var grid = EXP_PHASE_PARAMS.map(function (pr) {
        return '<label class="exp-f"><span>' + pr.label + (pr.unit ? ' <i>' + pr.unit + '</i>' : '') + '</span>' +
          '<input type="number" step="any" value="' + p[pr.k] + '" onchange="FIP.expSet(\'process.exp.' + ph.k + '.' + pr.k + '\', parseFloat(this.value))"></label>';
      }).join("");
      if (ph.k === "batch") {
        grid += '<label class="exp-f"><span>碳源消耗速率 (g/L/h)</span><input type="number" step="any" value="' + (p.c_depletion_rate != null ? p.c_depletion_rate : "") + '" onchange="FIP.expSet(\'process.exp.batch.c_depletion_rate\', parseFloat(this.value))"></label>';
      }
      var feedBody = "";
      if (p.feed_mode === "constant") {
        feedBody = '<label class="exp-f" style="margin-top:8px"><span>恒定流速 (mL/kg/h)</span><input type="number" step="any" value="' + p.feed_const + '" onchange="FIP.expSet(\'process.exp.' + ph.k + '.feed_const\', parseFloat(this.value))"></label>';
      } else if (p.feed_mode === "stepwise") {
        feedBody = '<div style="margin-top:8px">' + expFeedStepTable(ph, e) + '</div>';
      }
      var _hk = hostKeyOf(e.host_type || "ecoli");
      var _medNote = (_hk === "pichia")
        ? '<div class="method-note" style="margin-bottom:6px"><b>BSM 基础盐全程恒定</b>；本阶段可变项为<b>碳源</b>：生长/流加期 = 甘油，诱导期 = 甲醇。请勿在此「换培养基」。</div>'
        : '<div class="method-note" style="margin-bottom:6px">该培养基全程恒定；碳源为葡萄糖。</div>';
      var med = _medNote + '<div class="sub" style="margin-top:8px">培养基信息（C:N）</div>' +
        (function () {
          var hk = hostKeyOf(e.host_type || "ecoli");
          var grp = (FermentSim.MEDIA && FermentSim.MEDIA[hk]) || {};
          var opts = Object.keys(grp).map(function (mk) {
            return '<option value="' + mk + '">' + esc(grp[mk].label) + "</option>";
          }).join("");
          return '<div style="display:flex;gap:8px;align-items:center;flex-wrap:wrap;margin-bottom:6px">' +
            '<span class="hint" style="margin:0">培养基备选项</span>' +
            '<select class="fip-select" onchange="FIP.expMediumFill(\'' + ph.k + "', this.value)\"><option value=\"\">— 选择 —</option>" + opts + "</select>" +
            "</div>";
        })() +
        '<div class="exp-grid">' +
        '<label class="exp-f"><span>碳源名称</span><input type="text" value="' + esc(p.c_source || "") + '" onchange="FIP.expSet(\'process.exp.' + ph.k + '.c_source\', this.value)"></label>' +
        '<label class="exp-f"><span>' + expCConcLabel(ph.k, _hk) + '</span><input type="number" step="any" value="' + esc(p.c_conc != null ? p.c_conc : "") + '" onchange="FIP.expSet(\'process.exp.' + ph.k + '.c_conc\', parseFloat(this.value))"></label>' +
        '<label class="exp-f"><span>氮源名称</span><input type="text" value="' + esc(p.n_source || "") + '" onchange="FIP.expSet(\'process.exp.' + ph.k + '.n_source\', this.value)"></label>' +
        '<label class="exp-f"><span>C:N 比例</span><input type="text" value="' + esc(p.cn_ratio || "") + '" onchange="FIP.expSet(\'process.exp.' + ph.k + '.cn_ratio\', this.value)"></label>' +
        '</div>' +
        '<div class="method-note" style="margin-top:6px">' + expCConcHint(ph.k, _hk) + '</div>' +
        '<div style="margin-top:8px;display:flex;gap:8px;align-items:center;flex-wrap:wrap">' +
        '<button class="btn2" onclick="FIP.expOpenRecipe(\'' + ph.k + '\')">📋 培养基配方（Excel 填写）</button>' +
        '<span class="hint" style="margin:0">已填 ' + ((p.recipe || []).length) + ' 组分</span></div>' +
        recipeLibCtl(ph);
      var header;
      if (ph.k === "batch") {
        header = ph.label + ' · 碳源耗尽判定 <span style="color:#2a8c6a;font-weight:700">≈ ' + fmt(mt.batchDur, 1) + ' h</span>（Fed-batch 起始）' +
          '<div style="margin-top:3px;font-size:12px;color:#9fb3c4">预期 OD600：起始 <b style="color:#cfe8ff">' + fmt(mt.od0, 1) + '</b> → 结束 <b style="color:#cfe8ff">' + fmt(mt.odBatchEnd, 1) + '</b></div>';
      } else if (ph.k === "induction") {
        header = ph.label + ' · 起始依据 ' + indTriggerChips(e) + ' <span style="color:#b07d2a;font-weight:700">起始 ≈ ' + fmt(mt.inductionStart, 1) + ' h</span> · 诱导时长 <input class="exp-dur" type="number" step="any" value="' + p.dur + '" onchange="FIP.expSet(\'process.exp.induction.dur\', parseFloat(this.value))"> h' +
          '<div style="margin-top:3px;font-size:12px;color:#9fb3c4">起始 OD600 ≈ <b style="color:#cfe8ff">' + fmt(mt.odFedEnd, 1) + '</b>（= Fed-batch 末值）</div>';
      } else {
        header = ph.label + ' · 时长 <input class="exp-dur" type="number" step="any" value="' + p.dur + '" onchange="FIP.expSet(\'process.exp.' + ph.k + '.dur\', parseFloat(this.value))"> h' +
          '<div style="margin-top:3px;font-size:12px;color:#9fb3c4">预期 OD600：起始 <b style="color:#cfe8ff">' + fmt(mt.odBatchEnd, 1) + '</b> → 结束 <b style="color:#cfe8ff">' + fmt(mt.odFedEnd, 1) + '</b>（Induction 起始）</div>';
      }
      var indBlock = (ph.k === "induction" && e.induction_trigger === "od600") ?
        '<label class="exp-f" style="margin-top:8px"><span>目标 OD600</span><input type="number" step="any" value="' + (e.target_od600 != null ? e.target_od600 : "") + '" onchange="FIP.expSet(\'process.exp.target_od600\', parseFloat(this.value))"></label>' : "";
      return '<div class="exp-phase"><div class="exp-phase-h">' + header + '</div>' +
        '<div class="exp-grid">' + grid + '</div>' +
        indBlock +
        '<div class="sub" style="margin-top:8px">补料模式</div>' +
        '<div style="display:flex;gap:8px;flex-wrap:wrap">' + expFeedModeChips(ph, e) + '</div>' +
        feedBody + feedStratCtl(ph) + med + '</div>';
    }).join("");
    var phaseToggles = EXP_PHASES.map(function (ph) {
      return '<label class="cmp-opt"><input type="checkbox" ' + (e.phases[ph.k] ? "checked" : "") + ' onchange="FIP.expTogglePhase(\'' + ph.k + '\')"> ' + ph.label + '</label>';
    }).join("");
    var doGrid = DO_CTRL_PARAMS.map(function (pr) {
      return '<label class="exp-f"><span>' + pr.label + (pr.unit ? ' <i>' + pr.unit + '</i>' : '') + '</span>' +
        '<input type="number" step="any" value="' + e.do_ctrl[pr.k] + '" onchange="FIP.expSet(\'process.exp.do_ctrl.' + pr.k + '\', parseFloat(this.value))"></label>';
    }).join("");
    var growthGrid = [
      ["od0", "初始 OD600", ""], ["mu_batch", "Batch 比生长速率 μ", "/h"], ["mu_fed", "Fed-batch 比生长速率 μ", "/h"], ["target_od600", "目标 OD600（诱导）", ""]
    ].map(function (g) {
      return '<label class="exp-f"><span>' + g[1] + (g[2] ? ' <i>' + g[2] + '</i>' : '') + '</span><input type="number" step="any" value="' + e[g[0]] + '" onchange="FIP.expSet(\'process.exp.' + g[0] + '\', parseFloat(this.value))"></label>';
    }).join("");
    var html =
      '<div class="sub" style="margin-top:6px">宿主与分子设计信息</div>' +
      '<div class="exp-grid">' +
      '<label class="exp-f"><span>宿主类型</span><select class="fip-select" onchange="FIP.expSet(\'process.exp.host_type\', this.value)">' +
      ['yeast', 'ecoli'].map(function (h) { return '<option value="' + h + '"' + (e.host_type === h ? " selected" : "") + '>' + (h === 'yeast' ? '🍶 酵母 Yeast' : '🦠 大肠杆菌 E.coli') + '</option>'; }).join("") + '</select></label>' +
      '<label class="exp-f"><span>发酵罐规模</span><select class="fip-select" onchange="FIP.expSet(\'process.exp.vessel\', this.value)">' +
      VESSELS.map(function (v) { return '<option value="' + v + '"' + (e.vessel === v ? " selected" : "") + '>' + v + '</option>'; }).join("") + '</select></label>' +
      '<label class="exp-f"><span>起始装液量 (kg)</span><input type="number" step="any" value="' + e.fill_kg + '" onchange="FIP.expSet(\'process.exp.fill_kg\', parseFloat(this.value))"></label>' +
      '</div>' +
      '<label class="exp-f" style="display:block;margin-top:8px"><span>分子设计 / 工程菌信息</span>' +
      '<input type="text" class="exp-strain" placeholder="如：PAOX1 启动子 + α-factor 信号肽；或 pET-28a + PelB 信号肽 + 目标蛋白" value="' + esc(e.strain) + '" onchange="FIP.expSet(\'process.exp.strain\', this.value)"></label>' +
      '<div class="sub" style="margin-top:12px">DO / 通气 / 搅拌 控制策略（跨阶段）</div>' +
      '<div class="method-note">空气保持固定通气量（= 设定 VVM × 起始装液量）；随菌体生长 DO 需求升高，搅拌转速由<b>最低</b>自动升至<b>最高</b>以跟踪目标 DO；若转速已达上限而 DO 仍低于<b>补氧阈值</b>，则自动提升氧气通气量。</div>' +
      '<div class="exp-grid">' + doGrid + '</div>' +
      '<div class="sub" style="margin-top:12px">生长 / 诱导模型参数</div>' +
      '<div class="method-note">Batch 时长由<b>碳源耗尽</b>推算（= 初始碳源浓度 ÷ 碳源消耗速率）；OD600 按指数模型 OD=OD₀·e^(μ·t) 增长；Induction 起始可设为<b>固定时间</b>（Fed-batch 结束）或<b>达到目标 OD600</b>。修改后下方曲线与时间轴自动重算。</div>' +
      '<div class="exp-grid">' + growthGrid + '</div>' +
      '<div class="sub" style="margin-top:12px">发酵阶段（勾选启用，分别设定 pH / 温度 / 补料）</div>' +
      '<div class="method-note">每阶段可分别设定 pH / 温度 / 补料，碳源与 C:N 见「培养基信息」卡。注意：以上 <b>pH / 温度 / 碳源浓度 / C:N 为实验设计记录</b>；当前仿真引擎按时间自动切换碳源（甘油→甲醇）、整程采用单一 pH/温度，<b>暂未耦合</b>这些阶段字段（见上方「生长/诱导模型参数」）。</div>' +
      '<div class="exp-phase-toggles">' + phaseToggles + '</div>' +
      '<div class="exp-phases" style="margin-top:8px">' + phaseHtml + '</div>' +
      '<div style="margin-top:10px;display:flex;gap:8px;flex-wrap:wrap">' +
      '<button class="btn2" onclick="FIP.expGenJSON()">生成实验方案 JSON</button>' +
      '<button class="btn2" onclick="FIP.expDownload()">下载 .json</button>' +
      '<button class="btn2" onclick="FIP.expSave()">💾 保存参数</button>' +
      '<button class="btn2" onclick="FIP.expLoad()">📂 载入已保存</button>' +
      '</div>' +
      '<div style="display:flex;gap:8px;flex-wrap:wrap;margin-top:10px;align-items:center">' +
        '<button class="btn" onclick="FIP.runExpSim()">▶ 用此设计运行机理仿真</button>' +
        '<span class="hint" style="margin:0">把三阶段 pH/温度/碳源浓度映射为分段曲线与碳源参数，在 What-if 仿真页跑出曲线（甲醇诱导流加速率沿用默认，可在仿真页微调）</span>' +
      '</div>' +
      '<div id="p_exp_savemsg" class="hint" style="margin-top:6px"></div>' +
      '<div class="sub" style="margin-top:14px">📉 补料流速曲线（实时 · 单位 mL/kg/h）</div>' +
      '<div id="p_feed_chart"></div>' +
      '<div id="p_feed_note" class="method-note"></div>' +
      '<pre id="p_exp_json" class="exp-json" style="display:none"></pre>';
    host.innerHTML = html;
    renderExpFeedChart();
    feedStratLibRender();
    recipeLibRender2();
  }
  // 生长 / 诱导时序模型：Batch 时长 = 碳源耗尽时刻；Induction 起始 = 固定时间 或 目标 OD600
  function expModelTimes(e) {
    var b = e.batch, fb = e.fedbatch, ind = e.induction;
    var cdr = b.c_depletion_rate > 0 ? b.c_depletion_rate : 0;
    var cc = b.c_conc > 0 ? b.c_conc : 0;
    var batchDur = (cdr > 0 && cc > 0) ? Math.max(1, cc / cdr) : Math.max(1, b.dur || 12);
    var od0 = e.od0 > 0 ? e.od0 : 0.5;
    var muB = e.mu_batch > 0 ? e.mu_batch : 0.18;
    var muF = e.mu_fed > 0 ? e.mu_fed : 0.10;
    var odBatchEnd = od0 * Math.exp(muB * batchDur);
    var fedDur = 0, inductionStart = batchDur;
    if (e.phases.fedbatch) {
      if (e.induction_trigger === "od600" && e.phases.induction) {
        var tgt = e.target_od600 > 0 ? e.target_od600 : 30;
        fedDur = (tgt > odBatchEnd && muF > 0) ? Math.max(0.5, Math.log(tgt / odBatchEnd) / muF) : Math.max(0.5, fb.dur || 12);
      } else {
        fedDur = Math.max(0.5, fb.dur || 12);
      }
      inductionStart = batchDur + fedDur;
    }
    var indDur = e.phases.induction ? Math.max(0, ind.dur || 0) : 0;
    var odFedEnd = e.phases.fedbatch ? odBatchEnd * Math.exp(muF * fedDur) : odBatchEnd; // Fed-batch 末 OD600（= Induction 起始）
    return { batchDur: batchDur, odBatchEnd: odBatchEnd, fedDur: fedDur, inductionStart: inductionStart, indDur: indDur, total: inductionStart + indDur, od0: od0, muB: muB, muF: muF, odFedEnd: odFedEnd, trigger: e.induction_trigger };
  }
  // 碳源浓度字段语义按阶段区分（避免"批式初始浓度"与"补料液浓度"混用导致机理引擎参数量级错误）
  // batch：初始发酵液中的碳源浓度 g/L；fedbatch / induction：补料液的碳源浓度 g/L
  function expCConcLabel(phk, host) {
    if (phk === "batch") return "初始碳源浓度 <i>g/L 发酵液</i>";
    if (phk === "induction") return host === "pichia" ? "甲醇补料液浓度 <i>g/L（纯甲醇 792）</i>" : "诱导期补料液浓度 <i>g/L</i>";
    return "补料液碳源浓度 <i>g/L 补料液</i>";
  }
  // 语义迁移：旧版把 fedbatch / induction 的 c_conc 记为「发酵液内碳源浓度」（典型 5 g/L），
  // 新版统一为「补料液碳源浓度」（g/L，常用 400–800）。若载入旧存档或初始化时仍为 <50 g/L，
  // 直接按补料液浓度解释会把体积流率放大百倍（装液量爆掉），故一次性纠正到新默认值。
  function expNormalizeCConc(e) {
    var changed = false;
    ["fedbatch", "induction"].forEach(function (k) {
      var v = parseFloat(e[k] && e[k].c_conc);
      if (isFinite(v) && v > 0 && v < 50) {
        e[k].c_conc = (k === "induction" && hostKeyOf(e.host_type || "ecoli") === "pichia") ? 792 : 500;
        changed = true;
      }
    });
    return changed;
  }
  // 补料液浓度合理性告警（不进入模型，仅提示）
  function expFeedConcWarn(e) {
    var out = [];
    var host = hostKeyOf(e.host_type || "ecoli");
    if (e.phases.fedbatch) {
      var v = parseFloat(e.fedbatch.c_conc);
      if (isFinite(v) && v > 0 && v < 200) out.push("Fed-batch 补料液浓度 " + v + " g/L 偏低：体积流率过大、稀释严重（常用 400–600 g/L）");
    }
    if (e.phases.induction && host === "pichia") {
      var m = parseFloat(e.induction.c_conc);
      if (isFinite(m) && m > 0 && m < 300) out.push("甲醇补料液浓度 " + m + " g/L 偏低：带入大量水分（纯甲醇 792 g/L）");
    }
    return out;
  }
  function expCConcHint(phk, host) {
    if (phk === "batch") return "批式起始加入发酵液的碳源量，直接决定 Batch 时长（浓度 ÷ 消耗速率）";
    if (phk === "induction") {
      return host === "pichia"
        ? "诱导期按甲醇补料液浓度换算体积流率（纯甲醇 792 g/L 为基准；浓度越低带入水分越多、稀释越大）"
        : "诱导期补料液浓度：与流加期不同时按此值计算体积流率（质量流率由 μ<sub>fed</sub> 决定）";
    }
    return "流加液碳源浓度：质量补料速率由 μ<sub>fed</sub> 决定，浓度只改变体积流率与稀释程度";
  }
  // 实验设计 → 机理仿真 recipe：把三阶段 ph/temp 映射为分段曲线；batch 碳源浓度 → 初始碳源，
  // fedbatch / induction 碳源浓度 → 补料液浓度（毕赤诱导期为甲醇补料液浓度，换算体积流率）。
  function expToSimRecipe(e) {
    var host = hostKeyOf(e.host_type || "ecoli");
    var mt = expModelTimes(e);
    var batchEnd = mt.batchDur, indStart = mt.inductionStart, total = Math.max(mt.total, 1);
    var phDef = host === "pichia" ? 5.5 : 7.0, tDef = host === "pichia" ? 28 : 37;
    var phPts = [], tempPts = [];
    function segPush(t0, p) {
      if (p && p.ph != null) phPts.push([t0, p.ph]);
      if (p && p.temp != null) tempPts.push([t0, p.temp]);
    }
    if (e.phases.batch) segPush(0, e.batch);
    if (e.phases.fedbatch) segPush(batchEnd, e.fedbatch);
    if (e.phases.induction) segPush(indStart, e.induction);
    var rec = { host: host, t_end_h: total, ph_curve: phPts, temp_curve: tempPts };
    var muFed = e.mu_fed > 0 ? e.mu_fed : 0.10;
    if (host === "pichia") {
      rec.s0_g_l = e.batch.c_conc > 0 ? e.batch.c_conc : 20;   // JS 引擎：gly0_g_l 读自 s0_g_l
      rec.gly_feed_start_h = batchEnd;
      rec.gly_feed_end_h = batchEnd + mt.fedDur;
      rec.meoh_start_h = indStart;
      rec.mu_set = muFed;                                      // JS 引擎：mu_set_gly 读自 mu_set
      rec.feed_s_g_l = e.fedbatch.c_conc > 0 ? Math.min(e.fedbatch.c_conc, 1200) : 500;  // 甘油补料液浓度
      // 诱导期甲醇：按补料液浓度换算体积流率（默认 0.012 L/h 对应纯甲醇 792 g/L）
      var mConc = e.induction.c_conc > 0 ? Math.max(50, Math.min(e.induction.c_conc, 792)) : 792;
      rec.meoh_feed_rate_lh = 0.012 * 792 / mConc;
      rec.temp_pre = e.batch.temp != null ? e.batch.temp : tDef;
      rec.ph_set = e.batch.ph != null ? e.batch.ph : phDef;
      rec.medium = "bsm_methanol";
    } else {
      rec.s0_g_l = e.batch.c_conc > 0 ? e.batch.c_conc : 20;
      rec.feed_s_g_l = e.fedbatch.c_conc > 0 ? Math.min(e.fedbatch.c_conc, 900) : 500;   // 流加液葡萄糖浓度
      rec.feed_s_g_l_post = e.induction.c_conc > 0 ? Math.min(e.induction.c_conc, 900) : 0; // 0 = 同流加期
      rec.feed_start_h = batchEnd;
      rec.induction_h = indStart;
      rec.mu_set = muFed;
      rec.temp_pre = e.batch.temp != null ? e.batch.temp : tDef;
      rec.temp_post = e.induction.temp != null ? e.induction.temp : 30;
      rec.ph_set = e.batch.ph != null ? e.batch.ph : phDef;
      rec.medium = "defined";
    }
    var vol = e.fill_kg > 0 ? e.fill_kg : 10;               // 近似：1 kg ≈ 1 L
    rec.batch_volume_l = vol;
    rec.x0_g_l = Math.max(0.05, (e.od0 > 0 ? e.od0 : 0.5) * 0.4);   // OD600 → g DCW/L（≈0.4 g/L per OD）
    rec.airflow_lmin = Math.max(1, vol);                    // 按 1 vvm 同步通气量（与基线 10 L / 10 L·min⁻¹ 一致）
    rec._cn_advisory = expCnAdvisory(e).concat(expFeedConcWarn(e));
    return rec;
  }
  // C:N 比软提示（不进入机理模型，仅作为仿真页提示）：过高 → 氮限制/溢流风险
  function expCnAdvisory(e) {
    var out = [];
    ["batch", "fedbatch", "induction"].forEach(function (ph) {
      if (!e.phases[ph]) return;
      var r = String(e[ph].cn_ratio || "");
      var m = r.match(/^\s*([\d.]+)\s*:\s*1\s*$/);
      if (!m) return;
      var v = parseFloat(m[1]);
      if (!isFinite(v)) return;
      var label = { batch: "Batch", fedbatch: "Fed-batch", induction: "Induction" }[ph];
      if (v > 20) out.push(label + " C:N=" + r + " 偏高（碳过剩/氮限制风险）");
      else if (v < 3) out.push(label + " C:N=" + r + " 偏低（碳限制、生物量不足风险）");
    });
    return out;
  }
  function runExpSim() {
    var e = state.process.exp;
    var rec = expToSimRecipe(e);
    state.sim = state.sim || {};
    state.sim.host = rec.host;
    state.sim.medium = rec.medium;
    state.sim.expRecipe = rec;
    state.sim.usingExp = true;
    FIP.nav("simanim");
  }

  // 历史批次 → 仿真 recipe：批次摘要 + 工艺设定点（build_interactive 已注入 DATA.batches）
  function batchToSimRecipe(b) {
    var host = hostKeyOf(b.host_organism || "");
    var isP = host === "pichia";
    var vol = num(b.working_volume_l, 10);
    var tEnd = num(b.duration_h, isP ? 96 : 24);
    var med = batchMedium(b) || (FermentSim.DEFAULT_MEDIUM && FermentSim.DEFAULT_MEDIUM[host]) || "defined";
    var rec = {
      host: host, batch_volume_l: vol, t_end_h: tEnd, medium: med,
      ph_set: num(b.ph_mean, isP ? 5.5 : 7.0),
      rpm: num(b.rpm_mean, 800), airflow_lmin: num(b.air_mean, Math.max(1, vol)),
      mu_set: isP ? 0.08 : 0.11,
      kla_scale: 1.0
    };
    if (isP) {
      rec.gly0_g_l = num(b.s0_g_l, 20);
      rec.temp_pre = 28;
      rec.temp_c = num((b.induc_temp_mean != null ? b.induc_temp_mean : b.temp_mean), 28);
      rec.meoh_start_h = num(b.induc_start_h, 34);
      rec.meoh_feed_rate_lh = 0.012;
      rec.gly_feed_start_h = 18; rec.gly_feed_end_h = 32;
      rec.feed_s_g_l = num(b.feed_carbon_glc, 500);
    } else {
      rec.temp_pre = 37;
      rec.feed_start_h = 7;
      rec.induction_h = num(b.induc_start_h, 16);
      rec.temp_post = num((b.induc_temp_mean != null ? b.induc_temp_mean : b.temp_mean), 30);
      rec.feed_s_g_l = num(b.feed_carbon_glc, 500);
    }
    return rec;
  }

  function simLoadBatch(bid) {
    if (!bid) return;
    var b = batchById(bid); if (!b) return;
    var host = hostKeyOf(b.host_organism || "");
    state.sim.host = host;
    state.sim.medium = batchMedium(b) || (FermentSim.DEFAULT_MEDIUM && FermentSim.DEFAULT_MEDIUM[host]) || "defined";
    state.sim.loadedRecipe = batchToSimRecipe(b);
    state.sim.loadedFrom = "batch";
    state.sim.srcBatch = bid;
    state.sim.expRecipe = null;
    if (state.sim._body) renderSimAnim(state.sim._body); else FIP.nav("simanim");
  }

  function simAdoptExp() {
    var e = state.process.exp;
    if (!e) { alert("请先在 Process Twin 页填写「本次实验参数输入」"); return; }
    var rec = expToSimRecipe(e);
    state.sim.host = rec.host;
    state.sim.medium = rec.medium;
    state.sim.loadedRecipe = rec;
    state.sim.loadedFrom = "exp";
    state.sim.srcBatch = null;
    state.sim.expRecipe = null;
    if (state.sim._body) renderSimAnim(state.sim._body); else FIP.nav("simanim");
  }

  function simSourceBadge() {
    if (state.sim.srcBatch) return "📚 历史批次 " + esc(state.sim.srcBatch);
    if (state.sim.loadedFrom === "exp") return "🧪 本次实验参数输入";
    if (state.sim.loadedRecipe) return "📥 外部加载";
    return "✋ 手动滑块";
  }
  function indTriggerChips(e) {
    var opts = [["time", "固定时间(批次结束)"], ["od600", "目标 OD600"]];
    return opts.map(function (o) {
      var on = e.induction_trigger === o[0];
      return '<span class="chip' + (on ? " chip-on" : "") + '" onclick="FIP.expSet(\'process.exp.induction_trigger\', \'' + o[0] + '\')">' + o[1] + '</span>';
    }).join("");
  }
  // 诱导起始占全程的比例（无 Induction 阶段时为 0 = 不门控，产物从 0 时刻累积）
  function expInductionFrac() {
    var e = state.process.exp;
    if (!e.phases.induction) return 0;
    var mt = expModelTimes(e);
    if (mt.total <= 0) return 0;
    return Math.max(0, Math.min(1, mt.inductionStart / mt.total));
  }
  // 补料流速曲线：根据实验设计中各阶段补料模式（无 / 恒定 / 阶梯）实时绘制 mL/kg/h 随时间变化
  // 横轴为累计发酵时间（各启用阶段时长相加），纵轴为补料流速，阶梯补料按 {t, rate} 绘制阶梯；背景带区分 Batch/Fed-batch/Induction
  function renderExpFeedChart() {
    var host = $("p_feed_chart"); if (!host) return;
    var note = $("p_feed_note");
    var e = state.process.exp, mt = expModelTimes(e);
    var pts = [], phases = [], t0 = 0, total = mt.total, maxR = 0;
    var seg = [];
    if (e.phases.batch) seg.push({ k: "batch", dur: mt.batchDur, p: e.batch });
    if (e.phases.fedbatch) seg.push({ k: "fedbatch", dur: mt.fedDur, p: e.fedbatch });
    if (e.phases.induction) seg.push({ k: "induction", dur: mt.indDur, p: e.induction });
    seg.forEach(function (s) {
      var p = s.p, dur = Math.max(0, s.dur), t1 = t0 + dur;
      var col = s.k === "batch" ? "#3b6ea5" : (s.k === "fedbatch" ? "#2a8c6a" : "#b07d2a");
      var plabel = s.k === "batch" ? "Batch" : (s.k === "fedbatch" ? "Fed" : "Ind");
      phases.push({ start: t0, end: t1, label: plabel, color: col });
      var mode = p.feed_mode || "none";
      if (mode === "none") {
        pts.push([t0, 0]); pts.push([t1, 0]);
      } else if (mode === "constant") {
        var r = Math.max(0, p.feed_const || 0); if (r > maxR) maxR = r;
        pts.push([t0, r]); pts.push([t1, r]);
      } else {
        var ss = (p.feed_steps || []).slice().sort(function (a, b) { return (a.t || 0) - (b.t || 0); });
        if (!ss.length) { pts.push([t0, 0]); pts.push([t1, 0]); }
        else {
          var prev = 0; pts.push([t0, prev]);
          ss.forEach(function (st) {
            var stt = t0 + Math.max(0, st.t || 0);
            if (stt > t0) { pts.push([stt, prev]); }
            pts.push([stt, st.rate]); if (st.rate > maxR) maxR = st.rate; prev = st.rate;
          });
          pts.push([t1, prev]);
        }
      }
      t0 = t1;
    });
    // OD600 轨迹（指数模型）：0..batchDur 用 μ_batch；此后用 μ_fed；并叠加目标 OD600 参考线
    var odPts = [], odMax = 0, odTarget = (mt.trigger === "od600" && e.target_od600 > 0) ? e.target_od600 : 0;
    if (seg.length && total > 0) {
      var N = 60;
      for (var i = 0; i <= N; i++) {
        var tt = total * i / N, od;
        if (tt <= mt.batchDur) od = mt.od0 * Math.exp(mt.muB * tt);
        else od = mt.odBatchEnd * Math.exp(mt.muF * (tt - mt.batchDur));
        odPts.push([tt, od]); if (od > odMax) odMax = od;
      }
    }
    if (!pts.length || total <= 0) {
      host.innerHTML = '<div class="method-note">启用至少一个发酵阶段并设置补料模式（恒定流速 / 阶梯补料）后，此处生成补料流速 / 累计流量曲线（单位 mL/kg/h），并叠加 OD600 轨迹与目标线。</div>';
      if (note) note.innerHTML = "";
      return;
    }
    var cum = 0, cumPts = [[0, 0]];
    for (var ci = 1; ci < pts.length; ci++) {
      var a = pts[ci - 1], bb = pts[ci], dt = bb[0] - a[0];
      if (dt > 0) cum += a[1] * dt;
      cumPts.push([bb[0], cum]);
    }
    var series = [
      { name: "补料流速", color: "#ff6b6b", unit: "mL/kg/h", ylo: 0, yhi: maxR > 0 ? maxR * 1.1 : 1, points: pts },
      { name: "累计流量", color: "#ffce4d", unit: "mL/kg", ylo: 0, yhi: cum > 0 ? cum * 1.1 : 1, points: cumPts }
    ];
    if (odPts.length) {
      series.push({ name: "OD600", color: "#4ad0c0", unit: "", ylo: 0, yhi: Math.max(odMax, odTarget, 10) * 1.1, points: odPts });
      if (odTarget > 0) series.push({ name: "目标 OD600", color: "#4ad0c0", unit: "", ylo: 0, yhi: Math.max(odMax, odTarget, 10) * 1.1, dash: true, points: [[0, odTarget], [total, odTarget]] });
      // 预期 OD600 相位边界参考线（Batch 末 = Fed-batch 起；Fed-batch 末 = Induction 起）
      series.push({ name: "OD₆₀₀ Batch末值", color: "#7d93a6", unit: "", ylo: 0, yhi: Math.max(odMax, odTarget, 10) * 1.1, dash: true, points: [[0, mt.odBatchEnd], [total, mt.odBatchEnd]] });
      if (e.phases.fedbatch) series.push({ name: "OD₆₀₀ Fed末值", color: "#9fb3c4", unit: "", ylo: 0, yhi: Math.max(odMax, odTarget, 10) * 1.1, dash: true, points: [[0, mt.odFedEnd], [total, mt.odFedEnd]] });
    }
    lineChartMulti(host, series, { phases: phases, xtitle: "时间 (h)" });
    if (note) {
      var enabled = EXP_PHASES.filter(function (ph) { return e.phases[ph.k]; }).map(function (ph) {
        var p = e[ph.k], m = p.feed_mode || "none";
        var txt = m === "none" ? "无补料" : (m === "constant" ? ("恒定 " + fmt(p.feed_const || 0, 2) + " mL/kg/h") : ("阶梯 " + (p.feed_steps || []).length + " 段"));
        return ph.label + "·" + txt;
      }).join(" ｜ ");
      var trigTxt = (e.phases.induction) ? (mt.trigger === "od600" ? ("目标 OD600 " + fmt(e.target_od600, 1) + "（届时 OD600≈" + fmt(mt.odBatchEnd * Math.exp(mt.muF * mt.fedDur), 1) + "）") : "固定时间（Fed-batch 结束）") : "—";
      note.innerHTML = "Batch 碳源耗尽 ≈ " + fmt(mt.batchDur, 1) + " h → Fed-batch 起始；Induction 起始 ≈ " + fmt(mt.inductionStart, 1) + " h（依据：" + trigTxt + "）；总时长 " + fmt(total, 1) + " h · 起始装液量 " + fmt(e.fill_kg, 2) + " kg。预期 OD600 末值：Batch ≈ <b>" + fmt(mt.odBatchEnd, 1) + "</b>，Fed-batch ≈ <b>" + fmt(mt.odFedEnd, 1) + "</b>。各阶段补料：" + enabled + "。累计流量末值 " + fmt(cum, 1) + " mL/kg（× 起始装液量 " + fmt(e.fill_kg, 2) + " kg ≈ " + fmt(cum * e.fill_kg, 1) + " mL 总补料量）。";
    }
  }
  // What-if 仿真采纳实验设计模块设定：目标DO / 末阶段温度 / 补料流速 / 总时长
  function ptAdoptExp() {
    var e = state.process.exp;
    if (e.do_ctrl && e.do_ctrl.do_set != null) state.process.do = e.do_ctrl.do_set;
    var en = EXP_PHASES.filter(function (ph) { return e.phases[ph.k]; });
    if (en.length) {
      var lastPh = en[en.length - 1];
      if (e[lastPh.k].temp != null) state.process.temp = e[lastPh.k].temp;
      var feedPh = en.filter(function (ph) { return (e[ph.k].feed_mode || "none") !== "none"; });
      var fp = feedPh.length ? feedPh[0] : null;
      state.process.feed = fp ? (e[fp.k].feed_const || 0) : 0;
      var dur = expModelTimes(e).total;
      if (dur > 0) state.process.dur = dur;
    }
    if (state.selBatch) viewProcess(state.selBatch);
  }
  function renderProcess(body) {
    var hostOpts = [
      ["all", "全部宿主"],
      ["yeast", "🍶 酵母 Yeast"],
      ["ecoli", "🦠 大肠杆菌 E.coli"]
    ];
    // 宿主过滤后可选批次
    var filt = DATA.batches.filter(function (b) { return state.process.host === "all" || ptHostType(b) === state.process.host; });
    if (!filt.length) filt = DATA.batches;
    if (!state.selBatch || !filt.some(function (b) { return b.batch_id === state.selBatch; })) state.selBatch = filt[0].batch_id;
    if (!state.process.cmp.batches.length) state.process.cmp.batches = filt.slice(0, Math.min(4, filt.length)).map(function (b) { return b.batch_id; });
    var linkHost = state.process.host === "yeast" ? "yeast" : (state.process.host === "ecoli" ? "ecoli" : "");
    var html = card(
      '<div class="sub">宿主类型（对应 Expression Twin）</div>' +
      '<div style="display:flex;gap:8px;flex-wrap:wrap;align-items:center">' +
      hostOpts.map(function (h) {
        var on = state.process.host === h[0];
        return '<span class="chip' + (on ? " chip-on" : "") + '" onclick="FIP.set(\'process.host\', \'' + h[0] + '\')">' + h[1] + '</span>';
      }).join("") +
      (linkHost ? ' <button class="btn2" style="margin-left:4px" onclick="FIP.nav(\'' + linkHost + '\')">↗ 打开对应 Expression Twin</button>' : '') +
      '</div>' +
      '<div class="sub" style="margin-top:10px">批次选择（' + filt.length + ' 批）</div>' +
      '<div style="display:flex;gap:10px;flex-wrap:wrap"><select id="p_batch" onchange="FIP.viewProcess(this.value)">' +
      filt.map(function (b) { return '<option value="' + b.batch_id + '">' + b.batch_id + ' · ' + b.host_organism.replace("_", " ") + ' · ' + b.scale_tag + '</option>'; }).join("") + '</select></div>'
    );
    html += card(`<div class="sub">📈 历史批次横向对比（单参数 · 多批次叠加）</div>` +
      ptCmpParamSelect() + ptCmpBatchSelect(filt) +
      `<div id="p_cmp" style="margin-top:10px"></div>` +
      `<div class="method-note">同一参数在多个历史批次上的轨迹叠加对比（数据源与批次数据库同步，筛选当前宿主）。下拉选择参数与批次即时刷新。</div>`);
    // 实验设计信息量反馈 (DOE)：按所选批次宿主/尺度展示工艺设定点覆盖（show_construct=false，对应 Streamlit）
    var _b0 = DATA.batches.filter(function (x) { return x.batch_id === state.selBatch; })[0];
    if (_b0) {
      html += card(doeSectionHtml({
        host: _b0.host_organism, showConstruct: false, scale: _b0.scale_tag, withHostSelect: false
      }));
    }
    html += '<div id="p_out"></div>';
    body.innerHTML = html;
    if (state.selBatch) viewProcess(state.selBatch);
  }

  // 终点预测的可信度披露：置信区间已升级为共形分位数回归（CQR），经概率校准
  // （保外覆盖率≈名义 95%，宽度随条件误差自适应）。下方同时给出留批次交叉验证的
  // 真实误差，避免用户把区间宽度误当成"准确度"。
  function ptModelQualityRow(target, cn) {
    var Q = DATA.model_quality;
    if (!Q || !Q.cv || !Q.cv[target]) return "";
    var v = Q.cv[target];
    var r2 = v.cv_r2_mean;
    var col = r2 >= 0.85 ? "#3a6" : (r2 >= 0.70 ? "#e0a94d" : "#e06c6c");
    var lvl = r2 >= 0.85 ? "高" : (r2 >= 0.70 ? "中" : "偏低");
    return `<tr><td style="color:#8fa3b8">${cn} 模型可信度</td>` +
      `<td colspan="2"><span style="color:${col};font-weight:600">CV R²=${fmt(r2, 3)}（${lvl}）</span>` +
      `　留批次交叉验证 MAE≈${fmt(v.cv_mae, 2)}（目标量标准差 ${fmt(v.y_std, 2)}，相对误差 ${fmt(v.cv_mae / v.y_std * 100, 1)}%）` +
      `　OOB R²=${v.oob_r2 == null ? "—" : fmt(v.oob_r2, 3)}</td></tr>`;
  }
  // 置信区间校准状态（路线图 #8）：共形分位数回归，保外覆盖率≈名义 95%。
  function ptCICalibrationRow(target, cn) {
    var Q = DATA.model_quality, ci = Q && Q.ci ? Q.ci[target] : null;
    if (!ci) return "";
    var cov = ci.empirical_coverage, old = ci.old_rfvar_coverage;
    var ok = cov != null && cov >= 0.90 && cov <= 0.99;
    var col = ok ? "#3a6" : "#e0a94d";
    var verb = ok ? "已校准" : "需复核";
    return `<tr><td style="color:#8fa3b8">${cn || target} 区间校准</td>` +
      `<td colspan="2"><span style="color:${col};font-weight:600">共形分位数回归 · ${verb}</span>` +
      `　保外覆盖率≈${(cov * 100).toFixed(0)}%（名义 95%）` +
      (old != null ? `　旧 RF 树方差法≈${(old * 100).toFixed(0)}%（未校准代理）` : "") +
      `</td></tr>`;
  }
  function viewProcess(bid) {
    state.selBatch = bid;
    var ts = DATA.timeseries[bid];
    var pred = DATA.predictions[bid];
    var b = DATA.batches.filter(function (x) { return x.batch_id === bid; })[0];
    var kpi = ts ? ptKPI(ts) : null;
    var SM = ts ? ptSoftMeta(bid) : ptSoftMeta(null);
    // 三法交叉校核：OUR 氧衡算 / CER 碳衡算（按 RQ 折算为等效耗氧）/ Logistic 生长律拟合（不依赖尾气）。
    var doOur = ts ? ptDcwFromGas(ts.hours, ts.OUR, SM.y_ox, SM.mo2, SM.x0) : [];
    var doCer = ts ? ptDcwFromGas(ts.hours, ts.CER, SM.y_ox / SS_RQ, SM.mo2, SM.x0) : [];
    var doExp = ts ? ptDcwLogistic(ts.hours, ts.DCW) : [];
    var refBid = state.process.ref;
    if (!refBid) {
      DATA.batches.forEach(function (x) { if (!refBid && x.success && x.batch_id !== bid) refBid = x.batch_id; });
      if (!refBid) DATA.batches.forEach(function (x) { if (!refBid && x.batch_id !== bid) refBid = x.batch_id; });
      if (!refBid) refBid = DATA.batches[0].batch_id;
    }
    var dev = (b && ts && refBid && refBid !== bid) ? ptGolden(b, ts, kpi, refBid) : null;
    // 情景基准：DO/温度取批次实际工况（诱导后中位），而非末值（末值受终点氧耗崩溃影响会失真）
    var base = b ? {
      titer: b.titer_g_l, dur: b.duration_h,
      doRef: ts ? ptDoRef(ts) : 30, tempRef: ts ? ptTempRef(ts) : 30
    } : null;
    // 切换批次时把 What-if 滑块同步为该批次的**实际工况**（DO/温度/补料/时长），
    // 使「未调节」的情景严格等于基线（modX=modQ=1），用户调节后才是真实的 delta。
    if (base && (state.process._lastBid !== bid || state.process.dur == null || !(state.process.dur > 0))) {
      state.process.dur = base.dur;
      state.process.do = Math.max(10, Math.min(60, Math.round(base.doRef)));
      state.process.temp = Math.max(20, Math.min(37, Math.round(base.tempRef * 10) / 10));
      state.process.feed = 1;
    }
    state.process._lastBid = bid;
    var wi = base ? ptWhatIf(state.process, base, bid) : null;

    var html = card(`<div class="sub">在线过程曲线（参数可选 · 分期标注）</div>` + ptParamChips() +
      `<div class="ysub" style="margin:6px 0 4px">背景带：<span style="color:#3b6ea5">■</span> batch 期 · <span style="color:#2a8c6a">■</span> fed-batch 期 · <span style="color:#b07d2a">■</span> 诱导/生产期（派生）。所有勾选参数合并于同一时间轴，<b>每个参数拥有独立纵坐标（真实量纲，左右两侧分列着色）</b>，颜色对应参数，可直接读各自绝对数值，跨量纲同步对比。</div>` +
      `<div id="p_lines"></div>`);
    html += card(`<div class="sub">🧫 软传感器 Biomass · 多方法一致性</div>` +
      `<div class="ysub" style="margin-bottom:6px">① OUR 氧衡算反演 / ② CER 碳衡算反演 / ③ Logistic 生长律拟合（不依赖尾气）三法对比，集成软测量为基准</div>` +
      `<div id="p_soft"></div><div id="p_softinfo"></div>`);
    html += card(`<div class="sub">📊 工艺 KPI 仪表盘</div><div id="p_kpi"></div>`);
    html += card(`<div class="sub">🎯 Golden Batch 偏离分析</div>` +
      `<div style="display:flex;gap:10px;align-items:center;margin-bottom:6px"><label class="ksub">参考标杆批</label>` +
      `<select id="p_ref" onchange="FIP.set('process.ref',this.value)">` +
      DATA.batches.map(function (x) { return `<option value="${x.batch_id}"${x.batch_id === refBid ? " selected" : ""}>${x.batch_id} · ${x.success ? "成功" : "失败"}</option>`; }).join("") +
      `</select></div><div id="p_dev"></div>`);
    // 实验设计（本次实验参数输入）——位于 What-if 上方，作为仿真基准
    html += card(`<div class="sub">🧪 本次实验参数输入（实验设计）· What-if 仿真基准</div><div id="p_exp"></div>` +
      `<div class="method-note">录入本次发酵实验的设计参数（宿主与分子设计、发酵罐规模、起始装液量，各发酵阶段的 DO/pH/温度/补料，以及培养基 C:N 与配方）。What-if 仿真模块的 DO、温度、补料速率、发酵时长可一键「采用本模块设定」作为仿真基准（见下方 What-if 卡片按钮）。可保存参数 / 生成 / 下载实验方案 JSON。</div>`);
    html += card(`<div class="sub">🔧 What-if 工艺仿真（简化情景模型 · 基准取自上方实验设计）</div>` +
      `<div style="display:flex;justify-content:space-between;align-items:center;gap:10px;margin-bottom:8px;flex-wrap:wrap">` +
      `<div class="ysub" style="margin:0">调参后实时仿真；或一键采纳上方「本次实验参数输入」模块设定作为基准</div>` +
      `<button class="btn2" onclick="FIP.ptAdoptExp()">↰ 采用上方实验设计模块设定</button>` +
      `</div>` +
      `<div class="whatif-grid">` +
      ptSliderRow("DO 设定 (%)", "process.do", state.process.do, 10, 60, 1) +
      ptSliderRow("温度 (℃)", "process.temp", state.process.temp, 20, 37, 0.5) +
      ptSliderRow("补料速率 (×)", "process.feed", state.process.feed, 0.5, 2.0, 0.1) +
      ptDurRow() +
      `</div><div id="p_whatif" style="margin-top:12px"></div>` +
      ptWhatIfParamChips() +
      `<div class="method-note">情景模型：DO/温度/补料经响应面拆分为<b>生物量因子 modX</b> 与<b>比产物生成速率因子 modQ</b>，效价 P=∫q_P·X dt ⇒ Δtiter ≈ modX·modQ（各计一次；旧实现的 DCW 与 titer 双乘同一因子 ⇒ mod² 错误放大已修正）。其中 DO 对生长呈饱和、对表达更敏感；温度按生长 34 ℃ / 表达 30 ℃ 双高斯；补料 X∝F<sup>0.65</sup>、q_P∝F<sup>−0.15</sup>（过量补料→溢流/乙酸使比产率下降）。生物量按宿主上限钳制（E. coli 200 / 毕赤酵母 400 g DCW/L）。情景参数轨迹（基线 vs 情景 至设定发酵时长）已合并为一张多轴图（勾选参数即时刷新）。<b>注：</b>为趋势演示级经验响应面，非机理预测；正式终点预测见下方 RF 模型。</div>`);
    if (pred && b) {
      var t = pred.titer_g_l;
      html += card(`<div class="sub">🎯 批次终点预测（RF · 截至 60% 进度）</div>` +
        `<table class="tbl"><tr><th>指标</th><th>预测均值</th><th>95% 置信区间</th></tr>` +
        `<tr><td>效价 titer (g/L)</td><td>${fmt(t.mean, 2)}</td><td>[${fmt(t.ci_low, 2)}, ${fmt(t.ci_high, 2)}]</td></tr>` +
        (pred.harvest_dcw_g_l ? `<tr><td>DCW (g/L)</td><td>${fmt(pred.harvest_dcw_g_l.mean, 2)}</td><td>[${fmt(pred.harvest_dcw_g_l.ci_low, 2)}, ${fmt(pred.harvest_dcw_g_l.ci_high, 2)}]</td></tr>` : "") +
        `</table>` + ptModelQualityRow("titer_g_l", "效价") + ptModelQualityRow("harvest_dcw_g_l", "DCW") +
        ptCICalibrationRow("titer_g_l") +
        `<div class="hint">实测效价 ${fmt(b.titer_g_l, 2)} g/L。软测量与方法见《5.2 模型结构》；终点预测对照《4.3 数据量门槛》（启动 50 批，本平台 ${DATA.stats.n_batches} 批）。</div>`);
    }
    $("p_out").innerHTML = html;
    renderExpForm();
    if (ts) {
      var hrs = ts.hours;
      var sel = PT_PARAMS.filter(function (p) { return state.process.params[p.key] && ts[p.key]; });
      if (!sel.length) sel = PT_PARAMS.filter(function (p) { return ts[p.key]; }).slice(0, 1);
      var multiSeries = sel.map(function (p) {
        var ys = ts[p.key];
        var lo = Math.min.apply(null, ys), hi = Math.max.apply(null, ys);
        return { name: p.label, color: p.color, unit: p.unit, ylo: lo, yhi: hi, points: hrs.map(function (h, i) { return [h, ys[i]]; }) };
      });
      $("p_lines").innerHTML = '<div id="p_multi"></div><div style="color:#7d93a6;font-size:12px;margin-top:4px">已选 ' + sel.length + ' 项参数 · 每个参数独立纵坐标（真实量纲，左右分列），曲线按各自量程映射于同一时间轴</div>';
      lineChartMulti($("p_multi"), multiSeries, { phases: ts.phases });
      lineChart($("p_soft"), [
        { name: "OUR反推", color: PRIMARY, points: hrs.map(function (h, i) { return [h, doOur[i]]; }), unit: "g/L" },
        { name: "CER反推", color: "#ffce4d", points: hrs.map(function (h, i) { return [h, doCer[i]]; }), unit: "g/L" },
        { name: "Logistic 拟合", color: ACCENT, points: hrs.map(function (h, i) { return [h, doExp[i]]; }), unit: "g/L" },
        { name: "集成软测量", color: "#80ed99", points: hrs.map(function (h, i) { return [h, ts.DCW[i]]; }), unit: "g/L" }
      ], { ymin: 0, ytitle: "DCW (g/L)" });
      var last = hrs.length - 1;
      if (last >= 0) {
        var mx = Math.max(doOur[last], doCer[last], doExp[last]), mn = Math.min(doOur[last], doCer[last], doExp[last]);
        var spread = (mx - mn) / Math.max(0.01, (doOur[last] + doCer[last] + doExp[last]) / 3);
        var cons = spread < 0.15 ? "三法高度一致 → 软测量高置信" : (spread < 0.35 ? "三法中度分歧，建议核查信号质量" : "三法显著分歧，软测量置信低");
        var calibTxt = SM.calibrated
          ? ("离线 harvest DCW " + fmt(SM.anchor, 2) + " g/L 标定（先验 " + fmt(SM.raw_end, 1) + " → 校正 " + fmt(SM.calib, 2) + "×；" +
            (SM.mode === "yield"
              ? "模式：标定生长氧耗 qO₂_g=" + fmt(SM.qo2_g, 1) + " mmol/gX，Y_OX=" + fmt(SM.y_ox, 4) + " g/mmol O₂"
              : "模式：氧耗能量不平衡，按 ODE 形状缩放至实测终点（非得率标定）") + "）")
          : "未标定（采用先验 Y_OX=" + fmt(SM.y_ox, 4) + " g/mmol O₂）";
        $("p_softinfo").innerHTML = '<div class="method-note">末值：OUR 反推 ' + fmt(doOur[last], 1) + ' · CER 反推 ' + fmt(doCer[last], 1) +
          ' · Logistic 拟合 ' + fmt(doExp[last], 1) + ' g/L DCW（集成软测量 ' + fmt(ts.DCW[last], 1) + '）。' + cons +
          '。<br><b>三法：</b>① OUR 氧衡算 ② CER 碳衡算（按 RQ=' + fmt(SS_RQ, 2) + ' 折算等效耗氧）③ Logistic 生长律拟合（不依赖尾气，仅校核形状）。' +
          'qO₂_g=' + fmt(SM.qo2_g, 2) + ' mmol/gX、mO₂=' + fmt(SM.mo2, 2) + ' mmol/gX/h、X₀=' + fmt(SM.x0, 2) + ' g/L。' +
          '<b>标定：</b>' + calibTxt +
          (SM.diag ? '<br><b>⚠ 诊断：</b>' + esc(SM.diag) : '') +
          '。<br><b>注：</b>旧版按「X=OUR/常数」逐点反演（隐含恒定比摄氧率），终点偏离实测 −61%~+158%，已废止。</div>';
      }
      $("p_kpi").innerHTML = ptKpiCards(kpi, b, ts, bid) +
        `<div class="kpi-grid" style="margin-top:10px"><div><div class="ksub" style="margin-bottom:4px">RQ 轨迹</div><div id="p_rq"></div></div><div><div class="ksub" style="margin-bottom:4px">μ 轨迹 1/h</div><div id="p_mu"></div></div></div>`;
      lineChart($("p_rq"), [{ name: "RQ", color: "#ffce4d", points: hrs.map(function (h, i) { return [h, kpi.rq[i]]; }) }], { ytitle: "RQ", ymin: 0 });
      lineChart($("p_mu"), [{ name: "μ", color: PRIMARY, points: hrs.map(function (h, i) { return [h, kpi.mu[i]]; }) }], { ytitle: "μ (1/h)" });
    }
    if (dev) $("p_dev").innerHTML = ptDevTable(dev);
    if (wi && ts && ts.DCW && ts.DCW.length) {
      var traj = ptWhatIfTraj(state.process, base, ts, bid);
      var wiHtml = '<div style="display:flex;gap:12px;flex-wrap:wrap">' +
        '<div style="flex:1;min-width:300px"><div class="ksub" style="margin-bottom:3px">DCW 趋势：基线 vs 情景</div><div id="p_wi_dcw"></div></div>' +
        '<div style="flex:1;min-width:300px"><div class="ksub" style="margin-bottom:3px">效价累积（情景）</div><div id="p_wi_titer"></div></div></div>';
      wiHtml += '<div class="predbox' + (traj.rel < -15 ? " risk" : "") + '" style="margin-top:12px"><div class="pv">' + fmt(traj.titer, 2) + ' <small>g/L</small></div>' +
        '<div class="pl">预测终点效价（基线 ' + fmt(base.titer, 2) + ' g/L，' + (traj.rel >= 0 ? "+" : "") + fmt(traj.rel, 0) + '%）</div>' +
        '<div class="pv" style="font-size:18px;margin-top:8px">' + fmt(traj.endDCW, 1) + ' <small>g DCW/L</small></div>' +
        '<div class="pl">情景终点生物量（基线 ' + fmt(traj.dcwBase[traj.dcwBase.length - 1], 1) + ' g/L）</div>' +
        '<div class="pv" style="font-size:18px;margin-top:8px">' + fmt(traj.dur, 0) + ' <small>h</small></div>' +
        '<div class="pl">设定发酵时长（基线 ' + fmt(base.dur, 0) + ' h' + (traj.dur > traj.lastH ? ' · 含外推区间' : '') + '）</div>' +
        '<div class="hint" style="margin-top:6px">情景因子：生物量 X ×' + fmt(traj.fa.modX, 2) + ' · 比产率 q_P ×' + fmt(traj.fa.modQ, 2) +
        ' ⇒ titer ∝ X·q_P（各计一次）</div>' +
        (traj.capped ? '<div class="hint" style="margin-top:4px;color:#ff6b6b">⚠ 情景生物量已达' + esc(ptHostOf(bid) === "yeast" ? "毕赤酵母" : "该宿主") +
          '可达上限 ' + fmt(traj.xMax, 0) + ' g DCW/L（氧传递/流变限制），已钳制，实际不可超越</div>' : '') +
        '</div>';
      $("p_whatif").innerHTML = wiHtml;
      lineChart($("p_wi_dcw"), [
        { name: "DCW 基线", color: "#6b7d8f", points: hrs.map(function (h, i) { return [h, traj.dcwBase[i]]; }) },
        { name: "DCW 情景", color: PRIMARY, points: traj.dcwScen }
      ], { ymin: 0, ytitle: "DCW (g/L)", phases: ts.phases, vline: { x: traj.dur, color: "#ff6b6b", label: "设定时长 " + fmt(traj.dur, 0) + "h" } });
      lineChart($("p_wi_titer"), [
        { name: "效价累积", color: "#80ed99", points: traj.titerScen }
      ], { ymin: 0, ytitle: "效价 (g/L)", phases: ts.phases, vline: { x: traj.dur, color: "#ff6b6b", label: "设定时长 " + fmt(traj.dur, 0) + "h" } });
      // 情景参数轨迹（可选）：合并为一张多轴图（每个参数独立纵坐标，基线 vs 情景，至设定发酵时长）
      var wiParams = ["TITER"].concat(PT_PARAMS.map(function (p) { return p.key; })).filter(function (k) { return state.process.wipar[k]; });
      renderWhatIfScenarioChart($("p_whatif_params"), wiParams, base, ts, traj);
    }
    renderCompare();
  }

  // ---- Expression Twin (参考复制 YeastExpress Pro) ----
  var EXAMPLES = {
    hsa: "DAHKSEVAHRFKDLGEENFKALVLIAFAQYLQQCPFEDHVKLVNEVTEFAKTCVADESAENCDKSLHTLFGDKLCTVATLRETYGEMADCCAKQEPERNECFLSHKDDSPDLPKLKPDPNTLCDEFKADEKKFWGKYLYEIARRHPYFYAPELLYYANKYNGVFQECCQAEDKGACLLPKIETMREKVLASSARQRLRCASIQKFGERALKAWSVARLSQKFPKAEFVEVTKLVTDLTKVHKECCHGDLLECADDRADLAKYICDNQDTISSKLKECCDKPLLEKSHCIAEVEKDAIPENLPPLTADFAEDKDVCKNYQEAKDAFLGSFLYEYSRRHPEYAVSVLLRLAKEYEATLEECCAKDDPHACYSTVFDKLKHLVDEPQNLIKQNCDQFEKLGEYGFQNALIVRYTRKVPQVSTPTLVEVSRSLGKVGTRCCTKPESERMPCTEDYLSLILNRLCVLHEKTPVSEKVTKCCTESLVNRRPCFSALTPDETYVPKAFDEKLFTFHADICTLPDTEKQIKKQTALVELVKHKPKATEEQLKTVMENFVAFVDKCCAADDKEACFAVEGPKLVVSTQTALA",
    gfp: "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK",
    insulin: "MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTGIVEQCCTSICSLYQLENYCN"
  };
  var PRED_LEVELS = [["none", "不表达"], ["low", "低表达"], ["mid", "中等表达"], ["high", "高表达"]];

  function layerCard(L, name, score) {
    var pct = Math.round(Math.max(0, Math.min(1, score)) * 100);
    var col = score >= 0.6 ? PRIMARY : (score >= 0.38 ? "#ffce4d" : "#ff5a5a");
    return '<div class="layer-card"><div class="ln">' + L + " " + esc(name) + '</div><div class="layer-bar"><div class="layer-fill" style="width:' + pct + '%;background:' + col + '"></div></div><div class="layer-score" style="color:' + col + '">' + fmt(score, 2) + '</div></div>';
  }
  function renderPred(level, yld) {
    var rows = PRED_LEVELS.map(function (p) {
      var yv = p[0] === level ? (p[0] === "none" ? "0" : fmt(yld, 0)) : "--";
      var cls = p[0] === level ? " pred-active" : "";
      return '<div class="predrow' + cls + '"><span>' + p[1] + '</span><span class="pv">' + yv + ' mg/L</span></div>';
    }).join("");
    var box = $("e_pred"); if (!box) return;
    box.innerHTML = rows + '<div class="ysub" style="margin-top:8px">输入蛋白序列并分析后，将显示预测的表达等级和产量估算</div>';
  }

  function renderExpressionYeast(body) {
    var e = state.expr;
    var hosts = [["P. pastoris GS115", "Pichia_pastoris_GS115"], ["P. pastoris X33", "Pichia_pastoris_X33"]];
    var promoters = [["PAOX1", "PAOX1"], ["PGAP", "PGAP"], ["GAP", "GAP"], ["TEF1", "TEF1"], ["GAL1", "GAL1"], ["T7lac", "T7lac"], ["tac", "tac"], ["trc", "trc"]];
    var sps = [["alpha-MF", "alpha-MF"], ["alpha-factor", "alpha-factor"], ["Ost1", "Ost1"], ["PHO5-SP", "PHO5-SP"], ["GlaA-SP", "GlaA-SP"], ["SUC2-SP", "SUC2-SP"], ["Native SP", "Native SP"], ["none", "none"]];
    var levelOpts = [["0 - 不表达", "none"], ["1 - 低表达", "low"], ["2 - 中等表达", "mid"], ["3 - 高表达", "high"]];
    function opt(arr, cur) { return arr.map(function (o) { return '<option value="' + o[1] + '"' + (cur === o[1] ? " selected" : "") + '>' + o[0] + '</option>'; }).join(""); }
    var html = "";
    // Hero
    html += '<div class="yhero"><div class="ytitle">YeastExpress Pro</div><div class="ysub">酵母重组蛋白表达预测评分系统</div>' +
      '<div class="ytrans">从 Transcription 到 Secretion 的全流程数字化评估，覆盖 7 个关键层级，整合 ExpressYeaself、AlphaFold、SignalP、NetNGlyc 等权威工具算法</div></div>';
    // Section 1 蛋白序列输入
    html += card(
      '<div class="sub">蛋白序列输入</div><div class="ysub" style="margin-bottom:10px">输入氨基酸序列，启动7层级全流程评估</div>' +
      '<div class="formgrid">' +
      '<label>蛋白名称</label><input id="e_name" value="' + esc(e.protein_name || "") + '" placeholder="如 GFP">' +
      '<label>宿主菌株</label><select id="e_host" onchange="FIP.calcExpr()">' + opt(hosts, e.host) + '</select>' +
      '<label>信号肽选择</label><select id="e_sp" onchange="FIP.calcExpr()">' + opt(sps, e.sp) + '</select>' +
      '</div>' +
      '<div class="ksub" style="margin:10px 0 2px">氨基酸序列 (Single-letter code)</div>' +
      '<textarea id="e_seq" rows="3" style="width:100%;font-family:monospace" placeholder="MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKT" oninput="if($(\'e_count\'))$(\'e_count\').textContent=(this.value.length)+\' aa\'">' + esc(e.seq || "") + '</textarea>' +
      '<div class="count"><span id="e_count">' + (e.seq ? e.seq.length : 0) + ' aa</span> · 示例：' +
      '<span class="exbtn" onclick="FIP.example(\'hsa\')">HSA</span> ' +
      '<span class="exbtn" onclick="FIP.example(\'gfp\')">GFP</span> ' +
      '<span class="exbtn" onclick="FIP.example(\'insulin\')">Insulin</span></div>' +
      '<div class="sub" style="margin-top:10px">转录层参数配置</div>' +
      '<div class="formgrid">' +
      '<label>启动子</label><select id="e_promoter" onchange="FIP.calcExpr()">' + opt(promoters, e.promoter) + '</select>' +
      '<label>CAI</label><input id="e_cai" value="' + esc(e.cai || "") + '" style="width:90px" onchange="FIP.calcExpr()">' +
      '<label>GC%</label><input id="e_gc" value="' + esc(e.gc || "") + '" style="width:90px" onchange="FIP.calcExpr()">' +
      '</div>' +
      '<button class="btn" onclick="FIP.calcExpr()">🚀 启动全流程评估</button>'
    );
    // Section 2+3 Pipeline + 表达量预测
    html += '<div style="display:flex;gap:14px;flex-wrap:wrap">';
    html += card('<div class="sub">评估流程 Pipeline</div><div class="ysub" style="margin-bottom:8px">7层级从基因到分泌蛋白的完整表达路径</div><div id="e_pipeline"></div>', 'style="flex:1.25;min-width:380px"');
    html += card('<div class="sub">表达量预测</div><div class="ysub" style="margin-bottom:8px">分析后自动预测蛋白表达水平和产量</div><div id="e_pred"></div>', 'style="flex:1;min-width:300px"');
    html += "</div>";
    // Section 4 实验数据回流飞轮
    html += card(
      '<div style="display:flex;justify-content:space-between;align-items:center"><div class="sub" style="margin:0">实验数据回流飞轮</div><button class="btn2" onclick="FIP.clearFeedback()">清空</button></div>' +
      '<div class="ysub" style="margin-bottom:8px">回填实验结果，校准预测模型（云端存储）</div>' +
      '<div class="formgrid">' +
      '<label>表达等级</label><select id="e_fblevel">' + opt(levelOpts, "low") + '</select>' +
      '<label>溶解度 (%)</label><input id="e_sol" placeholder="0-100">' +
      '<label>包涵体 (%)</label><input id="e_inc" placeholder="0-100">' +
      '<label>产量 (mg/L)</label><input id="e_fbyield" placeholder="如 1200">' +
      '<label>培养方式 ⚠</label><select id="e_fbmode">' +
      opt([["补料分批（高密）· 可并入校准", "fed_batch"], ["摇瓶 · 仅存档不并入", "shake_flask"],
           ["分批 · 仅存档", "batch"], ["连续 · 仅存档", "continuous"],
           ["未声明 · 不并入", "unspecified"]], "fed_batch") + '</select>' +
      '</div>' +
      '<div class="hint" style="margin-top:6px">⚠ 培养方式决定产量量级：摇瓶效价通常比高密低 <b>1–2 个数量级</b>，' +
      '直接混入 g/L 训练集会把模型拉偏。只有「补料分批（高密）」回填计入校准阈值并用于再校准。</div>' +
      '<label style="display:block;margin-top:8px">备注 (可选)</label><textarea id="e_note" rows="2" placeholder="实验条件、培养基、诱导策略…"></textarea>' +
      '<button class="btn" id="e_save" disabled onclick="FIP.saveFeedback()">保存实验结果</button>' +
      '<div class="hint" id="e_calib"></div>' +
      '<div class="sub" style="margin-top:10px">已保存的反馈数据 · <span id="e_fbcount"></span></div>' +
      '<div id="e_fblist" class="fb-list"></div>'
    );
    // 补充：FIP 高表达构架排序 (PEFM)
    html += card('<div class="sub">📊 高表达构架排序 (FIP · PEFM)</div><div style="max-height:460px;overflow:auto"><table class="tbl" id="e_rank">' +
      '<tr><th>蛋白</th><th>宿主</th><th>级别</th><th>效价(mg/L)</th><th>证据</th></tr></table></div>' +
      '<div class="hint">证据等级 high/medium 来自酵母文献实证集（evidence_level / source_ref），可溯源；排序由 Expression Twin 7 层评分生成。</div>', 'style="margin-top:0"');
    // 实验设计信息量反馈 (DOE)：设计期构型覆盖 + 批次/工艺设定点覆盖（show_construct=true，对应 Streamlit Expression 页）
    html += card(doeSectionHtml({ host: "Pichia_pastoris_X33", showConstruct: true, scale: null, withHostSelect: true }));
    body.innerHTML = html;
    renderExprRank();
    renderFeedbackList();
    calcExpr();
  }

  function calcExpr() {
    var seq = $("e_seq") ? $("e_seq").value : "";
    var name = $("e_name") ? $("e_name").value : "";
    var host = $("e_host") ? $("e_host").value : "Pichia_pastoris_X33";
    var promoter = $("e_promoter") ? $("e_promoter").value : "PAOX1";
    var sp = $("e_sp") ? $("e_sp").value : "alpha_factor";
    var cai = $("e_cai") ? $("e_cai").value : "";
    var gc = $("e_gc") ? $("e_gc").value : "";
    state.expr = { protein_name: name, host: host, promoter: promoter, sp: sp, seq: seq, cai: cai, gc: gc };
    var r = sevenLayer(host, promoter, sp, seq, cai, gc, null);
    state.expr._result = r; state.expr.analyzed = true;
    if ($("e_count")) $("e_count").textContent = (seq ? seq.length : 0) + " aa";
    // Pipeline 7 卡
    var cards = r.levels.map(function (l, i) { return layerCard("L" + (i + 1), l.name, l.score); });
    cards.push(layerCard("L7", "全流程预测", r.overall));
    if ($("e_pipeline")) $("e_pipeline").innerHTML = '<div class="pipeline">' + cards.join("") + "</div>";
    renderPred(r.expression_level, r.yield_mg_l);
    var sb = $("e_save"); if (sb) sb.disabled = false;
  }

  function renderExprRank() {
    var t = $("e_rank"); if (!t) return;
    var ranked = DATA.constructs.slice().sort(function (a, b) { return (b.overall || 0) - (a.overall || 0); }).slice(0, 40);
    t.innerHTML = '<tr><th>蛋白</th><th>宿主</th><th>级别</th><th>效价(mg/L)</th><th>证据</th></tr>' +
      ranked.map(function (c) {
        var lvl = { high: "🔥", mid: "🟢", low: "🟡", none: "⚪" }[c.expression_level] || "—";
        return "<tr><td>" + esc(c.protein_name) + "</td><td style='color:" + hostColor(c.host) + "'>" + (c.host || "").replace("_", " ") + "</td><td>" + lvl + "</td><td>" + fmt(c.yield_mg_l, 0) + "</td><td>" + (c.evidence_level || "—") + (c.source_ref ? " *" : "") + "</td></tr>";
      }).join("");
  }

  // 统一存储抽象 fipStore：优先 localStorage，沙箱(无 allow-same-origin)下透明降级为内存存储
  // ——避免预览面板 iframe 中访问 localStorage 抛 SecurityError 导致载入/保存/库操作失败
  var _fipMem = {};
  var _fipLsOK = (function () {
    try { var k = "__fip_ls_probe__"; localStorage.setItem(k, "1"); localStorage.removeItem(k); return true; }
    catch (e) { return false; }
  })();
  function fipGet(k) {
    try { if (_fipLsOK) return localStorage.getItem(k); } catch (e) {}
    return (Object.prototype.hasOwnProperty.call(_fipMem, k)) ? _fipMem[k] : null;
  }
  function fipSet(k, v) {
    try { if (_fipLsOK) { localStorage.setItem(k, v); return; } } catch (e) {}
    _fipMem[k] = v;
  }
  function fipDel(k) {
    try { if (_fipLsOK) { localStorage.removeItem(k); return; } } catch (e) {}
    delete _fipMem[k];
  }
  function fipStorageMode() { return _fipLsOK ? "local" : "memory"; } // 供 UI 提示

  // 实验反馈飞轮 (localStorage 持久化，复刻"云端存储"语义)
  function fbKey() { return "yep_feedback_v1"; }
  function getFeedback() { try { return JSON.parse(fipGet(fbKey()) || "[]"); } catch (e) { return []; } }
  function exampleFill(key) {
    var s = EXAMPLES[key] || "";
    var ta = $("e_seq"); if (ta) { ta.value = s; if ($("e_count")) $("e_count").textContent = s.length + " aa"; }
    if (key === "hsa" && $("e_name")) $("e_name").value = "HSA";
    if (key === "gfp" && $("e_name")) $("e_name").value = "GFP";
    if (key === "insulin" && $("e_name")) $("e_name").value = "Insulin";
    FIP.calcExpr();
    if (FIP.yeastRecompute) FIP.yeastRecompute();
  }
  function saveFeedback() {
    if (!state.expr.analyzed) { alert("请先输入蛋白序列并启动全流程评估"); return; }
    var name = $("e_name") ? $("e_name").value : "";
    var level = $("e_fblevel") ? $("e_fblevel").value : "low";
    var sol = $("e_sol") ? $("e_sol").value : "";
    var inc = $("e_inc") ? $("e_inc").value : "";
    var yld = $("e_fbyield") ? $("e_fbyield").value : "";
    var note = $("e_note") ? $("e_note").value : "";
    var mode = $("e_fbmode") ? $("e_fbmode").value : "unspecified";
    if (yld === "") { alert("请填写产量 (mg/L)"); return; }
    var yv = Number(yld);
    if (!isFinite(yv) || yv < 0) { alert("产量必须是非负数值 (mg/L)"); return; }
    // 量级合理性提示（不阻断保存，仅提醒复核单位）
    var CAP = { shake_flask: 5000, batch: 20000, fed_batch: 200000, continuous: 200000, unspecified: 200000 };
    if (CAP[mode] && yv > CAP[mode]) {
      if (!confirm(yv + " mg/L 超出「" + mode + "」的常规量级（≤" + CAP[mode] + " mg/L）。\n" +
        "摇瓶常为 mg/L 级、高密常为 g/L 级，请确认单位无误。\n仍要保存吗？")) return;
    }
    var arr = getFeedback();
    arr.push({ id: Date.now(), protein: name || "(未命名)", level: level, sol: Number(sol) || null, inc: Number(inc) || null, yld: yv, culture_mode: mode, note: note, date: new Date().toISOString().slice(0, 10) });
    fipSet(fbKey(), JSON.stringify(arr));
    renderFeedbackList(); updateCalibHint();
  }
  function clearFeedback() { fipDel(fbKey()); renderFeedbackList(); updateCalibHint(); }
  function renderFeedbackList() {
    var arr = getFeedback();
    var c = $("e_fbcount"); if (c) c.textContent = arr.length + " 条记录";
    var box = $("e_fblist"); if (!box) return;
    if (!arr.length) { box.innerHTML = '<div class="ysub">暂无反馈数据</div>'; return; }
    var MODE_TXT = { fed_batch: "高密", shake_flask: "摇瓶", batch: "分批", continuous: "连续", unspecified: "未声明" };
    box.innerHTML = arr.slice().reverse().map(function (f) {
      var lvlTxt = { none: "不表达", low: "低表达", mid: "中等表达", high: "高表达" }[f.level] || f.level;
      var pname = (f.protein || "").length > 14 ? f.protein.slice(0, 14) + "..." : f.protein;
      var m = f.culture_mode || "unspecified";
      var okM = (m === "fed_batch");
      return '<div class="fb-item"><span class="fb-i">' + esc(pname) + '</span><span class="fb-lvl">' + lvlTxt + '</span>' +
        '<span class="fb-y">' + fmt(f.yld, 0) + ' mg/L</span>' +
        '<span class="fb-lvl" style="color:' + (okM ? "#1f9e89" : "#c05621") + '">' + (MODE_TXT[m] || m) + (okM ? " ✅" : " ⚠") + '</span>' +
        '<span class="fb-d">' + f.date + '</span></div>';
    }).join("");
  }
  function updateCalibHint() {
    var arr = getFeedback();
    var h = $("e_calib"); if (!h) return;
    var mergeable = arr.filter(function (f) { return (f.culture_mode || "unspecified") === "fed_batch"; });
    var others = arr.length - mergeable.length;
    if (mergeable.length < 3) {
      h.innerHTML = "校准状态：需要 <b>3 条高密同口径</b>回填（当前 " + mergeable.length + " 条"
        + (others ? "，另有 " + others + " 条非高密/未声明，<b>不并入 g/L 训练集</b>" : "") + "）。"
        + "摇瓶效价通常比高密低 1–2 个数量级，混入会把模型拉偏。数据存储在云端，多设备可同步。";
    } else {
      h.innerHTML = "校准状态：已具备 " + mergeable.length + " 条高密同口径数据，模型校准已启用（本机回放）。"
        + (others ? " 另有 " + others + " 条非高密回填仅存档、不参与校准。" : "")
        + "后续预测将向您的实际表达分布对齐。";
    }
  }

  // =========================================================================
  // 酵母 Expression Twin · L1-L7 独立详情页（严格参考 yeast-express-pro skill）
  // 移植 scoring-template.ts / scoring-algorithms.md 的 7 层级启发式算法
  // =========================================================================
  var YEAST_LEVELS = [
    { L: "L1", name: "转录层", en: "Transcription", w: 0.15, eval: "启动子强度、基因拷贝数、整合位点、mRNA稳定性", tools: ["ExpressYeaself", "Promoter Calculator"] },
    { L: "L2", name: "翻译层", en: "Translation", w: 0.15, eval: "CAI、GC含量、5'端发夹结构、稀有密码子", tools: ["GenSmart Codon Optimization", "IDT Codon Optimization", "JCat"] },
    { L: "L3", name: "蛋白折叠", en: "Folding", w: 0.15, eval: "pLDDT、无序区、聚集热点、可溶性", tools: ["AlphaFold", "ColabFold", "Camsol"] },
    { L: "L4", name: "ER加工", en: "ER Processing", w: 0.15, eval: "二硫键、分子量、糖基化负荷、UPR风险", tools: ["DiANNA", "UPR Predictor"] },
    { L: "L5", name: "糖基化", en: "Glycosylation", w: 0.10, eval: "N/O-糖基化位点、超糖基化风险", tools: ["NetNGlyc", "NetOGlyc", "GlycoEP"] },
    { L: "L6", name: "分泌", en: "Secretion", w: 0.20, eval: "信号肽、跨膜区、亚细胞定位", tools: ["SignalP 6.0", "Phobius", "WoLF PSORT"] },
    { L: "L7", name: "全流程预测", en: "Full Pipeline", w: 0.10, eval: "分泌倾向、表达量、可溶性（DeepLoc/DeepSec/ESM2）", tools: ["DeepLoc 2.0", "DeepSec", "ProtTrans", "ESM2"] }
  ];
  var YEAST_SP_OPTIONS = [["alpha-MF", "alpha-MF"], ["alpha-factor", "alpha-factor"], ["Ost1", "Ost1"], ["PHO5-SP", "PHO5-SP"], ["GlaA-SP", "GlaA-SP"], ["SUC2-SP", "SUC2-SP"], ["Native SP", "Native SP"], ["none", "none"]];
  var YEAST_PROMOTERS = [["PAOX1", "PAOX1"], ["PGAP", "PGAP"], ["GAP", "GAP"], ["TEF1", "TEF1"], ["GAL1", "GAL1"], ["T7lac", "T7lac"], ["tac", "tac"], ["trc", "trc"]];

  function yp_grade(s) { if (s >= 90) return "A"; if (s >= 75) return "B"; if (s >= 60) return "C"; if (s >= 40) return "D"; return "F"; }
  function yp_status(v, th) { if (v >= th[1]) return "good"; if (v >= th[0]) return "warning"; return "danger"; }
  function yp_clean(seq) { return (seq || "").replace(/[^ACDEFGHIKLMNPQRSTVWY]/gi, "").toUpperCase(); }
  function yp_cnt(seq, aa) { var n = 0; for (var i = 0; i < seq.length; i++) if (seq[i] === aa) n++; return n; }
  function yp_mw(seq) { var m = { A: 89.09, R: 174.20, N: 132.12, D: 133.10, C: 121.16, E: 147.13, Q: 146.15, G: 75.03, H: 155.16, I: 131.17, L: 131.17, K: 146.19, M: 149.21, F: 165.19, P: 115.13, S: 105.09, T: 119.12, W: 204.23, Y: 181.19, V: 117.15 }; var t = 0; for (var i = 0; i < seq.length; i++) t += (m[seq[i]] || 0); return (t - (seq.length - 1) * 18.02) / 1000; }
  function yp_oSites(seq) { var c = 0; for (var i = 0; i < seq.length; i++) { if (seq[i] === "S" || seq[i] === "T") { var ctx = seq.slice(Math.max(0, i - 5), i + 6); var k = 0; for (var j = 0; j < ctx.length; j++) if ("PST".indexOf(ctx[j]) >= 0) k++; if (k >= 3) c++; } } return c; }

  // 序列驱动自动预估（确定性，无随机；移植自 autoPredictFromSequence）
  function yp_auto(seq, spKey) {
    var clean = yp_clean(seq); var len = clean.length || 1;
    var spInfo = SIGNAL_PEPTIDE_SCORE[spKey]; var hasSel = spInfo && spInfo > 0 && spKey !== "none";
    var gc = (clean.match(/[GC]/g) || []).length / len;
    var cai = Math.max(0.2, Math.min(0.95, 0.5 + (1 - Math.abs(gc - 0.45)) * 0.6 + (len > 200 ? 0.1 : 0)));
    var head = clean.slice(0, 30); var hg = 0; for (var i = 0; i < head.length; i++) hg += (KD[head[i]] || 0); hg = hg / Math.max(head.length, 1);
    var dG5 = -(0.5 + 3 * Math.max(0, hg));
    var rare = Math.floor((1 - cai) * 30);
    var cys = yp_cnt(clean, "C"), pro = yp_cnt(clean, "P"), gly = yp_cnt(clean, "G");
    var disProp = (pro + gly) / len;
    var plddt = Math.max(40, Math.min(95, 80 - disProp * 100 + (len < 300 ? 5 : -5)));
    var nSites = nGlyco(clean);
    var mw = yp_mw(clean);
    var hydroSet = "AILMFWV".split(""); var hc = 0; for (var i = 0; i < clean.length; i++) if (hydroSet.indexOf(clean[i]) >= 0) hc++;
    var hydroRatio = hc / len;
    var agg = Math.min(100, hydroRatio * 200);
    var sol = 1 - agg / 100;
    var disBonds = Math.floor(cys / 2);
    var upr = Math.min(100, (cys > 6 ? 30 : 0) + (mw > 80 ? 25 : 0) + (nSites.length > 5 ? 25 : 0) + (len > 500 ? 20 : 0));
    var nTerm = clean.slice(0, 25); var ntH = 0; for (var i = 0; i < nTerm.length; i++) if (hydroSet.indexOf(nTerm[i]) >= 0) ntH++;
    var nTermHydro = ntH / Math.max(nTerm.length, 1);
    var nativeHasSP = nTermHydro > 0.45; var nativeProb = Math.min(0.99, nTermHydro * 1.5);
    var hasSP = hasSel ? true : nativeHasSP;
    var spProb = hasSel ? Math.min(0.99, spInfo * 0.95) : nativeProb;
    var secProb = hasSP ? 0.65 + spProb * 0.2 : 0.15;
    var cytoProb = hasSP ? 0.10 : 0.55;
    var deepSec = hasSP ? 0.70 : 0.25;
    var pTexpr = Math.max(0.2, Math.min(0.9, 0.5 + (cai - 0.5) * 0.8 + (plddt - 70) * 0.005));
    var pTsol = Math.max(0.2, Math.min(0.9, sol));
    var pTsec = hasSP ? 0.65 : 0.20;
    return {
      cai: cai, gc: gc, dG5: dG5, rare: rare, len: len,
      plddt: plddt, disorderRegions: Math.floor(disProp * 10), agg: agg, sol: sol,
      cys: cys, mw: mw, nSites: nSites.length, disBonds: disBonds, upr: upr,
      hasSP: hasSP, spProb: spProb, secProb: secProb, cytoProb: cytoProb,
      deepSec: deepSec, pTexpr: pTexpr, pTsol: pTsol, pTsec: pTsec,
      loc: { Secreted: secProb, Cytoplasm: cytoProb, Membrane: hasSP ? 0.08 : 0.15, Vacuole: 0.05, Nucleus: 0.03, Mitochondria: 0.04 }
    };
  }

  function yp_L1(promoter, auto) {
    var ps = { GAP: 85, AOX1: 95, TEF1: 80, PHO5: 70, CUP1: 75, custom: 50, PGAP: 85, GAL1: 70, T7lac: 50, tac: 50, trc: 50 };
    var promoterScore = ps[promoter] || 50;
    var copy = 2, copyScore = Math.min(100, copy * 20);
    var integ = "AOX1_locus"; var is = { AOX1_locus: 90, HIS4: 75, random: 45, episomal: 60 }; var integScore = is[integ] || 50;
    var mrna = 70, term = 70;
    var weighted = promoterScore * 0.30 + copyScore * 0.20 + integScore * 0.20 + mrna * 0.15 + term * 0.15;
    var sub = [
      { name: "启动子强度", value: promoterScore, max: 100, desc: "Promoter: " + promoter, status: yp_status(promoterScore, [60, 80]) },
      { name: "拷贝数", value: copyScore, max: 100, desc: copy + " copies", status: yp_status(copyScore, [40, 80]) },
      { name: "整合位点", value: integScore, max: 100, desc: "Site: " + integ, status: yp_status(integScore, [50, 80]) },
      { name: "mRNA稳定性", value: mrna, max: 100, desc: "mRNA半寿期估计", status: yp_status(mrna, [50, 75]) },
      { name: "终止子效率", value: term, max: 100, desc: "Terminator efficiency", status: yp_status(term, [60, 80]) }
    ];
    var rec = [];
    if (promoterScore < 80) rec.push("考虑使用 AOX1 或 GAP 强启动子以提高转录水平");
    if (copy < 2) rec.push("增加基因拷贝数可线性提升 mRNA 产量");
    if (integ === "random") rec.push("随机整合可能导致位置效应，建议 AOX1 位点定向整合");
    if (mrna < 60) rec.push("优化 5'/3' UTR 提升 mRNA 稳定性");
    return { score: Math.round(weighted), grade: yp_grade(weighted), subScores: sub, tools: ["ExpressYeaself", "Promoter Calculator"], recommendations: rec };
  }
  function yp_L2(auto, caiN, gcN) {
    var cai = (caiN != null && !isNaN(caiN)) ? Math.max(0, Math.min(1, caiN)) : auto.cai;
    var gc = (gcN != null && !isNaN(gcN)) ? Math.max(0, Math.min(1, gcN)) : auto.gc;
    var caiScore = cai * 100;
    var gcPen = Math.abs(gc - 0.5) > 0.15 ? 20 : 0; var gcScore = Math.max(0, 100 - gcPen * 3);
    var hpPen = auto.dG5 < -2 ? 30 : (auto.dG5 < -1 ? 15 : 0); var fpScore = Math.max(0, 100 - hpPen);
    var rareScore = Math.max(0, 100 - auto.rare * 5);
    var weighted = caiScore * 0.35 + gcScore * 0.15 + fpScore * 0.25 + rareScore * 0.25;
    var sub = [
      { name: "CAI 密码子适应指数", value: Math.round(caiScore), max: 100, desc: "CAI = " + cai.toFixed(3), status: yp_status(caiScore, [60, 80]) },
      { name: "GC含量合理性", value: Math.round(gcScore), max: 100, desc: "GC = " + (gc * 100).toFixed(1) + "%", status: yp_status(gcScore, [60, 80]) },
      { name: "5'端发夹结构", value: Math.round(fpScore), max: 100, desc: "ΔG = " + auto.dG5.toFixed(1) + " kcal/mol", status: yp_status(fpScore, [60, 85]) },
      { name: "稀有密码子", value: Math.round(rareScore), max: 100, desc: auto.rare + " 个稀有密码子", status: yp_status(rareScore, [60, 80]) }
    ];
    var rec = [];
    if (cai < 0.6) rec.push("CAI偏低，建议用酵母偏好密码子优化编码序列");
    if (auto.dG5 < -2) rec.push("5'端发夹过强，可能阻碍核糖体扫描，建议优化编码序列前段");
    if (auto.rare > 10) rec.push("稀有密码子过多，可能导致翻译暂停与错误折叠");
    if (Math.abs(gc - 0.5) > 0.15) rec.push("GC含量偏离理想范围(40-60%)，可能影响 mRNA 稳定性");
    return { score: Math.round(weighted), grade: yp_grade(weighted), subScores: sub, tools: ["GenSmart Codon Optimization", "IDT Codon Optimization", "JCat"], recommendations: rec };
  }
  function yp_L3(auto) {
    var plddt = auto.plddt, disorderPen = auto.disorderRegions * 15, disorderScore = Math.max(0, 100 - disorderPen);
    var aggScore = 100 - auto.agg, solScore = auto.sol * 100;
    var weighted = plddt * 0.30 + disorderScore * 0.20 + aggScore * 0.25 + solScore * 0.25;
    var sub = [
      { name: "pLDDT 结构置信度", value: Math.round(plddt), max: 100, desc: "Avg pLDDT = " + plddt.toFixed(1), status: yp_status(plddt, [60, 80]) },
      { name: "无序区评估", value: Math.round(disorderScore), max: 100, desc: auto.disorderRegions + " 个无序区(>30aa)", status: yp_status(disorderScore, [60, 80]) },
      { name: "聚集倾向 Aggregation", value: Math.round(aggScore), max: 100, desc: "Aggregation = " + auto.agg.toFixed(1), status: yp_status(aggScore, [60, 80]) },
      { name: "可溶性 Camsol", value: Math.round(solScore), max: 100, desc: "Solubility = " + auto.sol.toFixed(3), status: yp_status(solScore, [50, 75]) }
    ];
    var rec = [];
    if (plddt < 70) rec.push("平均 pLDDT 较低，蛋白可能存在大量无序区或错误折叠倾向");
    if (auto.disorderRegions > 2) rec.push("存在多个长无序区，建议截短体设计或融合稳定结构域");
    if (auto.agg > 50) rec.push("聚集倾向较高，易形成包涵体，可降培养温度或共表达分子伴侣");
    if (auto.sol < 0.4) rec.push("预测可溶性较差，建议表面残基突变或融合标签");
    return { score: Math.round(weighted), grade: yp_grade(weighted), subScores: sub, tools: ["AlphaFold", "ColabFold", "Camsol"], recommendations: rec };
  }
  function yp_L4(auto) {
    var cys = auto.cys;
    var cysScore = cys <= 4 ? 90 : cys <= 6 ? 70 : cys <= 10 ? 50 : 30;
    var disScore = auto.disBonds <= 3 ? 90 : auto.disBonds <= 6 ? 65 : 40;
    var mwScore = auto.mw < 50 ? 95 : auto.mw < 80 ? 80 : auto.mw < 120 ? 60 : 40;
    var glyScore = auto.nSites <= 3 ? 90 : auto.nSites <= 5 ? 70 : 45;
    var uprScore = 100 - auto.upr;
    var weighted = cysScore * 0.15 + disScore * 0.20 + mwScore * 0.20 + glyScore * 0.20 + uprScore * 0.25;
    var sub = [
      { name: "半胱氨酸数量", value: cysScore, max: 100, desc: cys + " Cys", status: yp_status(cysScore, [50, 75]) },
      { name: "二硫键配对", value: disScore, max: 100, desc: auto.disBonds + " 对预测", status: yp_status(disScore, [50, 75]) },
      { name: "分子量", value: mwScore, max: 100, desc: auto.mw.toFixed(1) + " kDa", status: yp_status(mwScore, [50, 80]) },
      { name: "糖基化负荷", value: glyScore, max: 100, desc: auto.nSites + " N-糖基化位点", status: yp_status(glyScore, [50, 75]) },
      { name: "UPR风险", value: Math.round(uprScore), max: 100, desc: "UPR risk = " + auto.upr + "/100", status: yp_status(uprScore, [50, 75]) }
    ];
    var rec = [];
    if (cys > 6) rec.push("半胱氨酸 > 6，二硫键负荷大，建议共表达 PDI（蛋白二硫键异构酶）");
    if (auto.mw > 80) rec.push("分子量 > 80 kDa，ER加工负担重，可考虑分段表达");
    if (auto.nSites > 5) rec.push("糖基化位点 > 5，可能超糖基化与 ER 负荷增加");
    if (auto.upr > 60) rec.push("UPR激活风险高，建议降培养温度、化学伴侣或减弱启动子");
    return { score: Math.round(weighted), grade: yp_grade(weighted), subScores: sub, tools: ["DiANNA", "UPR Predictor"], recommendations: rec };
  }
  function yp_L5(auto, seqClean) {
    var n = auto.nSites; var oCount = yp_oSites(seqClean);
    var nScore = n === 0 ? 95 : n <= 2 ? 85 : n <= 4 ? 65 : 40;
    var oScore = oCount <= 5 ? 85 : oCount <= 10 ? 65 : 45;
    var hyper = Math.min(100, n * 18); var erLoad = Math.min(100, n * 15 + (auto.mw > 80 ? 20 : 0));
    var probScore = n > 0 ? 70 : 10;
    var weighted = nScore * 0.25 + oScore * 0.15 + (100 - hyper) * 0.25 + (100 - erLoad) * 0.20 + probScore * 0.15;
    var sub = [
      { name: "N-糖基化位点数", value: nScore, max: 100, desc: n + " N-X-S/T", status: yp_status(nScore, [50, 75]) },
      { name: "O-糖基化位点数", value: oScore, max: 100, desc: oCount + " O-糖基", status: yp_status(oScore, [50, 75]) },
      { name: "超糖基化风险", value: Math.round(100 - hyper), max: 100, desc: "Risk = " + hyper + "/100", status: yp_status(100 - hyper, [50, 75]) },
      { name: "ER负荷", value: Math.round(100 - erLoad), max: 100, desc: "ER load = " + erLoad + "/100", status: yp_status(100 - erLoad, [50, 75]) }
    ];
    var rec = [];
    if (n > 4) rec.push("N-糖基化位点过多(>4)，可能超糖基化与分泌下降");
    if (hyper > 60) rec.push("超糖基化风险高，建议定点突变消除非关键位点(N→Q)");
    if (erLoad > 60) rec.push("糖基化致 ER 负荷重，建议糖基化缺陷型宿主(如 Och1 突变)");
    return { score: Math.round(weighted), grade: yp_grade(weighted), subScores: sub, tools: ["NetNGlyc", "NetOGlyc", "GlycoEP"], recommendations: rec };
  }
  function yp_L6(auto, sp) {
    var exists = auto.hasSP; var prob = auto.spProb;
    var spScore = exists ? Math.round(prob * 100) : 20;
    var tmScore = 100; var secProb = auto.secProb * 100;
    var weighted = spScore * 0.35 + tmScore * 0.25 + secProb * 0.40;
    var sub = [
      { name: "信号肽", value: spScore, max: 100, desc: exists ? ("SP P=" + prob.toFixed(2)) : "No signal peptide", status: yp_status(spScore, [50, 75]) },
      { name: "跨膜区", value: tmScore, max: 100, desc: "0 跨膜区(简化估算)", status: yp_status(tmScore, [60, 85]) },
      { name: "分泌定位概率", value: Math.round(secProb), max: 100, desc: "P(secreted) = " + (auto.secProb * 100).toFixed(1) + "%", status: yp_status(secProb, [50, 75]) },
      { name: "胞质定位概率", value: Math.round(auto.cytoProb * 100), max: 100, desc: "P(cyto) = " + (auto.cytoProb * 100).toFixed(1) + "%", status: auto.cytoProb > 0.3 ? "danger" : "good" }
    ];
    var rec = [];
    if (!exists) rec.push("未检测到信号肽，蛋白可能无法进入分泌途径，建议添加 α-factor 或 MFα1 信号肽");
    if (auto.secProb < 0.5) rec.push("分泌概率较低，建议优化信号肽或 Kex2 切割位点");
    return { score: Math.round(weighted), grade: yp_grade(weighted), subScores: sub, tools: ["SignalP 6.0", "Phobius", "WoLF PSORT"], recommendations: rec };
  }
  function yp_L7(auto) {
    var deepSec = auto.deepSec * 100, expr = auto.pTexpr * 100, sol = auto.pTsol * 100, sec = auto.pTsec * 100;
    var locs = auto.loc; var top = Object.keys(locs).sort(function (a, b) { return locs[b] - locs[a]; })[0];
    var locScore = top === "Secreted" ? locs[top] * 100 : (1 - locs[top]) * 100;
    var weighted = deepSec * 0.25 + expr * 0.25 + sol * 0.20 + sec * 0.20 + locScore * 0.10;
    var sub = [
      { name: "DeepSec 分泌倾向", value: Math.round(deepSec), max: 100, desc: "DeepSec = " + auto.deepSec.toFixed(3), status: yp_status(deepSec, [50, 75]) },
      { name: "ProtTrans 表达量", value: Math.round(expr), max: 100, desc: "Expression = " + auto.pTexpr.toFixed(3), status: yp_status(expr, [50, 75]) },
      { name: "ProtTrans 可溶性", value: Math.round(sol), max: 100, desc: "Solubility = " + auto.pTsol.toFixed(3), status: yp_status(sol, [50, 75]) },
      { name: "ProtTrans 分泌能力", value: Math.round(sec), max: 100, desc: "Secretion = " + auto.pTsec.toFixed(3), status: yp_status(sec, [50, 75]) },
      { name: "DeepLoc 定位", value: Math.round(locScore), max: 100, desc: "Top: " + top + " (" + (locs[top] * 100).toFixed(1) + "%)", status: yp_status(locScore, [50, 75]) }
    ];
    var rec = [];
    if (auto.deepSec < 0.5) rec.push("DeepSec 分泌倾向低，建议综合优化信号肽与表面电荷");
    if (auto.pTexpr < 0.4) rec.push("pLM 表达量偏低，尝试多拷贝整合或强启动子");
    if (auto.pTsol < 0.4) rec.push("pLM 可溶性差，建议融合标签(SUMO/MBP)或共表达折叠辅助因子");
    return { score: Math.round(weighted), grade: yp_grade(weighted), subScores: sub, tools: ["DeepLoc 2.0", "DeepSec", "ProtTrans", "ESM2"], recommendations: rec };
  }
  function yeastAssess(seq, promoter, sp, auto, caiN, gcN) {
    var clean = yp_clean(seq);
    return [yp_L1(promoter, auto), yp_L2(auto, caiN, gcN), yp_L3(auto), yp_L4(auto), yp_L5(auto, clean), yp_L6(auto, sp), yp_L7(auto)];
  }

  function yp_bar(v, max, color) {
    var pct = Math.max(0, Math.min(100, v / max * 100));
    var col = color || (pct >= 75 ? "#1f9e89" : pct >= 50 ? "#ffce4d" : "#ff5a5a");
    return '<div style="display:flex;align-items:center;gap:8px"><div class="layer-bar" style="flex:1;min-width:120px"><div class="layer-fill" style="width:' + pct + '%;background:' + col + '"></div></div><span style="font-size:12px;color:#cdd9e5;min-width:38px;text-align:right">' + Math.round(v) + '</span></div>';
  }
  function yp_subTable(sub) {
    return '<table class="tbl"><tr><th>指标</th><th style="width:55%">评分</th><th>状态</th></tr>' +
      sub.map(function (s) {
        var st = { good: "🟢", warning: "🟡", danger: "🔴" }[s.status] || "";
        return "<tr><td>" + esc(s.name) + '<div class="ksub">' + esc(s.desc) + '</div></td><td>' + yp_bar(s.value, s.max) + '</td><td style="color:' + (s.status === "good" ? "#3ddc97" : s.status === "warning" ? "#ffce4d" : "#ff5a5a") + '">' + st + "</td></tr>";
      }).join("") + "</table>";
  }
  function yp_kv(label, value) { return '<div class="ec-kv"><span>' + esc(label) + "</span><b>" + esc(String(value)) + "</b></div>"; }
  function levelNarrative(L, auto) {
    if (L === "L1") return "转录层评估启动子强度(PAOX1≈95 / GAP≈85 / TEF1≈80)、基因拷贝数、整合位点(AOX1_locus/HIS4/random)与 mRNA 稳定性：启动子×0.30 + 拷贝×0.20 + 整合×0.20 + mRNA×0.15 + 终止子×0.15。";
    if (L === "L2") return "翻译层评估 CAI、GC含量(理想40-60%)、5'端发夹(ΔG<-4强)与稀有密码子：CAI×0.35 + GC×0.15 + 5'发夹×0.25 + 稀有×0.25。参考工具 GenSmart/IDT/JCat 仅作设计参考，非实际调用。";
    if (L === "L3") return "折叠层评估 pLDDT(结构置信度)、无序区、聚集倾向(疏水簇)与可溶性(CamSol)：pLDDT×0.30 + 有序×0.20 + 抗聚集×0.25 + 可溶×0.25。参考 AlphaFold/ColabFold。";
    if (L === "L4") return "ER加工层评估半胱氨酸数、二硫键、分子量(<80kDa佳)、糖基化负荷与 UPR风险：半胱×0.15 + 二硫×0.20 + 分子量×0.20 + 糖基化×0.20 + 抗UPR×0.25。参考 DiANNA。";
    if (L === "L5") return "糖基化层预测 N-糖基化 N-X-S/T 与 O-糖基化位点，评估超糖基化风险与 ER 负荷：N×0.25 + O×0.15 + 抗超糖×0.25 + 抗ER负荷×0.20 + 概率×0.15。参考 NetNGlyc/NetOGlyc（非实际调用）。";
    if (L === "L6") return "分泌层评估信号肽(α-MF/α-factor等)、跨膜区与亚细胞定位(WoLF PSORT)：信号肽×0.35 + 跨膜×0.25 + 分泌定位×0.40。参考 SignalP 6.0/Phobius。";
    if (L === "L7") return "全流程预测整合 DeepSec 分泌倾向、ProtTrans 表达量/可溶性/分泌、DeepLoc 定位：DeepSec×0.25 + 表达×0.25 + 可溶×0.20 + 分泌×0.20 + 定位×0.10。参考 DeepLoc 2.0/DeepSec/ProtTrans/ESM2（非实际调用）。";
    return "";
  }
  function renderLevelDetail(m, res, auto, r7) {
    var html = "";
    html += '<div style="display:flex;gap:14px;flex-wrap:wrap">';
    html += card('<div class="sub">' + m.L + ' 评分</div><div style="font-size:34px;font-weight:800;color:' + (res.score >= 75 ? "#1f9e89" : res.score >= 60 ? "#3ddc97" : res.score >= 40 ? "#ffce4d" : "#ff5a5a") + '">' + res.score + '<span style="font-size:14px;color:#9fb3c8"> /100 · ' + res.grade + '级</span></div>' + yp_bar(res.score, 100), 'style="flex:1;min-width:200px"');
    html += card('<div class="sub">详细子项评分</div>' + yp_subTable(res.subScores), 'style="flex:2;min-width:340px"');
    html += "</div>";
    if (res.recommendations && res.recommendations.length) html += card('<div class="sub">优化建议</div><ul class="notes">' + res.recommendations.map(function (r) { return "<li>" + esc(r) + "</li>"; }).join("") + "</ul>");
    html += card('<div class="sub">算法说明（' + m.en + '）</div><div class="hint">' + esc(levelNarrative(m.L, auto)) + "</div>");
    if (m.L === "L7") {
      html += card('<div class="sub">全流程综合（7层加权）</div>' +
        yp_kv("L1-L7 加权总评", (r7.overall * 100).toFixed(1) + "/100") +
        yp_kv("预测表达等级", { high: "高", mid: "中等", low: "低", none: "不表达" }[r7.expression_level]) +
        yp_kv("预测产量", r7.yield_mg_l + " mg/L") +
        '<div class="hint">' + esc(r7.evidence) + "</div>");
    } else {
      html += card('<div class="sub">对总评的贡献</div><div class="hint">本层权重 ' + fmt(m.w, 2) + "，当前层评分 " + res.score + "/100，对 7 层加权总评贡献约 " + (m.w * res.score / 100 * 100).toFixed(1) + " 分。</div>");
    }
    return html;
  }
  function yeastLevelParamCard() {
    var e = state.expr || {};
    function opt(a, cur) { return a.map(function (o) { return '<option value="' + o[1] + '"' + (cur === o[1] ? " selected" : "") + ">" + o[0] + "</option>"; }).join(""); }
    return card(
      '<div class="sub">参数输入（与酵母总览共享）</div>' +
      '<div class="formgrid">' +
      '<label>宿主菌株</label><select id="e_host" onchange="FIP.yeastRecompute()">' + opt([["P. pastoris GS115", "Pichia_pastoris_GS115"], ["P. pastoris X33", "Pichia_pastoris_X33"]], e.host || "Pichia_pastoris_GS115") + '</select>' +
      '<label>信号肽</label><select id="e_sp" onchange="FIP.yeastRecompute()">' + opt(YEAST_SP_OPTIONS, e.sp || "alpha-factor") + '</select>' +
      '<label>启动子</label><select id="e_promoter" onchange="FIP.yeastRecompute()">' + opt(YEAST_PROMOTERS, e.promoter || "PAOX1") + '</select>' +
      '<label>CAI</label><input id="e_cai" value="' + esc(e.cai || "") + '" style="width:90px" oninput="FIP.yeastRecompute()">' +
      '<label>GC%</label><input id="e_gc" value="' + esc(e.gc || "") + '" style="width:90px" oninput="FIP.yeastRecompute()">' +
      "</div>" +
      '<div class="ksub" style="margin:6px 0 2px">氨基酸序列</div>' +
      '<textarea id="e_seq" rows="3" style="width:100%;font-family:monospace" placeholder="MALWMRLLPLLALLALWGPDPAAA..." oninput="FIP.yeastRecompute()">' + esc(e.seq || "") + "</textarea>" +
      '<div class="count"><span id="e_count">' + (e.seq ? e.seq.length : 0) + " aa</span> · 示例：" +
      '<span class="exbtn" onclick="FIP.example(\'hsa\')">HSA</span> <span class="exbtn" onclick="FIP.example(\'gfp\')">GFP</span> <span class="exbtn" onclick="FIP.example(\'insulin\')">Insulin</span></div>'
    );
  }
  function yeastLevelRecompute(idx) {
    var seq = $("e_seq") ? $("e_seq").value : "";
    var host = $("e_host") ? $("e_host").value : "Pichia_pastoris_GS115";
    var promoter = $("e_promoter") ? $("e_promoter").value : "PAOX1";
    var sp = $("e_sp") ? $("e_sp").value : "alpha-factor";
    var cai = $("e_cai") ? $("e_cai").value : "";
    var gc = $("e_gc") ? $("e_gc").value : "";
    state.expr = state.expr || {};
    state.expr.host = host; state.expr.promoter = promoter; state.expr.sp = sp; state.expr.seq = seq; state.expr.cai = cai; state.expr.gc = gc;
    if ($("e_count")) $("e_count").textContent = (seq ? seq.length : 0) + " aa";
    var box = $("yl_detail"); if (!box) return;
    var caiN = cai === "" ? null : parseFloat(cai); if (caiN != null && caiN > 1.5) caiN = caiN / 100;
    var gcN = gc === "" ? null : parseFloat(gc); if (gcN != null && gcN > 1.5) gcN = gcN / 100;
    var auto = yp_auto(seq, sp);
    var det = yeastAssess(seq, promoter, sp, auto, caiN, gcN);
    var r7 = sevenLayer(host, promoter, sp, seq, cai, gc, null);
    box.innerHTML = renderLevelDetail(YEAST_LEVELS[idx], det[idx], auto, r7);
  }
  function renderYeastLevel(idx, body) {
    var m = YEAST_LEVELS[idx];
    var html = "";
    html += '<div style="margin-bottom:8px"><span class="exbtn" onclick="FIP.nav(\'yeast\')">← 返回酵母总览</span></div>';
    html += yeastLevelParamCard();
    html += card(
      '<div style="display:flex;justify-content:space-between;align-items:center;flex-wrap:wrap;gap:8px">' +
      '<div class="sub" style="margin:0">' + m.L + " · " + m.name + ' <span class="ksub" style="color:#9fb3c8">' + m.en + "</span></div>" +
      '<span class="ec-badge ec-info">基于规则预估 · 非深度学习推理</span>' +
      '<span class="ksub">权重 ' + fmt(m.w, 2) + "</span></div>" +
      '<div class="ysub" style="margin:6px 0">核心评估：' + esc(m.eval) + "</div>" +
      '<div class="ksub" style="margin-bottom:6px">参考工具（标注"非实际调用"）：' + m.tools.map(function (t) { return '<span class="ec-badge ec-good" style="margin:1px">' + esc(t) + "</span>"; }).join("") + "</div>" +
      '<div id="yl_detail"></div>'
    );
    body.innerHTML = html;
    yeastLevelRecompute(idx);
  }

  // ---- E. coli Expression Advisor (参考 ecoli-expression-advisor skill) ----
  var ECOLI_SAMPLE = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAVGVCIKVEKGDAPTLDVITLE";
  var ECOLI_STRAIN_ORDER = ["K12", "BL21", "SHUFFLE", "ROSETTA"];
  var ECOLI_CLASS_COLORS = ["#ff6b6b", "#ffd166", "#80ed99", "#4cc9f0"];

  function ecBadge(text, type) {
    return '<span class="ec-badge ec-' + type + '">' + esc(text) + '</span>';
  }
  function ecBar(v, color) {
    var pct = Math.round(Math.max(0, Math.min(1, v)) * 100);
    return '<div style="display:flex;align-items:center;gap:8px"><div class="ec-bar"><div class="ec-fill" style="width:' + pct + '%;background:' + (color || "#4cc9f0") + '"></div></div><span class="ec-barval">' + fmt(v, 2) + '</span></div>';
  }
  function ecKV(label, value) {
    return '<div class="ec-kv"><span>' + esc(label) + '</span><b>' + esc(String(value)) + '</b></div>';
  }

  function renderExpressionEcoli(body) {
    var e = state.ecoli;
    var html = "";
    // Hero
    html += '<div class="yhero"><div class="ytitle" style="color:#4cc9f0">E. coli 异源蛋白表达顾问</div><div class="ysub">大肠杆菌 E. coli · 纯前端离线版</div>' +
      '<div class="ytrans">基于 Jiang et al. 2024 (<em>Biotechnology Advances</em>) 根因诊断→设计框架。输入蛋白/DNA 序列，在 K12 / BL21 / SHuffle / Rosetta 间切换，按菌株各自密码子表重算表达潜力。纯前端、序列不出浏览器。本模块额外提供周质空间定位与分泌表达预测评分。</div></div>';
    // 根因图例
    html += '<div style="display:flex;gap:10px;flex-wrap:wrap;margin-bottom:14px">' +
      '<div class="card" style="flex:1;min-width:200px"><b style="color:#4cc9f0">根因1 蛋白毒性</b><br><span class="ksub">破坏宿主生理 → 严控基础表达 / 分泌表达</span></div>' +
      '<div class="card" style="flex:1;min-width:200px"><b style="color:#4cc9f0">根因2 密码子偏置</b><br><span class="ksub">稀有密码子→tRNA 耗竭 → 密码子优化 / Rosetta</span></div>' +
      '<div class="card" style="flex:1;min-width:200px"><b style="color:#4cc9f0">根因3 mRNA 结构</b><br><span class="ksub">5\' 发夹阻碍起始 → 去结构化设计</span></div></div>';
    // Section 1 输入
    html += '<div style="display:flex;gap:14px;flex-wrap:wrap">';
    html += card(
      '<div class="sub">1. 输入序列</div>' +
      '<textarea id="e2_seq" rows="4" style="width:100%;font-family:monospace" placeholder="粘贴氨基酸序列（如 MKT...）或 DNA 编码序列（ATG...）。也可点击「载入示例」。">' + esc(e.seq || "") + '</textarea>' +
      '<label style="display:inline-flex;gap:6px;align-items:center;margin-top:8px;font-size:13px;color:#9fb3c8"><input type="checkbox" id="e2_tox" ' + (e.tox ? "checked" : "") + '> 疑似毒性/膜蛋白</label> ' +
      '<label style="display:inline-flex;gap:6px;align-items:center;margin-top:8px;font-size:13px;color:#9fb3c8"><input type="checkbox" id="e2_secret" ' + (e.secret ? "checked" : "") + '> 走分泌表达</label>' +
      '<div style="display:flex;gap:10px;align-items:center;margin-top:8px"><span class="ksub">目标菌株：</span><select id="e2_strain">' + ECOLI_STRAIN_ORDER.map(function (k) { return '<option value="' + k + '"' + (e.strain === k ? " selected" : "") + '>' + ECOLI.STRAINS[k].label + '</option>'; }).join("") + '</select></div>' +
      '<div style="display:flex;gap:8px;margin-top:10px"><button class="btn" style="background:#4cc9f0;color:#06212b" onclick="FIP.ecoliAnalyze()">▶ 分析并给出方案</button>' +
      '<button class="btn2" onclick="FIP.ecoliExample()">载入示例</button><button class="btn2" onclick="FIP.ecoliClear()">清空</button></div>' +
      '<div class="hint" style="margin-top:8px">启发式按所选菌株密码子表重算稀有密码子/CAI/疏水性；5\' 发夹为序列固有。Nussinov DP 本地计算，建议单次 ≤600 aa。</div>',
      'style="flex:1.1;min-width:360px"');
    html += '<div id="e2_diag" style="flex:1.2;min-width:360px"></div>';
    html += "</div>";
    // Section 3 跨菌株对比
    html += '<div id="e2_compare"></div>';
    // Section 4 5' 结构
    html += '<div id="e2_fold"></div>';
    // Section 5 表达预测
    html += '<div id="e2_ml"></div>';
    // Section 6 周质空间 & 分泌表达
    html += card(
      '<div class="sub">周质空间 & 分泌表达预测评分</div><div class="ysub" style="margin-bottom:8px">在酵母/胞内表达之外，额外评估分泌定位（参考 E. coli 信号肽池）。</div>' +
      '<div class="formgrid"><label>信号肽</label><select id="e2_sp" onchange="FIP.ecoliSetSp(this.value)">' + Object.keys(ECOLI.SIGNAL_PEPTIDES_EC).map(function (k) { return '<option value="' + k + '"' + (e.sp === k ? " selected" : "") + '>' + k + '（' + ECOLI.SIGNAL_PEPTIDES_EC[k].route + '）</option>'; }).join("") + '</select></div>' +
      '<div id="e2_peri" style="margin-top:8px"></div>',
      'style="margin-top:14px"');
    // Section 7 反馈飞轮
    html += card(
      '<div class="sub">实验数据回流飞轮（localStorage，离线持久化）</div>' +
      '<div class="ksub" style="margin-bottom:6px">把真实实验结果回填，系统据此校准预测（演示数据飞轮）。数据仅存于本浏览器。</div>' +
      '<div style="display:flex;gap:8px;flex-wrap:wrap;align-items:center">' +
      '<span class="ksub">表达等级</span><select id="e2_fbo"><option value="0">0 不表达</option><option value="1">1 低</option><option value="2">2 中</option><option value="3">3 高</option></select>' +
      '<span class="ksub">溶解度%</span><input id="e2_fbs" value="' + e.fbSol + '" style="width:64px">' +
      '<span class="ksub">包涵体%</span><input id="e2_fbi" value="' + e.fbInc + '" style="width:64px">' +
      '<span class="ksub">产量mg/L</span><input id="e2_fby" value="' + e.fbYield + '" style="width:72px">' +
      '<button class="btn" style="background:#4cc9f0;color:#06212b" onclick="FIP.ecoliSaveFb()">+ 保存实验结果</button></div>' +
      '<div style="display:flex;gap:8px;margin-top:8px;align-items:center">' +
      '<button class="btn2" style="background:#ffd166;color:#3a2c00;border-color:#ffd166" onclick="FIP.ecoliTrainFb()">用反馈数据校准模型</button>' +
      '<button class="btn2" onclick="FIP.ecoliClearFb()">清空反馈</button>' +
      '<span id="e2_fbstat" class="ksub"></span></div>' +
      '<div class="hint" style="margin-top:8px">校准逻辑：用已存反馈的「序列特征→实际表达等级」做最小二乘拟合，得到偏置项作用于启发式表达潜力；≥3 条反馈即可校准。</div>',
      'style="margin-top:14px"');
    // 实验设计信息量反馈 (DOE)：设计期构型覆盖 + 批次/工艺设定点覆盖（show_construct=true，对应 Streamlit Expression 页）
    html += card(doeSectionHtml({ host: "E.coli_BL21(DE3)", showConstruct: true, scale: null, withHostSelect: true }));
    body.innerHTML = html;
    ecoliRecompute();
    renderEcoliFbStat();
  }

  function ecoliRecompute() {
    var e = state.ecoli;
    var seq = $("e2_seq") ? $("e2_seq").value : "";
    var strain = $("e2_strain") ? $("e2_strain").value : e.strain;
    var tox = $("e2_tox") ? $("e2_tox").checked : e.tox;
    var secret = $("e2_secret") ? $("e2_secret").checked : e.secret;
    var sp = $("e2_sp") ? $("e2_sp").value : e.sp;
    e.seq = seq; e.strain = strain; e.tox = tox; e.secret = secret; e.sp = sp;
    var r = ECOLI.analyze(seq, strain, tox, secret);
    var diag = $("e2_diag");
    if (!diag) return;
    if (!r) {
      diag.innerHTML = card('<div class="sub">2. 诊断与推荐（当前菌株）</div><div class="ksub" style="margin-top:6px">等待输入…含稀有密码子簇的序列会触发「密码子优化 + Rosetta(pRARE)」；勾选毒性会触发 C41/C43 与自诱导建议。切换菌株可即时重算。</div>');
      $("e2_compare").innerHTML = ""; $("e2_fold").innerHTML = ""; $("e2_ml").innerHTML = ""; $("e2_peri").innerHTML = "";
      return;
    }
    var rec = ECOLI.recommend(r, strain);
    // 诊断面板
    var dh = '<div class="sub">2. 诊断与推荐（' + ECOLI.STRAINS[strain].label + '）</div>' + ecKV("序列类型", r.fromProtein ? "蛋白(已回译为CDS)" : "DNA(编码链)") +
      ecKV("蛋白长度", r.proteinLen + " aa") + ecKV("GC 含量(全长)", (r.gcFull * 100).toFixed(1) + "%") +
      ecKV("CAI(相对本菌株)", r.cai.toFixed(3)) + ecKV("5' 区 GC", (r.gc5 * 100).toFixed(1) + "%") +
      ecKV("5' 配对比例", (r.paired5Frac * 100).toFixed(1) + "%") + ecKV("疏水性占比", (r.hydroFrac * 100).toFixed(1) + "%") +
      ecKV("半胱氨酸(二硫键)", String(r.cysteine));
    dh += '<div style="margin:8px 0">' + rec.badges.map(function (b) { return ecBadge(b.text, b.type); }).join("") + '</div>';
    if (r.rareList.length) {
      dh += '<div class="ksub" style="margin:6px 0 2px">本菌株判定的稀有密码子（W<' + ECOLI.RARE_THR + '）</div>';
      var rset = {}; r.rareList.forEach(function (o) { rset[o.c] = (rset[o.c] || 0) + 1; });
      dh += '<table class="tbl"><tr><th>密码子</th><th>AA</th><th>次数</th><th>W</th><th>pRARE</th></tr>' +
        Object.keys(rset).sort().map(function (c) {
          var sup = ECOLI.STRAINS[strain].supp.indexOf(c) >= 0;
          return "<tr><td style='color:#ff8b8b;font-weight:700'>" + c + "</td><td>" + ECOLI.AA_MAP[c] + "</td><td>" + rset[c] + "</td><td>" + ECOLI.STRAINS[strain].table.w[c].toFixed(3) + "</td><td>" + (sup ? "✅" : "否") + "</td></tr>";
        }).join("") + "</table>";
    } else {
      dh += '<div class="ksub" style="margin-top:6px">该序列在 ' + ECOLI.STRAINS[strain].label + ' 下无判定为稀有的密码子。</div>';
    }
    dh += '<div style="margin-top:8px"><div class="ec-kv"><span>宿主菌株</span><b>' + esc(rec.host) + '</b></div>' +
      '<div class="ec-kv"><span>诱导策略</span><b>' + esc(rec.induction) + '</b></div>' +
      '<div class="ec-kv"><span>融合标签</span><b>' + esc(rec.tag) + '</b></div>' +
      (rec.signalPeptide ? '<div class="ec-kv"><span>分泌信号肽</span><b>' + esc(rec.signalPeptide) + '</b></div>' : '') +
      '<div class="ec-kv"><span>序列重设计</span><b>' + esc(rec.optimization.join("；")) + '</b></div></div>';
    dh += '<div class="sub" style="margin-top:8px">实验 SOP 清单</div><ul class="notes">' + rec.sopItems.map(function (it) { return "<li>" + esc(it) + "</li>"; }).join("") + '</ul>';
    dh += '<div class="hint" style="margin-top:6px">' + esc(ECOLI.STRAINS[strain].note) + '</div>';
    diag.innerHTML = card(dh);

    // 跨菌株对比
    var rows = ECOLI.crossStrainCompare(r.dna, r.gc5, r.hydroFrac, r.paired5Frac);
    var ch = '<div class="sub">3. 跨菌株对比（同一序列，按各菌株密码子表重算）</div><div style="overflow:auto"><table class="tbl"><tr><th>指标</th>' +
      ECOLI_STRAIN_ORDER.map(function (k) { return "<th>" + ECOLI.STRAINS[k].label.replace("E. coli ", "") + "</th>"; }).join("") + "</tr>" +
      rows.map(function (row) {
        return "<tr><td>" + esc(row.label) + "</td>" + row.values.map(function (v, j) {
          return "<td" + (row.bestIdx === j ? ' style="color:#80ed99;font-weight:700"' : "") + ">" + esc(v) + "</td>";
        }).join("") + "</tr>";
      }).join("") + "</table></div>" +
      '<div class="hint" style="margin-top:6px">高亮列为该指标最优菌株。序列固有指标（5\'发夹/5\'GC/疏水性）所有菌株相同；Rosetta 因 pRARE 回补不再判 AGA/AGG/ATA/CTA/CCC/GGA/ACA/TCA 为稀有，但 <b style="color:#80ed99">TTA、CGA、CGG 仍稀有</b>。</div>';
    $("e2_compare").innerHTML = card(ch);

    // 5' 结构
    var fh = '<div class="sub">4. 5\' 起始区 mRNA 二级结构（Nussinov 动态规划折叠）</div>' +
      ecKV("5' 区长度", r.headLen + " nt") + ecKV("近似 MFE", r.mfe.toFixed(1) + " kcal/mol") + ecKV("配对碱基比例", (r.paired5Frac * 100).toFixed(1) + "%") +
      '<div class="ec-rna">' + esc(r.headRna) + '</div>' +
      '<div class="ec-rna">' + r.hairpinDot.split("").map(function (ch2) { return '<span style="color:' + (ch2 === "(" || ch2 === ")" ? "#ff8b8b" : "#6b7686") + '">' + ch2 + "</span>"; }).join("") + '</div>' +
      '<div class="hint" style="margin-top:4px">红/灰为点括号：<b style="color:#ff8b8b">( )</b> 配对茎区，<b style="color:#6b7686">·</b> 未配对环区。MFE 为 Nussinov 最大配对近似（每对 -2.0 kcal/mol），仅用于 5\' 起始区结构风险判断。</div>';
    $("e2_fold").innerHTML = card(fh);

    // 表达预测 ML
    var feats = { cai: r.cai, paired5: r.paired5Frac, rare: r.rareClusters, hyd: r.hydroFrac, tox: r.tox };
    var ml = ECOLI.predictML(feats, r.dna, e.model, ECOLI.STRAINS[strain].table.w, e.bias);
    e._r = r; e._ml = ml;
    var labels = ["不表达", "低表达", "中表达", "高表达"];
    var mh = '<div class="sub">5. 表达预测（可切换备选模型）</div><div style="display:flex;gap:8px;align-items:center;margin-bottom:6px"><span class="ksub">选择模型：</span><select id="e2_model" onchange="FIP.ecoliSetModel(this.value)">' +
      '<option value="builtin"' + (e.model === "builtin" ? " selected" : "") + '>内置启发式（规则+sigmoid）</option>' +
      '<option value="mpepe"' + (e.model === "mpepe" ? " selected" : "") + '>MPEPE 风格代理（高表达概率）</option>' +
      '<option value="deeptesr"' + (e.model === "deeptesr" ? " selected" : "") + '>DeepTESR 风格代理（TESR 短斜坡）</option></select>' +
      '<span class="ksub">当前：' + esc(ml.modelName) + '</span></div>';
    mh += ecKV("表达潜力（0-1）", ml.potential.toFixed(3)) + ecKV("原生输出（" + ml.nativeLabel + "）", ml.nativeValue + (ml.nativeExtra ? " (" + ml.nativeExtra + ")" : ""));
    mh += '<div style="margin:8px 0"><div class="ksub">分类概率</div>';
    labels.forEach(function (lab, i) {
      mh += '<div style="display:flex;justify-content:space-between;font-size:13px;margin:2px 0"><span>' + lab + '</span><b>' + (ml.classProbs[i] * 100).toFixed(1) + '%</b></div>' + ecBar(ml.classProbs[i], ECOLI_CLASS_COLORS[i]).replace('<span', '<span style="display:none"') ;
    });
    mh += '</div>';
    mh += '<div style="display:inline-block;border-radius:14px;padding:3px 12px;font-size:12px;font-weight:700;margin:4px 0;background:' + (ml.classIdx >= 2 ? "rgba(128,237,153,.16)" : ml.classIdx === 1 ? "rgba(255,209,102,.16)" : "rgba(255,107,107,.18)") + ';color:' + (ml.classIdx >= 2 ? "#a7f3c0" : ml.classIdx === 1 ? "#ffe0a3" : "#ff9b9b") + '">预测等级：' + ml.className + '</div>';
    mh += ecKV("预测溶解度", ml.solubility.toFixed(1) + "%") + ecKV("预测产量", ml.yield.toFixed(0) + " mg/L") + ecKV("预测包涵体比例", ml.inclusion.toFixed(1) + "%");
    mh += '<div class="hint" style="margin-top:6px">' + esc(ml.modelNote) + '</div>';
    $("e2_ml").innerHTML = card(mh);

    // 周质空间 & 分泌表达
    var hasDisulfide = r.cysteine > 0;
    var solProxy = ml.solubility / 100;
    var peri = ECOLI.scorePeriplasm(sp, hasDisulfide, strain, r.hydroFrac, solProxy);
    var sec = ECOLI.scoreSecretion(peri.score, sp, strain, r.proteinLen, r.hydroFrac);
    e._peri = peri; e._sec = sec;
    var ph = ecKV("周质空间定位评分", peri.score.toFixed(2)) + ecBar(peri.score, "#4cc9f0") +
      '<div class="ksub" style="margin:4px 0">途径：<b style="color:#4cc9f0">' + peri.route + '</b> · ' + (peri.disulfideOK === null ? "" : (peri.disulfideOK ? "二硫键可正确氧化 ✓" : "二硫键易聚集 ✗")) + '</div>' +
      '<div class="hint">' + esc(peri.note) + '</div>';
    ph += '<div style="margin-top:10px"></div>' + ecKV("分泌表达评分", sec.score.toFixed(2)) + ecBar(sec.score, "#80ed99") +
      ecKV("预测胞外产量", (ml.yield * sec.score).toFixed(0) + " mg/L（周质 " + (ml.yield).toFixed(0) + " × 分泌 " + sec.score.toFixed(2) + "）") +
      '<div class="hint">' + esc(sec.note) + '</div>';
    if (secret) ph += '<div class="ec-badge ec-info" style="margin-top:6px">已勾选「走分泌表达」：建议融合 ' + esc(rec.signalPeptide || "信号肽") + '。</div>';
    $("e2_peri").innerHTML = ph;
  }

  function renderEcoliFbStat() {
    var arr = ECOLI.loadFeedback();
    var counts = [0, 1, 2, 3].map(function (k) { return arr.filter(function (d) { return d.outcome === k; }).length; });
    var st = $("e2_fbstat"); if (st) st.textContent = "已存反馈：" + arr.length + " 条（不/低/中/高 = " + counts.join("/") + "）" + (state.ecoli.bias ? " · 已校准 BIAS=" + state.ecoli.bias.toFixed(2) : "");
  }

  function ecoliSaveFb() {
    var e = state.ecoli;
    var fbSeq = $("e2_seq") ? $("e2_seq").value : "";
    if (!fbSeq && !e._r) { alert("请先分析序列，再保存对应反馈。"); return; }
    var r = e._r || ECOLI.analyze(fbSeq, e.strain, e.tox, e.secret);
    var arr = ECOLI.loadFeedback();
    var entry = {
      ts: Date.now(), outcome: Number($("e2_fbo") ? $("e2_fbo").value : 0),
      solubility: Number($("e2_fbs") ? $("e2_fbs").value : 70),
      inclusion: Number($("e2_fbi") ? $("e2_fbi").value : 30),
      yield: Number($("e2_fby") ? $("e2_fby").value : 120),
      features: { cai: r ? r.cai : 0, paired5: r ? r.paired5Frac : 0, rare: r ? r.rareClusters : 0, hyd: r ? r.hydroFrac : 0, tox: e.tox, length: r ? r.proteinLen : 0, strain: e.strain }
    };
    arr.push(entry); ECOLI.saveFeedback(arr);
    renderEcoliFbStat();
    alert("已保存 " + arr.length + " 条反馈（本浏览器本地，菌株=" + e.strain + "）。点击「校准模型」应用。");
  }
  function ecoliTrainFb() {
    var arr = ECOLI.loadFeedback();
    if (arr.length < 3) { renderEcoliFbStat(); alert("反馈样本不足（需 ≥3 条）暂不校准；当前 " + arr.length + " 条。"); return; }
    var bias = ECOLI.calibrateFromFeedback(arr);
    state.ecoli.bias = bias;
    ecoliRecompute();
    renderEcoliFbStat();
  }
  function ecoliClearFb() { if (confirm("确认清空所有反馈数据？")) { ECOLI.saveFeedback([]); state.ecoli.bias = 0; ecoliRecompute(); renderEcoliFbStat(); } }

  // ---- Scale-up Twin ----
  function renderScaleup(body) {
    var fers = state.scale.fermenters;
    var tags = Object.keys(fers);
    var crit = [["pv", "等体积功率 P/V"], ["kla", "等 kLa（同 P/V）"], ["tip", "等桨尖速度"], ["tmix", "等混合时间"]];
    var s = state.scale;
    if (!fers[s.src]) s.src = tags[0];
    if (!fers[s.tgt] || s.tgt === s.src) s.tgt = tags[tags.length > 1 ? 1 : 0];
    var html = card('<div id="s_fermenters"></div>', 'style="flex:1;min-width:100%"');
    html += '<div style="display:flex;gap:14px;flex-wrap:wrap">';
    html += card('<div class="sub">⚖️ 几何相似放大换算（基于发酵罐真实几何）</div>' +
      '<div style="display:flex;gap:10px;flex-wrap:wrap">' +
      '源尺度<select onchange="FIP.set(\'scale.src\',this.value)">' + tags.map(function (x) { return '<option' + (s.src === x ? " selected" : "") + '>' + x + '</option>'; }).join("") + '</select>' +
      '目标尺度<select onchange="FIP.set(\'scale.tgt\',this.value)">' + tags.map(function (x) { return '<option' + (s.tgt === x ? " selected" : "") + '>' + x + '</option>'; }).join("") + '</select>' +
      '准则<select onchange="FIP.set(\'scale.criterion\',this.value)">' + crit.map(function (x) { return '<option value="' + x[0] + '"' + (s.criterion === x[0] ? " selected" : "") + '>' + x[1] + '</option>'; }).join("") + '</select>' +
      '<button class="btn" onclick="FIP.calcScale()">▶ 换算 + 风险评估</button></div>' +
      '<div id="s_out"></div>', 'style="flex:1.2;min-width:360px"');
    html += card('<div class="sub">📚 放大知识库</div><div id="s_kb" style="max-height:460px;overflow:auto"></div>', 'style="flex:1;min-width:320px"');
    html += '</div>';
    html += card('<div id="s_fchart"></div>', 'style="flex:1.4;min-width:420px"');
    html += card('<div id="s_sweep"></div>', 'style="flex:1;min-width:100%"');
    html += card('<div id="s_risk"></div>', 'style="flex:1;min-width:100%"');
    body.innerHTML = html;
    renderFermenters($("s_fermenters"));
    renderFermenterChart($("s_fchart"));
    renderSweep($("s_sweep"));
    $("s_kb").innerHTML = DATA.knowledge.filter(function (k) { return k.category === "scaleup" || !k.category; }).map(function (k) {
      return '<div class="kb"><b>' + esc(k.title) + '</b><div class="ksub">' + esc(k.content) + (k.source ? '<br>来源：' + esc(k.source) : '') + '</div></div>';
    }).join("");
    renderScaleRisk($("s_risk"));
    calcScale();
  }

  function calcScale() {
    var s = state.scale, fers = s.fermenters;
    var HP = scaleHostProfile(s);
    var host = $("s_out"); if (!host) return;
    if (!fers[s.src] || !fers[s.tgt]) { host.innerHTML = '<div class="hint">请先在上方发酵罐结构数据中添加源/目标尺度。</div>'; return; }
    var fS = fers[s.src], fT = fers[s.tgt];
    var dS = fermenterDerived(fS, s.do_set, s.mu);
    // 几何相似放大：目标转速由准则决定（Di_s/Di_t = 线性放大因子倒数）
    var ratio = (fT.Di > 0 && fS.Di > 0) ? fS.Di / fT.Di : 1;   // Di_s / Di_t
    var Nt;
    if (s.criterion === "tip") Nt = fS.N * ratio;               // 等桨尖速度：N_t = N_s·(Di_s/Di_t)
    else if (s.criterion === "tmix") Nt = fS.N;                 // 等混合时间：N_t = N_s（几何相似下 t_mix∝1/N）
    else Nt = fS.N * Math.pow(ratio, 2 / 3);                    // 等 P/V / 等 kLa：N_t = N_s·(Di_s/Di_t)^(2/3)
    var fT2 = Object.assign({}, fT, { N: Math.round(Nt) });
    var dT = fermenterDerived(fT2, s.do_set, s.mu);
    var f = Math.pow((fT.V_L || 1) / (fS.V_L || 1), 1 / 3);     // 线性放大因子
    // 比值
    var rPv = dS.Pv_kW > 0 ? dT.Pv_kW / dS.Pv_kW : 0;
    var rTip = dS.tip > 0 ? dT.tip / dS.tip : 0;
    var rKla = dS.kla_h > 0 ? dT.kla_h / dS.kla_h : 0;
    var rTmix = dS.tmix > 0 ? dT.tmix / dS.tmix : 0;
    var rOTR = dS.OTR > 0 ? dT.OTR / dS.OTR : 0;
    // 风险
    var risks = [];
    if (dT.tip > HP.shearTipThr) risks.push({ lvl: "red", t: "目标桨尖速度 " + fmt(dT.tip, 2) + " m/s > " + fmt(HP.shearTipThr, 1) + "，剪切损伤（" + (HP.key === "yeast" ? "毕赤酵母虽耐剪切仍" : "菌丝/融合蛋白") + "）风险" });
    else if (dT.tip > HP.shearTipThr * 0.8) risks.push({ lvl: "yellow", t: "目标桨尖速度 " + fmt(dT.tip, 2) + " m/s 偏高，需关注剪切" });
    if (dT.tmix > HP.tmixThr * 4) risks.push({ lvl: "red", t: "目标混合时间 " + fmt(dT.tmix, 0) + " s 过长，混合/传质不均" });
    else if (dT.tmix > HP.tmixThr * 2) risks.push({ lvl: "yellow", t: "目标混合时间 " + fmt(dT.tmix, 0) + " s 偏长" });
    if (dT.OTR < s.otr_target) risks.push({ lvl: "red", t: "目标 OTR " + fmt(dT.OTR, 0) + " < 目标需氧 " + fmt(s.otr_target, 0) + " mmol/L/h，氧限制风险" });
    else if (dT.OTR < s.otr_target * 1.3) risks.push({ lvl: "yellow", t: "目标 OTR 余量不足（" + fmt(dT.OTR, 0) + " vs " + fmt(s.otr_target, 0) + "）" });
    if (rKla < 0.6) risks.push({ lvl: "red", t: "目标 kLa 相对源下降明显，传氧能力弱化" });
    var lvl = "green";
    risks.forEach(function (r) { if (r.lvl === "red") lvl = "red"; else if (r.lvl === "yellow" && lvl !== "red") lvl = "yellow"; });
    if (!risks.length) risks.push({ lvl: "green", t: "各工程判据在推荐区间内" });
    // 推荐（氧限制时估算所需通气 / 氧气配比）
    var rec = "";
    if (dT.OTR < s.otr_target) {
      var dC = 0.008 * (1 - s.do_set / 100);
      // 空气底通（21% O₂）单靠提高 vvm 所需值
      var vvmAir = Math.pow(s.otr_target / (0.032 * Math.pow(dT.Pv_kW, 0.4) * dC * 112500), 2);
      if (vvmAir <= 1.5) {
        rec = "氧限制：建议目标尺度空气底通 vvm 提升至 ≈ " + fmt(vvmAir, 2) + "（≤1.5 VVM 内可解），并复核桨尖剪切。";
      } else {
        // 总通气封顶 1.5 VVM 时，达成所需气体 O₂ 分数（纯氧共通抬高 y_O2）
        var yO2need = Math.min(0.95, s.otr_target / (0.032 * Math.pow(dT.Pv_kW, 0.4) * Math.sqrt(1.5) * 0.008 / 0.21 * dC * 112500));
        rec = "氧限制：空气底通 1.5 VVM 仍不足，建议纯氧共通将气体 O₂ 抬至 ≈ " + fmt(yO2need * 100, 0) + "%（降低空气底通、提高纯氧占比，总通气≤1.5 VVM），并复核桨尖剪切 / 冷却余量。";
      }
    }
    var rows = [
      ["源体积 / 目标体积 (L)", fmt(fS.V_L, 0) + " / " + fmt(fT.V_L, 0)],
      ["线性放大因子 f", fmt(f, 2)],
      ["目标转速 (rpm)", fmt(Nt, 0) + "  （源 " + fmt(fS.N, 0) + "）"],
      ["P/V (kW/m³)", fmt(dS.Pv_kW, 2) + " → " + fmt(dT.Pv_kW, 2) + "  （×" + fmt(rPv, 2) + "）"],
      ["桨尖速度 (m/s)", fmt(dS.tip, 2) + " → " + fmt(dT.tip, 2) + "  （×" + fmt(rTip, 2) + "）"],
      ["kLa (1/h)", fmt(dS.kla_h, 1) + " → " + fmt(dT.kla_h, 1) + "  （×" + fmt(rKla, 2) + "）"],
      ["OTR (mmol/L/h)", fmt(dS.OTR, 0) + " → " + fmt(dT.OTR, 0) + "  （×" + fmt(rOTR, 2) + "）"],
      ["混合时间 (s)", fmt(dS.tmix, 0) + " → " + fmt(dT.tmix, 0) + "  （×" + fmt(rTmix, 2) + "）"]
    ].map(function (x) { return "<tr><td>" + x[0] + "</td><td>" + x[1] + "</td></tr>"; }).join("");
    var html = '<table class="tbl">' + rows + '</table>' +
      (rec ? '<div class="risk-yellow" style="margin:6px 0">💡 ' + esc(rec) + '</div>' : '') +
      '<div class="' + ({ red: "risk-red", yellow: "risk-yellow", green: "risk-green" }[lvl]) + '" style="font-size:18px;margin:8px 0">综合风险：' + ({ red: "🔴 高风险", yellow: "🟡 需关注", green: "🟢 可控" }[lvl]) + '</div>' +
      '<ul class="notes">' + risks.map(function (x) { return '<li class="' + ({ red: "risk-red", yellow: "risk-yellow", green: "risk-green" }[x.lvl]) + '">' + esc(x.t) + '</li>'; }).join("") + '</ul>' +
      '<div class="hint">几何相似放大：H/D=3，V=(3π/4)D³；P=Np·ρ·N³·Di⁵，P/V∝N³·Di²。等 P/V 与等 kLa（同 vvm）同源 → N_t=N_s·(Di_s/Di_t)^{2/3}；等桨尖速度 → N_t=N_s·(Di_s/Di_t)；等混合时间 → N_t=N_s。风险判据为工程代理，非工艺承诺。</div>';
    host.innerHTML = html;
    if ($("s_risk")) renderScaleRisk($("s_risk"));
  }

  // -------------------------------------------------------------------------
  // 工艺层面放大风险评估（宿主可切换：E. coli / 毕赤酵母 Pichia pastoris）
  //   覆盖：① OTR/OUR 供氧能力  ② 比生长速率 μ（氧可支撑 vs 代谢溢流阈值）
  //        ③ 碳源溢流 / 副产（E.coli→乙酸；Yeast→乙醇）  ④ CO₂ 积累（传质-吹脱平衡）
  //        ⑤ 混合  ⑥ 剪切  ⑦ 控制（冷却）
  //   均为筛选级工程估算，非工艺承诺。
  // -------------------------------------------------------------------------
  // 宿主特异性参数档案：生化常数 + 代谢/阈值默认（传氧、产热、抑制阈、剪切耐受等）
  var SCALE_HOSTS = {
    ecoli: {
      key: "ecoli", label: "🦠 大肠杆菌 E. coli",
      Y_OX: 0.015, H_CO2: 34, M_OVER: 60, DH_O2: 455,
      overflowName: "乙酸", overflowShort: "乙酸",
      qO2: 8, mu: 0.20, muCrit: 0.20, shearTipThr: 7, tmixThr: 30, overflowThr: 1.5, pco2Thr: 0.4,
      desc: "好氧呼吸 + 葡萄糖溢流 → 乙酸积累；μ_crit 为代谢溢流阈值；无细胞壁、剪切敏感；OTR/OUR 为放大首要瓶颈。"
    },
    yeast: {
      key: "yeast", label: "🍶 毕赤酵母 Pichia pastoris",
      Y_OX: 0.020, H_CO2: 30, M_OVER: 46, DH_O2: 460,
      overflowName: "乙醇", overflowShort: "乙醇",
      qO2: 5, mu: 0.15, muCrit: 0.15, shearTipThr: 10, tmixThr: 30, overflowThr: 10, pco2Thr: 0.6,
      desc: "毕赤酵母（Pichia pastoris，甲醇营养型）：甘油/葡萄糖过量→乙醇（非 Crabtree 严格型）；甲醇补料且氧受限→甲醇积累、比产率下降；有细胞壁、耐剪切；可达 400+ g DCW/L 高密度；CO₂ 抑制阈较高，放大瓶颈偏混合/氧传递（本档案针对毕赤酵母，不含酿酒酵母 S. cerevisiae）。"
    }
  };
  function scaleHostProfile(s) { return SCALE_HOSTS[(s && s.host) || "ecoli"] || SCALE_HOSTS.ecoli; }

  // 宿主可切换的工艺放大风险引擎（s.host 决定生化常数与代谢/阈值）
  function scaleProcRisk(f, s) {
    var d = fermenterDerived(f, s.do_set, s.mu);          // s.mu = 发酵液黏度 (Pa·s)
    var HP = scaleHostProfile(s);
    var P = s.proc || (s.proc = { X: 100, mu: 0.20, muCrit: 0.20, qO2: 8, overflowThr: 1.5, pco2Thr: 0.4, shearTipThr: 7, tmixThr: 30, coolU: 500, coolDT: 15, presMax: 1.5, corr: { kla: 1, tmix: 1, heat: 1 } });
    var X = P.X || 100, mu = P.mu || 0.2, muCrit = P.muCrit || 0.2, qO2 = (P.qO2 != null ? P.qO2 : HP.qO2);
    var overflowThr = (P.overflowThr != null ? P.overflowThr : HP.overflowThr), pco2Thr = (P.pco2Thr != null ? P.pco2Thr : HP.pco2Thr);
    var shearTipThr = (P.shearTipThr != null ? P.shearTipThr : HP.shearTipThr), tmixThr = (P.tmixThr != null ? P.tmixThr : HP.tmixThr), coolU = P.coolU || 500, coolDT = P.coolDT || 15;
    var corr = P.corr || { kla: 1, tmix: 1, heat: 1 };
    var klaC = corr.kla || 1, tmixC = corr.tmix || 1, heatC = corr.heat || 1;
    // 模型校正因子：kLa、混合时间、热量移除（可由 DO-stat/示踪/换热实验标定）
    d.kla_h = d.kla_h * klaC;
    d.OTR = d.OTR * klaC;                 // OTR ∝ kLa
    var tmix = d.tmix * tmixC;            // 校正后混合时间 (s)
    var Y_OX = HP.Y_OX;        // g DCW / mmol O₂
    var H_CO2 = HP.H_CO2;      // mmol/L/atm（亨利常数，随温度）
    var M_OVER = HP.M_OVER;    // g/mol（溢流副产摩尔质量：乙酸 60 / 乙醇 46）
    var DH_O2 = HP.DH_O2;      // J/mmol O₂（好氧氧化）
    // ① 供氧：OUR 需氧 vs OTR 供氧
    var OUR = qO2 * X;                         // mmol O₂/L/h
    var OTR = d.OTR;                           // mmol O₂/L/h（已含 kLa 校正）
    var otrMargin = OUR > 0 ? OTR / OUR : 0;   // OTR 相对需氧余量
    var o2def = Math.max(0, OUR - OTR);        // mmol/L/h 氧缺口
    var o2defPct = OUR > 0 ? o2def / OUR : 0;
    var muSupport = X > 0 ? OTR * Y_OX / X : 0;  // h⁻¹ 氧可支撑 μ
    // ③ 碳源溢流 / 副产（乙酸 or 乙醇）
    var muExcess = Math.max(0, mu - muCrit);
    var overflowRate_mmol = o2def * 0.6 + 2.0 * muExcess * X;
    var overflowRate = overflowRate_mmol / 1000 * M_OVER;     // g/L/h
    var t_hd = 15;
    var overflowAccum = overflowRate * t_hd;
    // ④ CO₂：RQ≈1 ⇒ CER≈OUR
    var CER = OUR;
    var kLaCO2 = d.kla_h * 0.9;
    var CO2_liq = kLaCO2 > 0 ? CER / kLaCO2 : (CER > 0 ? 1e9 : 0);
    var pCO2 = CO2_liq / H_CO2;
    // ⑤ 混合（t_mix / feed 点分散）
    var mixingMargin = tmix > 0 ? tmixThr / tmix : 0;   // >1 越充裕
    // ⑥ 剪切（桨尖 / 局部剪切应力）
    var tip = d.tip, gdot = d.gdot, tau = d.shear;
    var shearMargin = tip > 0 ? shearTipThr / tip : 0;
    // ⑦ 控制：冷却余量（A/V↓ 放大换热瓶颈）
    var Qgen = OUR * DH_O2 / 3600;             // W/L 产热
    var AoV = d.D > 0 ? 4 / d.D : 0;           // 夹套 A/V ≈ 4/D (1/m)
    var Qrem = (coolU * AoV * coolDT / 1000) * heatC;   // W/L 可移除（含校正）
    var coolMargin = Qgen > 0 ? Qrem / Qgen : 0;
    // 判级
    var risks = [];
    if (otrMargin < 1) risks.push({ lvl: "red", t: "供氧 OTR " + fmt(OTR, 0) + " < 需氧 OUR " + fmt(OUR, 0) + "（余量 " + fmt(otrMargin, 2) + "×），氧限制 → " + HP.overflowName + "溢流且 μ 被压低" });
    else if (otrMargin < 1.5) risks.push({ lvl: "yellow", t: "OTR 余量偏紧（" + fmt(otrMargin, 2) + "×，<1.5），高密期易逼近氧极限" });
    if (mu > muSupport) risks.push({ lvl: "red", t: "设定 μ " + fmt(mu, 2) + " > 氧可支撑 μ_cap " + fmt(muSupport, 3) + "，实际生长受氧限制" });
    else if (mu > muCrit) risks.push({ lvl: "yellow", t: "μ " + fmt(mu, 2) + " > μ_crit " + fmt(muCrit, 2) + "（代谢溢流阈值），有" + HP.overflowName + "形成风险" });
    if (overflowAccum > overflowThr) risks.push({ lvl: "red", t: "预计" + HP.overflowName + "累计 " + fmt(overflowAccum, 1) + " g/L > 抑制阈值 " + fmt(overflowThr, 1) + " g/L" });
    else if (overflowRate > 0.05) risks.push({ lvl: "yellow", t: HP.overflowName + "净产率 " + fmt(overflowRate, 2) + " g/L/h，存在溢流" });
    if (pCO2 > pco2Thr) risks.push({ lvl: "red", t: "稳态 pCO₂ ≈ " + fmt(pCO2, 2) + " atm > 抑制阈值 " + fmt(pco2Thr, 2) + " atm，CO₂ 抑制生长" });
    else if (pCO2 > 0.25) risks.push({ lvl: "yellow", t: "pCO₂ ≈ " + fmt(pCO2, 2) + " atm 偏高，需加强吹脱" });
    if (tmix > tmixThr) risks.push({ lvl: "red", t: "混合时间 " + fmt(tmix, 0) + " s > " + fmt(tmixThr, 0) + " s，feed 点局部底物过量（糖分脉冲→" + HP.overflowName + "）风险" });
    else if (tmix > tmixThr * 0.6) risks.push({ lvl: "yellow", t: "混合时间 " + fmt(tmix, 0) + " s 偏长，feed 分散不充分" });
    if (shearMargin < 1) risks.push({ lvl: "red", t: "桨尖速度 " + fmt(tip, 2) + " m/s > 阈值 " + fmt(shearTipThr, 1) + "，局部剪切损伤风险（Rushton 最显著）" });
    else if (shearMargin < 1.5) risks.push({ lvl: "yellow", t: "剪切余量 " + fmt(shearMargin, 2) + "× 偏紧，桨尖 " + fmt(tip, 2) + " m/s" });
    if (coolMargin < 1) risks.push({ lvl: "red", t: "冷却余量 " + fmt(coolMargin, 2) + "× 不足：产热 " + fmt(Qgen, 0) + " W/L > 可移除 " + fmt(Qrem, 0) + " W/L（A/V↓ 放大瓶颈）" });
    else if (coolMargin < 1.5) risks.push({ lvl: "yellow", t: "冷却余量 " + fmt(coolMargin, 2) + "× 偏紧，大尺度换热受限" });
    var level = "green";
    risks.forEach(function (r) { if (r.lvl === "red") level = "red"; else if (r.lvl === "yellow" && level !== "red") level = "yellow"; });
    if (!risks.length) risks.push({ lvl: "green", t: "各工艺判据在推荐窗口内" });
    return { d: d, host: HP.key, HP: HP, X: X, mu: mu, muCrit: muCrit, qO2: qO2, overflowThr: overflowThr, pco2Thr: pco2Thr,
      OUR: OUR, OTR: OTR, otrMargin: otrMargin, o2def: o2def, o2defPct: o2defPct, muSupport: muSupport,
      muExcess: muExcess, overflowRate: overflowRate, overflowAccum: overflowAccum, CER: CER, kLaCO2: kLaCO2, CO2_liq: CO2_liq, pCO2: pCO2,
      tmix: tmix, mixingMargin: mixingMargin, tip: tip, gdot: gdot, tau: tau, shearMargin: shearMargin,
      Qgen: Qgen, Qrem: Qrem, coolMargin: coolMargin,
      risks: risks, level: level };
  }

  // 氧气调节建议：给定目标需氧 OUR，求达成所需 OTR 的空气底通 / 纯氧共通配比。
  //   原则：① 默认 1 VVM 为空气底通；② DO 不足时先用纯氧共通抬升气体 O₂ 分数（y_O2↑ → C*↑ → OTR↑）；
  //        ③ 若总通气 > 1.5 VVM，则降低空气底通、提高纯氧占比（总通气封顶 1.5，保留 ≥0.3 VVM 空气用于 CO₂ 吹脱）。
  //   数值解法（二分），返回 {needed, air_vvm, o2_vvm, vvm_tot, yO2, capped, maxed, OTR_goal, achieved}
  function o2EnrichPlan(f, s, OUR) {
    var margin = 1.3, VMAX = 1.5, AIR_MIN = 0.3;
    var OTR_goal = Math.max(0, OUR) * margin;
    var air0 = (f.vvm != null && f.vvm > 0) ? f.vvm : 1.0;
    var doSet = s.do_set || 30;
    function otrWith(air, o2) { return fermenterDerived(Object.assign({}, f, { vvm: air, o2_vvm: o2 }), doSet, s.mu).OTR; }
    var achieved0 = otrWith(air0, 0);
    if (achieved0 >= OTR_goal || OUR <= 0) {
      return { needed: false, air_vvm: air0, o2_vvm: 0, vvm_tot: air0, yO2: 0.21, capped: false, maxed: false, OTR_goal: OTR_goal, achieved: achieved0 };
    }
    // ① 空气底通不变，仅加纯氧共通
    var lo = 0, hi = 4;
    for (var i = 0; i < 44; i++) { var m = (lo + hi) / 2; if (otrWith(air0, m) < OTR_goal) lo = m; else hi = m; }
    var o2a = (lo + hi) / 2, tota = air0 + o2a;
    if (tota <= VMAX) {
      return { needed: true, air_vvm: air0, o2_vvm: o2a, vvm_tot: tota, yO2: (0.21 * air0 + o2a) / tota, capped: false, maxed: false, OTR_goal: OTR_goal, achieved: otrWith(air0, o2a) };
    }
    // ② 总通气 > 1.5：封顶 1.5，降低空气底通、提高纯氧占比
    var o2max = VMAX - AIR_MIN, otrMax = otrWith(AIR_MIN, o2max);
    if (otrMax < OTR_goal) {
      return { needed: true, air_vvm: AIR_MIN, o2_vvm: o2max, vvm_tot: VMAX, yO2: (0.21 * AIR_MIN + o2max) / VMAX, capped: true, maxed: true, OTR_goal: OTR_goal, achieved: otrMax };
    }
    var aLo = AIR_MIN, aHi = air0;
    for (var j = 0; j < 44; j++) { var am = (aLo + aHi) / 2; if (otrWith(am, VMAX - am) < OTR_goal) aLo = am; else aHi = am; }
    var airB = (aLo + aHi) / 2;
    return { needed: true, air_vvm: airB, o2_vvm: VMAX - airB, vvm_tot: VMAX, yO2: (0.21 * airB + (VMAX - airB)) / VMAX, capped: true, maxed: false, OTR_goal: OTR_goal, achieved: otrWith(airB, VMAX - airB) };
  }

  // 放大策略对比：同一对源/目标尺度，按不同准则（等 P/V、等桨尖、等 kLa）求目标操作点并评估
  function scaleStrategyCompare(s) {
    var fS = s.fermenters[s.src], fT = s.fermenters[s.tgt];
    if (!fS || !fT) return [];
    var ratio = (fT.Di > 0 && fS.Di > 0) ? fS.Di / fT.Di : 1;   // Di_s / Di_t
    var P = s.proc || {};
    var X = P.X || 100, qO2 = P.qO2 || 8, OUR = qO2 * X, pco2Thr = P.pco2Thr || 0.4;
    var defs = [["constant_PV", "pv", "等 P/V"], ["constant_tip", "tip", "等桨尖速度"], ["constant_kLa", "kla", "等 kLa"]];
    return defs.map(function (dd) {
      var c = dd[1];
      var Nt = c === "tip" ? fS.N * ratio : fS.N * Math.pow(ratio, 2 / 3);   // 等 P/V / 等 kLa；等桨尖
      var fT2 = Object.assign({}, fT, { N: Math.round(Nt) });
      var r = scaleProcRisk(fT2, s);
      return { name: dd[0], label: dd[2], Nt: Nt, Pv: r.d.Pv_kW, tip: r.d.tip, kla: r.d.kla_h, OTR: r.d.OTR,
        otrMargin: r.otrMargin, pco2: r.pCO2, pco2Margin: r.pCO2 / pco2Thr, coolMargin: r.coolMargin, level: r.level };
    });
  }

  function renderScaleRisk(host) {
    if (!host) return;
    var s = state.scale, fers = s.fermenters, P = s.proc;
    var HP = scaleHostProfile(s);
    var link = s.link || null;
    var t_hd = 15;                              // 高密期时长 h（副产累计估算）
    var tags = Object.keys(fers).sort(function (a, b) { return fers[a].V_L - fers[b].V_L; });
    var tgt = fers[s.tgt] ? s.tgt : tags[0];
    if (!P.corr) P.corr = { kla: 1, tmix: 1, heat: 1 };
    function inp(lbl, key, w) {
      return '<label style="margin:0 8px 4px 0;display:inline-block">' + lbl + ' <input type="number" step="any" value="' + P[key] + '" style="width:' + (w || 66) + 'px" onchange="FIP.scaleSetProc(\'' + key + '\',parseFloat(this.value))"></label>';
    }
    function corrInp(lbl, key, w) {
      return '<label style="margin:0 6px 4px 0;display:inline-block">' + lbl + ' <input type="number" step="any" value="' + P.corr[key] + '" style="width:' + (w || 54) + 'px" onchange="FIP.scaleSetCorr(\'' + key + '\',parseFloat(this.value))"></label>';
    }
    function riskBox(title, main, sub, lvl, desc) {
      var bg = lvl === "red" ? "#2a1212" : lvl === "yellow" ? "#2a2410" : "#10241a";
      var bd = lvl === "red" ? "#a33" : lvl === "yellow" ? "#aa3" : "#3a6";
      return '<div style="flex:1;min-width:200px;background:' + bg + ';border:1px solid ' + bd + ';border-radius:10px;padding:10px;margin-bottom:6px">' +
        '<div style="font-weight:600;margin-bottom:3px">' + esc(title) + '</div>' +
        (main != null ? '<div style="font-size:20px;line-height:1.1">' + main + '</div>' : '') +
        (sub != null ? '<div class="hint" style="margin:2px 0 0">' + esc(sub) + '</div>' : '') +
        '<div class="hint" style="margin:4px 0 0">' + esc(desc) + '</div></div>';
    }
    function cl(v, red, yel) { return v < red ? ' class="risk-red"' : (v < yel ? ' class="risk-yellow"' : ''); }
    // 宿主切换条
    var hostBar = '<div style="display:flex;gap:8px;flex-wrap:wrap;align-items:center;margin-bottom:8px">' +
      '<span class="ksub">宿主类型：</span>' +
      Object.keys(SCALE_HOSTS).map(function (k) {
        var on = (s.host || "ecoli") === k;
        return '<span class="chip' + (on ? " chip-on" : "") + '" onclick="FIP.scaleSetHost(\'' + k + '\')">' + SCALE_HOSTS[k].label + '</span>';
      }).join('') +
      ' <span class="hint" style="margin-left:4px">' + esc(HP.desc) + '</span></div>';
    // Process Twin 关联条（实时读取当前选中批次的仿真峰值，可一键采纳）
    var bid = state.selBatch;
    var srcInfo = '';
    if (bid && DATA.timeseries[bid]) {
      var sts = DATA.timeseries[bid];
      var pkDCW = sts.DCW.reduce(function (a, c) { return Math.max(a, c); }, 0);
      var pkOUR = sts.OUR.reduce(function (a, c) { return Math.max(a, c); }, 0);
      var kk = ptKPI(sts);
      var ota = pkDCW > 0 ? pkOUR / pkDCW : 0;
      srcInfo = '来源 Process Twin：批次 ' + esc(bid) + ' · 峰值 DCW≈' + fmt(pkDCW, 1) + ' g/L · OUR≈' + fmt(pkOUR, 0) + ' · 比摄氧 OTA≈' + fmt(ota, 2);
    }
    var linkBar = '<div style="display:flex;gap:10px;flex-wrap:wrap;align-items:center;margin:2px 0 8px">' +
      '<button class="btn2" onclick="FIP.scaleLinkProcess()">🔗 采纳 Process Twin 参数</button>' +
      (link ? '<span class="risk-green" style="font-size:12px">✅ 已关联 ' + esc(link.batch) + '（DCW ' + fmt(link.dcw, 0) + ' · OUR ' + fmt(link.our, 0) + ' · qO₂ ' + fmt(link.ota, 2) + ' · DO ' + fmt(link.do, 0) + '%）</span> <button class="btn2" onclick="FIP.scaleUnlinkProcess()">✕ 解除</button>'
        : (srcInfo ? '<span class="hint">' + srcInfo + '</span>' : '<span class="hint">未选择 Process Twin 批次（先在 Process Twin 选批次）</span>')) +
      '</div>';
    var controls = '<div style="display:flex;gap:4px;flex-wrap:wrap;align-items:center;margin-bottom:6px">' +
      inp('Biomass X (g/L)', 'X', 72) + inp('比生长速率 μ (h⁻¹)', 'mu', 64) + inp('代谢溢流阈值 μ_crit (h⁻¹)', 'muCrit', 92) +
      inp('比摄氧率 qO₂ (mmol/g/h)', 'qO2', 84) + inp(HP.overflowName + ' 抑制 (g/L)', 'overflowThr', 64) + inp('pCO₂抑制 (atm)', 'pco2Thr', 60) +
      '</div>' +
      '<div style="display:flex;gap:4px;flex-wrap:wrap;align-items:center;margin-bottom:6px">' +
      inp('剪切阈值 tip (m/s)', 'shearTipThr', 96) + inp('混合阈值 t_mix (s)', 'tmixThr', 92) +
      inp('换热 U (W/m²K)', 'coolU', 76) + inp('ΔT冷 (K)', 'coolDT', 64) + inp('罐压上限 (bar)', 'presMax', 72) +
      '</div>' +
      '<div style="display:flex;gap:4px;flex-wrap:wrap;align-items:center;margin-bottom:6px">' +
      '<span class="ksub">模型校正因子：</span>' + corrInp('kLa×', 'kla', 54) + corrInp('t_mix×', 'tmix', 54) + corrInp('热量×', 'heat', 60) +
      '</div>';
    var rows = tags.map(function (t) {
      var r = scaleProcRisk(fers[t], s);
      var hl = (t === tgt) ? ' style="outline:2px solid #4cc9f0;outline-offset:-2px"' : '';
      return '<tr' + hl + '>' +
        '<td>' + esc(t) + '</td>' +
        '<td>' + fmt(r.OUR, 0) + '</td>' +
        '<td>' + fmt(r.OTR, 0) + '</td>' +
        '<td' + cl(r.otrMargin, 1, 1.5) + '>' + fmt(r.otrMargin, 2) + '</td>' +
        '<td' + cl(r.muSupport, r.mu, r.muCrit) + '>' + fmt(r.muSupport, 4) + '</td>' +
        '<td' + (r.overflowAccum > r.overflowThr ? ' class="risk-red"' : (r.overflowRate > 0.05 ? ' class="risk-yellow"' : '')) + '>' + fmt(r.overflowRate, 2) + '</td>' +
        '<td' + (r.pCO2 > r.pco2Thr ? ' class="risk-red"' : (r.pCO2 > 0.25 ? ' class="risk-yellow"' : '')) + '>' + fmt(r.pCO2, 2) + '</td>' +
        '<td' + cl(r.mixingMargin, 1, 1.5) + '>' + fmt(r.tmix, 0) + '</td>' +
        '<td' + cl(r.shearMargin, 1, 1.5) + '>' + fmt(r.shearMargin, 2) + '</td>' +
        '<td' + cl(r.coolMargin, 1, 1.5) + '>' + fmt(r.coolMargin, 2) + '</td>' +
        '<td class="' + (r.level === "red" ? "risk-red" : r.level === "yellow" ? "risk-yellow" : "risk-green") + '">' + ({ red: "🔴", yellow: "🟡", green: "🟢" }[r.level]) + '</td>' +
        '</tr>';
    }).join('');
    var thead = '<tr><th>尺度</th><th>OUR<br>需氧</th><th>OTR<br>供氧</th><th>OTR余量<br>×</th><th>μ氧可<br>支撑</th><th>' + HP.overflowName + '净产<br>g/L/h</th><th>pCO₂<br>atm</th><th>t_mix<br>混合s</th><th>剪切<br>余量×</th><th>冷却<br>余量×</th><th>综合</th></tr>';
    var tbl = '<div style="overflow:auto"><table class="tbl"><thead>' + thead + '</thead><tbody>' + rows + '</tbody></table></div>';
    var rT = scaleProcRisk(fers[tgt], s);
    var dT = rT.d, fT = fers[tgt];
    var gasFlow = dT.vvm_tot * fT.V_L;                    // 总通气 L/min
    var airFlow = dT.air_vvm * fT.V_L;                    // 空气底通 L/min
    var o2Flow = dT.o2_vvm * fT.V_L;                      // 纯氧共通 L/min
    var gasO2pct = dT.yO2 * 100;                          // 气体 O₂ 分数 %
    // 氧气调节建议（空气底通 + 纯氧共通；总通气 ≤1.5 VVM 约束）
    var o2plan = o2EnrichPlan(fT, s, rT.OUR);
    var hostTitle = (rT.host === "yeast") ? "毕赤酵母高密度发酵" : "E. coli 高密发酵";
    var verdict = '<div class="sub" style="margin-top:8px">🎯 目标尺度 ' + esc(tgt) + ' · ' + hostTitle + '工艺风险明细（7 维）</div>' +
      '<div style="display:flex;gap:8px;flex-wrap:wrap">' +
      riskBox('① OTR / OUR 供氧', fmt(rT.otrMargin, 2) + '×', '需氧 ' + fmt(rT.OUR, 0) + ' / 供氧 ' + fmt(rT.OTR, 0),
        rT.otrMargin < 1 ? "red" : rT.otrMargin < 1.5 ? "yellow" : "green",
        rT.otrMargin < 1 ? "供氧<需氧：氧限制，" + HP.overflowName + "溢流且 μ 被压低（纯氧共通抬升气体 O₂ 分数可增 OTR）" : rT.otrMargin < 1.5 ? "OTR 余量偏紧（<1.5×），高密期易逼近氧极限（可纯氧共通挖潜）" : "供氧充足，有余量") +
      riskBox('② μ 比生长速率', fmt(rT.mu, 2) + ' / ' + fmt(rT.muSupport, 4), '设定 / 氧可支撑 (h⁻¹)',
        rT.mu > rT.muSupport ? "red" : rT.mu > rT.muCrit ? "yellow" : "green",
        rT.mu > rT.muSupport ? "设定 μ 超过氧可支撑上限，实际 μ 被压低" : rT.mu > rT.muCrit ? "μ>μ_crit：代谢溢流风险（" + HP.overflowName + "）" : "μ 在代谢溢流阈值内") +
      riskBox('③ ' + HP.overflowName + ' / 碳源溢流', fmt(rT.overflowAccum, 1) + ' g/L', '净产率 ' + fmt(rT.overflowRate, 2) + '（阈值 ' + fmt(rT.overflowThr, 1) + '）',
        rT.overflowAccum > rT.overflowThr ? "red" : rT.overflowRate > 0.05 ? "yellow" : "green",
        "氧缺口 " + fmt(rT.o2defPct * 100, 0) + "%" + (rT.muExcess > 0 ? " · μ 超阈值 " + fmt(rT.muExcess, 2) : "") + " → " + HP.overflowName + "形成") +
      riskBox('④ CO₂ 积累', fmt(rT.pCO2, 2) + ' atm', 'CER≈' + fmt(rT.CER, 0) + '（阈值 ' + fmt(rT.pco2Thr, 2) + '）',
        rT.pCO2 > rT.pco2Thr ? "red" : rT.pCO2 > 0.25 ? "yellow" : "green",
        "溶解 CO₂≈" + fmt(rT.CO2_liq, 1) + " mmol/L，吹脱平衡 pCO₂") +
      riskBox('⑤ 混合 (t_mix)', fmt(rT.tmix, 0) + ' s', '余量 ' + fmt(rT.mixingMargin, 2) + '× / 阈值 ' + fmt(P.tmixThr || 30, 0) + 's',
        rT.mixingMargin < 1 ? "red" : rT.mixingMargin < 1.5 ? "yellow" : "green",
        rT.mixingMargin < 1 ? "混合过慢：feed 点局部底物过量（糖分脉冲→" + HP.overflowName + "）" : rT.mixingMargin < 1.5 ? "混合偏慢，feed 分散不充分" : "混合充分，feed 快速分散") +
      riskBox('⑥ 剪切 (tip/τ)', fmt(rT.tip, 2) + ' m/s', '余量 ' + fmt(rT.shearMargin, 2) + '× · τ ' + fmt(rT.tau, 3) + ' Pa',
        rT.shearMargin < 1 ? "red" : rT.shearMargin < 1.5 ? "yellow" : "green",
        rT.shearMargin < 1 ? "桨尖超阈值，局部剪切损伤（Rushton 最显著）" : rT.shearMargin < 1.5 ? "剪切余量偏紧" : (rT.host === "yeast" ? "剪切温和，毕赤酵母（有细胞壁）耐受" : "剪切温和，E. coli 耐受")) +
      riskBox('⑦ 控制 / 冷却', fmt(rT.coolMargin, 2) + '×', '产热 ' + fmt(rT.Qgen, 0) + ' / 可移除 ' + fmt(rT.Qrem, 0) + ' W/L',
        rT.coolMargin < 1 ? "red" : rT.coolMargin < 1.5 ? "yellow" : "green",
        rT.coolMargin < 1 ? "冷却不足：A/V↓ 放大换热瓶颈" : rT.coolMargin < 1.5 ? "冷却余量偏紧，大尺度受限" : "换热余量充足") +
      '</div>';
    var opWin = '<div class="method-note" style="margin-top:6px"><b>建议可行操作域（目标 ' + esc(tgt) + ' · ' + HP.overflowName + '阈值 ' + fmt(P.overflowThr != null ? P.overflowThr : HP.overflowThr, 1) + ' g/L）：</b> ' +
      'RPM≈' + fmt(fT.N, 0) + '（范围 ' + fmt(fT.n_min, 0) + '–' + fmt(fT.n_max, 0) + '） · ' +
      '空气底通≈' + fmt(fT.vvm, 2) + ' VVM（' + fmt(airFlow, 1) + ' L/min） · ' +
      '纯氧共通≈' + fmt(fT.o2_vvm || 0, 2) + ' VVM（' + fmt(o2Flow, 1) + ' L/min） · ' +
      '总通气≈' + fmt(dT.vvm_tot, 2) + ' VVM（' + fmt(gasFlow, 1) + ' L/min） · ' +
      '气体 O₂≈' + fmt(gasO2pct, 0) + '% · 压力≤' + fmt(P.presMax || 1.5, 1) + ' bar' +
      '<br><b>目标指标：</b> kLa≈' + fmt(dT.kla_h, 0) + ' h⁻¹ · OTR≈' + fmt(dT.OTR, 0) + ' mmol/L/h · OTR/OUR 余量 ' + fmt(rT.otrMargin, 2) + '× · P/V≈' + fmt(dT.Pv_kW, 2) + ' kW/m³ · tip≈' + fmt(dT.tip, 2) + ' m/s · CER≈' + fmt(rT.CER, 0) + ' mmol/L/h（RQ=1.0）</div>';
    // 氧气调节建议卡（空气底通 + 纯氧共通；总通气 ≤1.5 VVM 约束）
    var o2box = '<div class="method-note" style="margin-top:6px;border-left:3px solid ' + (o2plan.needed ? (o2plan.maxed ? '#a33' : '#aa3') : '#3a6') + '">';
    if (!o2plan.needed) {
      o2box += '<b>💨 氧气调节：</b> 当前空气底通 ' + fmt(fT.vvm, 2) + ' VVM（气体 O₂≈' + fmt(gasO2pct, 0) + '%）下 OTR 余量充足（' + fmt(rT.otrMargin, 2) + '×），无需纯氧共通；纯氧富氧可作高密期冗余。';
    } else if (o2plan.maxed) {
      o2box += '<b>💨 氧气调节（已封顶仍不足）：</b> 即使 1.5 VVM 纯氧富氧（气体 O₂≈' + fmt(o2plan.yO2 * 100, 0) + '%）仍不足覆盖 OUR ' + fmt(rT.OUR, 0) + '（达成仅 ' + fmt(o2plan.achieved, 0) + '）；需提高转速/功率（注意剪切）或增大罐径、降低峰值 DCW。建议空气底通 ' + fmt(o2plan.air_vvm, 2) + ' + 纯氧 ' + fmt(o2plan.o2_vvm, 2) + ' VVM。';
    } else if (o2plan.capped) {
      o2box += '<b>💨 氧气调节（降空气底通、提纯氧占比）：</b> 总通气封顶 1.5 VVM 时，建议空气底通 ' + fmt(o2plan.air_vvm, 2) + ' + 纯氧 ' + fmt(o2plan.o2_vvm, 2) + ' VVM（气体 O₂≈' + fmt(o2plan.yO2 * 100, 0) + '%），可使 OTR 达 ' + fmt(o2plan.achieved, 0) + '（目标 ' + fmt(o2plan.OTR_goal, 0) + '）。';
    } else {
      o2box += '<b>💨 氧气调节（纯氧共通）：</b> 空气底通保持 ' + fmt(o2plan.air_vvm, 2) + ' VVM，纯氧共通 ' + fmt(o2plan.o2_vvm, 2) + ' VVM（气体 O₂≈' + fmt(o2plan.yO2 * 100, 0) + '%），即可使 OTR 达 ' + fmt(o2plan.achieved, 0) + '（目标 ' + fmt(o2plan.OTR_goal, 0) + '）。';
    }
    o2box += '</div>';
    var verify = '<div class="method-note" style="margin-top:6px"><b>验证计划：</b> ' +
      '混合验证 — 用 pH/盐示踪法实测 t_mix，确认 feed 点无局部底物过量；' +
      'CO₂ 验证 — 实测排气 pCO₂ 轨迹，确认大尺度 CO₂ 吹脱余量，避免生长/表达抑制；' +
      '模型治理 — 版本 v1.0.0-alpha | 状态 pilot | 适用域 ' + esc(HP.label) + ' | 搅拌釜/高密度补料发酵。</div>';
    var strat = scaleStrategyCompare(s);
    var stratRows = strat.map(function (x) {
      return '<tr><td>' + x.name + '<br><span class="ksub">' + x.label + '</span></td>' +
        '<td>' + fmt(x.Nt, 0) + '</td><td>' + fmt(x.Pv, 2) + '</td><td>' + fmt(x.tip, 2) + '</td><td>' + fmt(x.kla, 0) + '</td>' +
        '<td' + cl(x.otrMargin, 1, 1.5) + '>' + fmt(x.otrMargin, 2) + '</td>' +
        '<td' + cl(x.pco2Margin, 1, 1.5) + '>' + fmt(x.pco2Margin, 2) + '</td>' +
        '<td' + cl(x.coolMargin, 1, 1.5) + '>' + fmt(x.coolMargin, 2) + '</td>' +
        '<td class="' + (x.level === "red" ? "risk-red" : x.level === "yellow" ? "risk-yellow" : "risk-green") + '">' + ({ red: "🔴", yellow: "🟡", green: "🟢" }[x.level]) + '</td></tr>';
    }).join('');
    var stratTbl = '<div class="sub" style="margin-top:8px">⚖️ 放大策略对比（源 ' + esc(s.src) + ' → 目标 ' + esc(tgt) + ' · ' + HP.overflowName + '阈值 ' + fmt(P.overflowThr != null ? P.overflowThr : HP.overflowThr, 1) + ' g/L）</div>' +
      '<div style="overflow:auto"><table class="tbl"><thead><tr><th>策略</th><th>RPM</th><th>P/V<br>kW/m³</th><th>Tip<br>m/s</th><th>kLa<br>h⁻¹</th><th>OTR/OUR<br>余量×</th><th>CO₂<br>余量×</th><th>冷却<br>余量×</th><th>整体</th></tr></thead><tbody>' + stratRows + '</tbody></table></div>';
    var stratNote = '<div class="method-note" style="margin-top:4px">' +
      '<b>策略同源说明：</b>「等 P/V」与「等 kLa」在本平台<b>等价</b>——放大时 vvm（通气/体积比）随工作体积同步缩放，而 kLa ∝ (P/V)<sup>0.4</sup>·vvm<sup>0.5</sup>，' +
      '故维持 kLa 即等价于维持 P/V，两行数值应完全一致（仅在 kla_safety 校正系数 ≠ 1.0 时才会出现差异，默认 1.0）。' +
      '真正独立的候选是「等桨尖速度」（保护剪切敏感菌株，N<sub>t</sub>=N<sub>s</sub>·(D<sub>s</sub>/D<sub>t</sub>)）。</div>';
    var ourSrc = link
      ? 'Process Twin 关联（批次 ' + esc(link.batch) + '）：峰值 DCW≈' + fmt(link.dcw, 1) + ' g/L、OTA≈' + fmt(link.ota, 2) + ' mmol/g/h、DO=' + fmt(link.do, 0) + '%；OUR = qO₂·X = ' + fmt(rT.qO2, 1) + '×' + fmt(rT.X, 0) + ' = ' + fmt(rT.OUR, 0) + ' mmol/L/h（仿真派生，非实测）'
      : '直接给定（OUR = qO₂·X = ' + fmt(rT.qO2, 1) + '×' + fmt(rT.X, 0) + ' = ' + fmt(rT.OUR, 0) + ' mmol/L/h）';
    var gov = '<div class="method-note" style="margin-top:6px">' +
      '<b>不确定性：</b> kLa 关联式基于典型搅拌釜、非实际测得，建议用 DO-stat 校准；热量衡算用经验 U·A/V，需夹套/盘管换热实验复核；混合时间基于功率数估算，需示踪验证。' +
      '<br><b>混合模型校正：</b> kLa×' + fmt(P.corr.kla, 3) + '、t_mix×' + fmt(P.corr.tmix, 3) + '、热量×' + fmt(P.corr.heat, 3) + '（在线校正因子，可由实验标定）。' +
      '<br><b>历史校准（源尺度 ' + esc(s.src) + ' · 宿主 ' + esc(HP.label) + '）：</b> ' + (link ? '已关联 Process Twin 批次 ' + esc(link.batch) + '，可据此在线校正偏置。' : '源尺度无历史批次，使用默认偏置。') +
      '<br><b>OUR 来源：</b> ' + ourSrc + '。' +
      '<br><b>审批记录：</b> 模型负责人 / 待审批 / 2026-09-11。</div>';
    var note = '<div class="method-note">模型假设（' + esc(HP.overflowName) + '溢流型）：纯好氧、RQ≈1 ⇒ CER≈OUR；OTR=kLa·ΔC，kLa=0.032·(P/V)^0.4·(总通气 vvm)^0.5，ΔC=C*·(1−DO/100)；<b>C*</b> 随气体氧分压升高：C*=0.008·y_O₂/0.21（y_O₂=(0.21·air+o2)/(air+o2)），默认 1 VVM 空气底通（y_O₂=0.21，C*≈8 mg/L），纯氧共通 o2_vvm 抬升 y_O₂ 进而提高 OTR；总通气 >1.5 VVM 时建议降空气底通、提纯氧占比。CO₂ 亨利常数 ' + fmt(HP.H_CO2, 0) + ' mmol/L/atm；Y_OX=' + fmt(HP.Y_OX, 3) + ' g DCW/mmol O₂；' + HP.overflowName + '净产率=0.6·氧缺口 + 2.0·(μ−μ_crit)·X；冷却 A/V≈4/D、产热 ' + fmt(HP.DH_O2, 0) + ' J/mmol O₂。均为筛选级估算，非工艺承诺。高密期按 ' + t_hd + ' h 估算' + HP.overflowName + '累计。</div>';
    host.innerHTML = hostBar + linkBar + '<div class="sub">🧫 ' + hostTitle + ' · 工艺放大风险评估（OTR/OUR · μ · ' + HP.overflowName + ' · CO₂ · 混合 · 剪切 · 控制）</div>' +
      controls + tbl + verdict + opWin + o2box + verify + stratTbl + stratNote + gov + note;
  }

  // Process Twin 关联：将当前选中批次的仿真峰值派生为 Scale-up 工艺参数
  function scaleLinkFromProcess() {
    var bid = state.selBatch;
    if (!bid) { alert("请先在 Process Twin 选择一个批次作为来源（左侧 Process Twin → 批次选择）。"); return; }
    var sts = DATA.timeseries[bid];
    var b = null; DATA.batches.forEach(function (x) { if (x.batch_id === bid) b = x; });
    if (!sts) { alert("该批次无时间序列仿真数据，无法关联。"); return; }
    var pkDCW = sts.DCW.reduce(function (a, c) { return Math.max(a, c); }, 0);
    var pkOUR = sts.OUR.reduce(function (a, c) { return Math.max(a, c); }, 0);
    var kk = ptKPI(sts);
    var ota = pkDCW > 0 ? pkOUR / pkDCW : 0;                 // 比摄氧率 qO₂
    var muPeak = kk ? Math.max.apply(null, kk.mu.concat([0.05])) : 0.05;
    // 宿主同步（批次 > Process Twin 过滤宿主）
    var ph = b ? ptHostType(b) : "other";
    if (ph === "yeast" || ph === "ecoli") state.scale.host = ph;
    else if (state.process.host === "yeast" || state.process.host === "ecoli") state.scale.host = state.process.host;
    var P = state.scale.proc;
    P.X = +pkDCW.toFixed(1);                 // 峰值 DCW → Biomass X
    P.qO2 = +ota.toFixed(2);                 // OTA → 比摄氧率
    P.mu = +(Math.max(muPeak, 0.05)).toFixed(2);
    state.scale.do_set = state.process.do || 30;
    state.scale.otr_target = Math.ceil(pkOUR);
    var tag = b ? b.scale_tag : null;
    if (tag && state.scale.fermenters[tag]) state.scale.src = tag;
    state.scale.link = { batch: bid, host: state.scale.host, dcw: +pkDCW.toFixed(1), ota: +ota.toFixed(2), our: +pkOUR.toFixed(0), do: state.scale.do_set, scaleTag: tag };
    if ($("s_risk")) renderScaleRisk($("s_risk"));
    if ($("s_out")) calcScale();
  }
  function scaleUnlinkProcess() {
    state.scale.link = null;
    if ($("s_risk")) renderScaleRisk($("s_risk"));
  }
  function scaleSetHost(host) {
    if (!SCALE_HOSTS[host]) return;
    state.scale.host = host;
    var HP = SCALE_HOSTS[host], P = state.scale.proc;
    // 应用宿主特异性代谢/阈值默认（保留 X、冷却、校正因子、DO/黏度等设定）
    P.qO2 = HP.qO2; P.mu = HP.mu; P.muCrit = HP.muCrit;
    P.shearTipThr = HP.shearTipThr; P.tmixThr = HP.tmixThr; P.overflowThr = HP.overflowThr; P.pco2Thr = HP.pco2Thr;
    if ($("s_risk")) renderScaleRisk($("s_risk"));
    if ($("s_out")) calcScale();
  }

  // ---- Copilot ----
  // 空状态引导（P3-⑬）：首次进入给出能力说明 + 可点击示例问题，避免「空白无引导」。
  var COPILOT_GUIDE = [
    "当前能力：基于知识库关键词检索（" + (DATA ? DATA.knowledge.length : 0) + " 条）与 " +
      (DATA ? DATA.stats.n_batches : 0) + " 个历史批次的问答顾问。问得越具体（带宿主/尺度/位号），命中越准。",
    "试试下面的问题，或直接在上方输入："
  ];
  var COPILOT_EXAMPLES = [
    "如何提高 Pichia 从 M 放大到 L 的成功率？",
    "软测量如何估计 Biomass / DCW？",
    "OUR 和 CER 分别代表什么，怎么用？",
    "大肠杆菌和酵母在放大时的关键风险差异？",
    "DO 跌破设定值通常说明什么？"
  ];
  function renderCopilot(body) {
    var chips = COPILOT_EXAMPLES.map(function (q, i) {
      return '<span class="copilot-ex" onclick="FIP.askWith(' + i + ')">' + esc(q) + '</span>';
    }).join(" ");
    var guide = COPILOT_GUIDE.map(function (g) { return esc(g); }).join("<br>");
    var html = card(
      '<div class="sub">💬 Fermentation Copilot</div>' +
      '<div class="copilot-cap">' + guide + '</div>' +
      '<textarea id="c_q" rows="3" style="width:100%" placeholder="问：如何提高 Pichia 从 M 放大到 L 的成功率？软测量如何估计 Biomass？">' + (state.copilot || "") + '</textarea>' +
      '<button class="btn" onclick="FIP.ask()">▶ 提问</button>' +
      '<div class="copilot-ex-wrap"><span class="copilot-ex-lbl">示例问题：</span>' + chips + '</div>' +
      '<div id="c_a" style="white-space:pre-wrap;margin-top:10px;line-height:1.6"></div>');
    body.innerHTML = html;
    if (state.copilot) $("c_a").textContent = copilotAnswer(state.copilot);
    else $("c_a").textContent = copilotAnswer("");  // 空状态兜底文案
  }

  // -------------------------------------------------------------------------
  // Expression ↔ Process 耦合（首试前 titer 先验 / 温融合）
  // -------------------------------------------------------------------------
  function couplingFlagColor(f) {
    if (f === "below_expression_potential") return "risk-red";
    if (f === "above_prior") return "risk-yellow";
    return "risk-green";
  }
  function couplingFlagCN(f) {
    if (f === "below_expression_potential") return "过程实现度低于表达潜力 ⚠";
    if (f === "above_prior") return "高于设计期先验";
    return "落在表达潜力带内 ✅";
  }
  function renderCoupling(body) {
    var cb = DATA.coupling || { cold_start: [], warm: [], host_prior: {}, calibration: {}, note: "" };
    var html = "";
    html += card('<div class="sub">🔗 Expression ↔ Process Twin 耦合</div>' +
      '<div class="copilot-cap">' + esc(cb.note || "Expression↔Process 耦合") + '</div>');

    // ---- 宿主级标定系数（Expression 量纲 → Process 量纲）----
    var cal = cb.calibration || {};
    var calKeys = Object.keys(cal).filter(function (k) { return cal[k] && cal[k].applied; });
    if (calKeys.length) {
      var crows = Object.keys(cal).map(function (k) {
        var c = cal[k];
        var sc = c.applied ? ("×" + fmt(c.scale, 3)) : "×1.000";
        var col = c.applied ? "#1f9e89" : "#888";
        return "<tr class='brow'><td style='color:" + hostColor(k) + "'>" + esc(k.replace(/_/g, " ")) + "</td>" +
          "<td><b style='color:" + col + "'>" + sc + "</b></td>" +
          "<td>" + (c.expression_center_g_l != null ? fmt(c.expression_center_g_l, 3) : "—") + "</td>" +
          "<td>" + (c.process_center_g_l != null ? fmt(c.process_center_g_l, 3) : "—") + "</td>" +
          "<td>" + (c.cap_g_l != null ? fmt(c.cap_g_l, 3) : "—") + "</td>" +
          "<td>" + (c.n_batches || 0) + " / " + (c.n_constructs || 0) + "</td>" +
          "<td class='ksub'>" + esc(c.note || "") + "</td></tr>";
      }).join("");
      html += card('<div class="sub">⚖ 宿主级标定系数（Expression 设计期量纲 → Process 实测量纲）</div>' +
        '<div class="hint">' + esc(cb.calibration_note || "") + '</div>' +
        '<div style="overflow:auto"><table class="tbl"><tr><th>宿主</th><th>标定系数</th><th>Expression 预测中位 (g/L)</th><th>Process 实测中位 (g/L)</th><th>截断上限 P90 (g/L)</th><th>样本 批次/构型</th><th>依据</th></tr>' + crows + '</table></div>');
    }

    // ---- 冷启动（首试前）----
    html += card('<div class="sub">❄ 冷启动：首试前 Expression 预测被 Process Twin 引用为 titer 先验</div>' +
      '<div class="hint">尚无过程历史时，Process Twin 直接引用 Expression Twin 设计期预测（source=expression_prior），不同分子第一次实验前即可给出 titer 参考。' +
      '表中「原口径」为标定前 Expression 输出，「标定后」为乘以宿主系数（并按需截断）后交给 Process Twin 的值；' +
      '排序按「未截断标定潜力」（潜力列出现时即标定后触顶被 P90 截断，潜力仍反映该构型相对位置）。</div>');
    if (!cb.cold_start.length) {
      html += card('<div class="hint">Expression Twin 暂无构型，无法派生冷启动先验。</div>');
    } else {
      var rows = cb.cold_start.map(function (c) {
        var scaled = c.calibration_scale != null && Math.abs(c.calibration_scale - 1) > 1e-9;
        var rawCell = c.expr_titer_raw != null
          ? fmt(c.expr_titer_raw, 2) + " <span class='ksub'>×</span>" + fmt(c.calibration_scale || 1, 2)
          : "—";
        // 未截断标定潜力：仅排序用，避免 P90 硬截断把多个 Top 构型压平到同一上限后丢失相对次序
        var potLine = (c.expr_potential != null && c.capped)
          ? " <div class='ksub' style='color:#2b6cb0'>潜力(未截断) " + fmt(c.expr_potential, 2) + " g/L</div>"
          : "";
        return "<tr class='brow'><td>" + esc(c.sequence_id) + "</td><td>" + esc(c.protein_name || "") + "</td>" +
          "<td style='color:" + hostColor(c.host) + "'>" + esc((c.host || "").replace(/_/g, " ")) + "</td>" +
          "<td>" + esc(c.predicted_level || "—") + "</td>" +
          "<td>" + rawCell + "</td>" +
          "<td><b style='color:#1f9e89'>" + fmt(c.expr_titer, 2) + "</b> <span class='ksub'>共形 95% PI " + fmt(c.ci[0], 2) + "–" + fmt(c.ci[1], 2) + "</span>" +
          (c.ci_tree ? " <span class='ksub'>（旧·树标准差 " + fmt(c.ci_tree[0], 2) + "–" + fmt(c.ci_tree[1], 2) + "）</span>" : "") +
          (c.capped ? " <span class='ec-badge' style='background:#f6e05e;color:#744210'>≥P90 截断</span>" : "") + potLine + "</td>" +
          "<td>" + fmt(c.confidence * 100, 0) + "%</td>" +
          "<td><span class='ec-badge ec-good'>Process Twin 引用</span></td></tr>";
      }).join("");
      html += card('<div style="overflow:auto"><table class="tbl"><tr><th>构型</th><th>蛋白</th><th>宿主</th><th>表达水平</th><th>原口径 (g/L)×系数</th><th>标定后 titer (g/L)<br><span class="ksub">按潜力↓排序</span></th><th>置信度</th><th>Process Twin</th></tr>' + rows + "</table></div>");
    }

    // ---- 温融合（已有批次）----
    html += card('<div class="sub">🔥 温融合：过程预测与 Expression 先验按成熟度加权</div>' +
      '<div class="hint">已有过程时序的批次：Process Twin 过程预测与 Expression 设计期先验融合；过程数据越薄越靠近先验（process_weight 越小）。</div>');
    if (!cb.warm.length) {
      html += card('<div class="hint">当前无已预测批次可温融合（需同时具备过程预测与同宿主 Expression 先验）。</div>');
    } else {
      var PS_TXT = { molecule: "<span class='ec-badge ec-good'>分子级 🧬</span>",
                     host_fallback: "<span class='ec-badge' style='background:#feebc8;color:#7b341e'>宿主代理</span>",
                     none: "<span class='ksub'>—</span>" };
      var wrows = cb.warm.map(function (w) {
        var pv = w.process_weight != null ? w.process_weight : 0;
        var pw = Math.round(pv * 100), iw = 100 - pw;
        var bar = "<div style='display:flex;gap:4px;align-items:center;min-width:160px'>" +
          "<div class='ec-bar'><div class='ec-fill' style='width:" + pw + "%;background:#2b6cb0'></div></div>" +
          "<span class='ksub'>过程 " + pw + "% / 先验 " + iw + "%</span></div>";
        return "<tr class='brow'><td>" + esc(w.batch_id) + "</td>" +
          "<td>" + (w.molecule ? esc(w.molecule) : "—") +
          "<div class='ksub'>" + (w.construct_id ? esc(w.construct_id) : "未关联构型") + "</div></td>" +
          "<td>" + (w.prior_titer_raw != null ? fmt(w.prior_titer_raw, 2) + " <span class='ksub'>→</span> " : "") +
          fmt(w.prior_titer, 2) +
          (w.calibration_scale != null && Math.abs(w.calibration_scale - 1) > 1e-9
            ? " <span class='ksub'>×" + fmt(w.calibration_scale, 2) + "</span>" : "") + "</td>" +
          "<td>" + fmt(w.process_titer != null ? w.process_titer : 0, 2) + "</td>" +
          "<td><b style='color:#1f9e89'>" + fmt(w.blended_mean != null ? w.blended_mean : 0, 2) + "</b> <span class='ksub'>CI " + fmt((w.ci && w.ci[0]) || 0, 2) + "–" + fmt((w.ci && w.ci[1]) || 0, 2) + "</span></td>" +
          "<td>" + (w.conformal_ci ? "<span class='ksub'>" + fmt(w.conformal_ci[0], 2) + "–" + fmt(w.conformal_ci[1], 2) + "</span> <span class='ec-badge " + (w.within_conformal_pi ? "ec-good" : "ec-warn") + "'>" + (w.within_conformal_pi ? "区间内 ✓" : "越界") + "</span>" : "—") + "</td>" +
          "<td>" + (w.maturity != null ? fmt(w.maturity, 2) : "—") + "</td>" +
          "<td>" + bar + "</td>" +
          "<td>" + (PS_TXT[w.prior_source] || "—") + "</td>" +
          "<td class='" + couplingFlagColor(w.flag) + "'>" + couplingFlagCN(w.flag) + "</td></tr>";
      }).join("");
      html += card('<div style="overflow:auto"><table class="tbl"><tr><th>批次</th><th>批次生产分子</th><th>Expression 先验 titer（原口径→标定后）</th><th>过程预测 titer</th><th>融合 titer (g/L)</th><th>表达潜力 95% PI</th><th>成熟度</th><th>权重</th><th>先验来源</th><th>偏差旗标</th></tr>' + wrows + "</table></div>");
      html += '<div class="hint">' + esc(cb.molecule_note || "") + '</div>';
      html += '<div class="hint">偏差旗标说明：<span class="risk-green">落在表达潜力带内</span>＝过程实现与同分子表达潜力一致（判据：过程 titer 在设计期潜力的 50%–150% 带内）；<span class="risk-red">过程实现度低于表达潜力</span>＝该分子欠表达，优先排查工艺（诱导/补料/溶氧）；<span class="risk-yellow">高于设计期先验</span>＝过程优于设计预期。<b>注意</b>：先验来源为「宿主代理」时用的是该宿主 Top1 构型，只作演示，不能作为分子级结论。</div>';
      html += '<div class="hint">「表达潜力 95% PI」＝ Expression Twin 经共形分位数回归（CQR）概率校准的预测区间，保外覆盖率≈名义 95%（不再是旧版 RF 树标准差那种 100% 过度覆盖）；「区间内 ✓」表示过程实现度落在该可信区间内。该区间与上方偏差旗标是<b>两件事</b>：区间管「潜力估计有多不确定」，比率带管「实现度是否达标」。</div>';
    }
    body.innerHTML = html;
  }

  // -------------------------------------------------------------------------
  // 批次数据 导入/导出（纯前端，无后端）
  // -------------------------------------------------------------------------
  var BATCH_COLS = ["batch_id", "host_organism", "scale_tag", "working_volume_l", "duration_h", "harvest_od", "titer_g_l", "success", "root_cause"];
  function ioNum(v) { var n = parseFloat(v); return isNaN(n) ? 0 : n; }
  function ioCoerceBatch(b) {
    ["working_volume_l", "duration_h", "harvest_od", "titer_g_l"].forEach(function (k) { if (b[k] !== undefined) b[k] = ioNum(b[k]); });
    if (b.success !== undefined) b.success = (b.success === true || b.success === "true" || b.success === "1" || b.success === 1);
    if (b.batch_id !== undefined) b.batch_id = String(b.batch_id);
    return b;
  }
  function ioSplitCSVLine(line) {
    var out = [], cur = "", q = false;
    for (var i = 0; i < line.length; i++) {
      var c = line[i];
      if (q) { if (c === '"') { if (line[i + 1] === '"') { cur += '"'; i++; } else q = false; } else cur += c; }
      else { if (c === '"') q = true; else if (c === ",") { out.push(cur); cur = ""; } else cur += c; }
    }
    out.push(cur); return out;
  }
  function ioParseCSV(text) {
    var lines = text.split(/\r\n|\n|\r/).filter(function (l) { return l.trim().length; });
    if (!lines.length) return [];
    var headers = ioSplitCSVLine(lines[0]).map(function (h) { return h.trim(); });
    return lines.slice(1).map(function (l) {
      var cells = ioSplitCSVLine(l), o = {};
      headers.forEach(function (h, i) { o[h] = cells[i] !== undefined ? cells[i] : ""; });
      return o;
    });
  }
  function ioBatchesToCSV() {
    var lines = [BATCH_COLS.join(",")];
    DATA.batches.forEach(function (b) {
      var row = BATCH_COLS.map(function (c) {
        var v = b[c]; if (v === undefined || v === null) v = ""; else v = String(v);
        if (/[",\n]/.test(v)) v = '"' + v.replace(/"/g, '""') + '"';
        return v;
      });
      lines.push(row.join(","));
    });
    return "﻿" + lines.join("\n");
  }
  function ioDownload(name, content, mime) {
    var blob = new Blob([content], { type: mime });
    var url = URL.createObjectURL(blob);
    var a = document.createElement("a"); a.href = url; a.download = name;
    document.body.appendChild(a); a.click(); a.remove(); URL.revokeObjectURL(url);
  }
  // 培养基配方库（本机 localStorage 持久化，保存后可直接复用）
  function recipeLibGet() { try { return JSON.parse(fipGet("fip_recipe_lib_v1") || "[]"); } catch (e) { return []; } }
  function recipeLibRender() {
    var sel = $("exp_recipe_libsel"); if (!sel) return;
    var lib = recipeLibGet();
    sel.innerHTML = '<option value="">— 选择已保存配方 —</option>' + lib.map(function (x) {
      return '<option value="' + esc(x.name) + '">' + esc(x.name) + '（' + (x.phase || "") + '）· ' + ((x.rows || []).length) + ' 组分</option>';
    }).join("");
  }
  // 培养基配方弹窗：首次调用时向 body 注入一次（不在 p_exp 内，避免重渲染丢失）
  function ensureExpRecipeModal() {
    if ($("exp_recipe_modal")) return;
    var m = document.createElement("div");
    m.id = "exp_recipe_modal";
    m.style.cssText = "display:none;position:fixed;inset:0;background:rgba(0,0,0,.55);z-index:60;align-items:center;justify-content:center";
    m.innerHTML =
      '<div style="background:#11202a;border:1px solid #24323c;border-radius:14px;padding:18px;width:600px;max-width:92vw;max-height:88vh;overflow:auto">' +
      '<div style="display:flex;justify-content:space-between;align-items:center"><div class="sub" id="exp_recipe_title">培养基配方</div>' +
      '<button class="btn2" onclick="FIP.expRecipeClose()">✕ 关闭</button></div>' +
      '<div class="ysub" style="margin:4px 0 8px">按组分填写培养基配方；可导出 .xls 在 Excel 中二次编辑</div>' +
      '<div class="sub" style="margin:10px 0 6px">配方库（保存到本机后可直接复用）</div>' +
      '<div style="display:flex;gap:8px;flex-wrap:wrap;align-items:center;margin-bottom:10px">' +
      '<input id="exp_recipe_libname" type="text" placeholder="配方名称，如 甘油复合培养基" style="flex:1;min-width:160px;background:#0c1922;border:1px solid #24323c;color:#e6edf3;border-radius:8px;padding:6px 8px">' +
      '<button class="btn2" onclick="FIP.expRecipeLibSave()">💾 存入配方库</button>' +
      '<select id="exp_recipe_libsel" class="fip-select" onchange="FIP.expRecipeLibApply(this.value)"><option value="">— 选择已保存配方 —</option></select>' +
      '<button class="btn2" onclick="FIP.expRecipeLibDel()">🗑 删除</button></div>' +
      '<table class="exp-step-tbl" id="exp_recipe_tbl" style="width:100%"></table>' +
      '<button class="btn2" style="margin-top:6px" onclick="FIP.expRecipeAddRow()">+ 添加组分</button>' +
      '<div style="margin-top:14px;display:flex;gap:8px;justify-content:flex-end;flex-wrap:wrap">' +
      '<button class="btn2" onclick="FIP.expRecipeAutoCN()">🔄 按配方自动计算 C源浓度/C:N</button>' +
      '<label class="btn2" style="cursor:pointer">📤 上传配方<input type="file" accept=".json,.csv,text/csv,application/json" onchange="FIP.expRecipeUpload(this.files[0])" style="display:none"></label>' +
      '<button class="btn2" onclick="FIP.expRecipeExport()">⬇ 导出 .xls</button>' +
      '<button class="btn" onclick="FIP.expRecipeSave()">保存并关闭</button></div></div>';
    document.body.appendChild(m);
    // 点击遮罩关闭
    m.addEventListener("click", function (ev) { if (ev.target === m) m.style.display = "none"; });
  }
  function expRecipeRowXml(cells) {
    var xml = "<Row>";
    cells.forEach(function (c) { xml += '<Cell><Data ss:Type="String">' + esc(String(c == null ? "" : c)) + '</Data></Cell>'; });
    return xml + "</Row>";
  }
  // 碳/氮质量分数（用于由配方自动推算 C:N 比）。按组分名粗估。
  function cnFractions(comp) {
    var c = (comp || "").toLowerCase();
    if (/甘油|glycer|葡萄糖|glucose|dextrose|蔗糖|麦芽糖|乳糖|果糖|fructose|糖|甲醇|methanol|淀粉|starch/.test(c)) return { C: 0.40, N: 0 };
    if (/蛋白胨|胨|peptone|酵母|yeast|大豆|soy|尿素|urea|铵|氨|ammon|硝酸|硝|玉米浆|大豆蛋白/.test(c)) return { C: 0.45, N: 0.11 };
    return { C: 0.40, N: 0.02 };
  }
  // 按组分名自动判定类别（碳源/氮源/其他），未显式标记时调用
  function expRecipeGuessCat(ph) {
    var rec = state.process.exp[ph].recipe || [];
    rec.forEach(function (r) {
      if (r.cat) return;
      var c = (r.comp || "").toLowerCase();
      if (/甘油|glycer|葡萄糖|glucose|dextrose|蔗糖|麦芽糖|乳糖|果糖|fructose|糖|甲醇|methanol|淀粉|starch/.test(c)) r.cat = "C";
      else if (/蛋白胨|胨|peptone|酵母|yeast|大豆|soy|尿素|urea|铵|氨|ammon|硝酸|硝|玉米浆|大豆蛋白/.test(c)) r.cat = "N";
      else r.cat = "O";
    });
  }
  // 补料策略库（本机 localStorage）
  function feedStratLibGet() { try { return JSON.parse(fipGet("fip_feed_lib_v1") || "[]"); } catch (e) { return []; } }
  function feedStratLibRender() {
    EXP_PHASES.forEach(function (ph) {
      var sel = $("feed_strat_sel_" + ph.k); if (!sel) return;
      var lib = feedStratLibGet();
      sel.innerHTML = '<option value="">— 选择已保存策略 —</option>' + lib.map(function (x) {
        var mode = ({ none: "无", constant: "恒定", stepwise: "阶梯" }[x.strat.feed_mode] || x.strat.feed_mode);
        return '<option value="' + esc(x.name) + '">' + esc(x.name) + ' · ' + mode + '</option>';
      }).join("");
    });
  }
  function applyFeedStrat(k, s) {
    if (!s || !s.feed_mode) return;
    var p = state.process.exp[k];
    p.feed_mode = s.feed_mode;
    p.feed_const = (s.feed_const != null) ? s.feed_const : p.feed_const;
    p.feed_steps = (s.feed_steps || []).map(function (st) { return { t: st.t, rate: st.rate }; });
  }

  // -------------------------------------------------------------------------
  // What-if 动态仿真 · 机理前向模型动画（与 Process / Scale-up Twin 同式）
  // -------------------------------------------------------------------------
  function renderSimAnim(body) {
    if (state.sim && state.sim.raf) { try { cancelAnimationFrame(state.sim.raf); } catch (e) {} state.sim.raf = null; }
    state.sim = state.sim || {};
    state.sim._body = body;
    var D = (DATA && DATA.sim) || {};
    if (!state.sim.host) state.sim.host = "ecoli";
    if (!state.sim.medium) state.sim.medium = (FermentSim.DEFAULT_MEDIUM && FermentSim.DEFAULT_MEDIUM[state.sim.host]) || "defined";
    if (!state.sim.speedMin) state.sim.speedMin = 4;        // 默认 24h 周期在 4 分钟内播完
    if (state.sim.baseline == null) state.sim.baseline = true;
    if (!state.sim.bubbles) state.sim.bubbles = [];
    // 阶段 FBA 耦合状态（默认关闭：手调默认值，不继承 FBA μ_het）
    if (state.sim.carbon_source == null) state.sim.carbon_source = "glucose";
    if (state.sim.strain_factor == null) state.sim.strain_factor = 1.0;
    if (state.sim.fba_coupled == null) state.sim.fba_coupled = false;
    if (state.sim.fba_temp == null) state.sim.fba_temp = 37.0;
    if (state.sim.expression_level == null) state.sim.expression_level = 0.15;
    // Pichia 表达策略（默认诱导型·纯甲醇，与原行为一致）
    if (state.sim.expression_mode == null) state.sim.expression_mode = "inducible";
    if (state.sim.meoh_feed_mode == null) state.sim.meoh_feed_mode = "methanol";
    if (state.sim.mixed_gly_feed_g_l == null) state.sim.mixed_gly_feed_g_l = 200.0;

    function clamp(v, a, b) { return v < a ? a : (v > b ? b : v); }
    function lerp(a, b, t) { return a + (b - a) * t; }
    function lerpColor(c0, c1, t) {
      return Math.round(lerp(c0[0], c1[0], t)) + "," + Math.round(lerp(c0[1], c1[1], t)) + "," + Math.round(lerp(c0[2], c1[2], t));
    }
    function roundRect(ctx, x, y, w, h, r) {
      ctx.beginPath(); ctx.moveTo(x + r, y);
      ctx.arcTo(x + w, y, x + w, y + h, r); ctx.arcTo(x + w, y + h, x, y + h, r);
      ctx.arcTo(x, y + h, x, y, r); ctx.arcTo(x, y, x + w, y, r); ctx.closePath();
    }
    var SERIES = {
      ecoli: [
        { key: "X", label: "菌浓 X", unit: "g/L", color: "#1f9e89" },
        { key: "P", label: "效价 P", unit: "g/L", color: "#80ed99" },
        { key: "DO", label: "DO", unit: "%", color: "#4cc9f0" },
        { key: "A", label: "乙酸", unit: "g/L", color: "#ff6b6b" }
      ],
      pichia: [
        { key: "X", label: "菌浓 X", unit: "g/L", color: "#1f9e89" },
        { key: "P", label: "效价 P", unit: "g/L", color: "#80ed99" },
        { key: "DO", label: "DO", unit: "%", color: "#4cc9f0" },
        { key: "M", label: "甲醇 M", unit: "g/L", color: "#ffce4d" }
      ]
    };

    function recipeFromUI() {
      if (state.sim && state.sim.expRecipe) return state.sim.expRecipe;
      var rec = FermentSim.defaultRecipe(state.sim.host, state.sim.medium);
      (FermentSim.paramRanges[state.sim.host] || []).forEach(function (pr) {
        var el = $("sim_" + pr.key); if (el) rec[pr.key] = parseFloat(el.value);
      });
      rec.medium = state.sim.medium;
      return rec;
    }
    function recompute(reset) {
      var rec;
      if (state.sim.loadedRecipe) {
        // 外部加载（历史批次 / 实验设计）：以加载 recipe 为 base，滑块字段用当前 UI 值覆盖，
        // 体积 / 诱导时刻 / 补料浓度等非滑块字段得以保留。
        rec = Object.assign({}, state.sim.loadedRecipe);
        (FermentSim.paramRanges[state.sim.host] || []).forEach(function (pr) {
          var el = $("sim_" + pr.key); if (el) rec[pr.key] = parseFloat(el.value);
        });
        rec.medium = state.sim.medium;
      } else {
        rec = recipeFromUI();
      }
      // 阶段 FBA 耦合：把界面 FBA 控件写入 recipe（E. coli 与 Pichia 均支持）
      if (state.sim.host === "ecoli") {
        rec.carbon_source = state.sim.carbon_source || "glucose";
        rec.strain_factor = (state.sim.strain_factor != null) ? state.sim.strain_factor : 1.0;
        rec.fba_coupled = !!state.sim.fba_coupled;
        rec.fba_temp = (state.sim.fba_temp != null) ? state.sim.fba_temp : 37.0;
        rec.expression_level = (state.sim.expression_level != null) ? state.sim.expression_level : 0.15;
        if (rec.fba_coupled) rec.expression_form = state.sim.form || null;  // FBA 查表需表达形式
      }
      // Pichia 表达策略 + FBA 耦合字段写入 recipe
      if (state.sim.host === "pichia") {
        rec.expression_mode = state.sim.expression_mode || "inducible";
        rec.meoh_feed_mode = state.sim.meoh_feed_mode || "methanol";
        rec.mixed_gly_feed_g_l = (state.sim.mixed_gly_feed_g_l != null) ? state.sim.mixed_gly_feed_g_l : 200.0;
        rec.strain_factor = (state.sim.strain_factor != null) ? state.sim.strain_factor : 1.0;
        rec.fba_coupled = !!state.sim.fba_coupled;
        rec.fba_temp = (state.sim.fba_temp != null) ? state.sim.fba_temp : 28.0;
        rec.expression_level = (state.sim.expression_level != null) ? state.sim.expression_level : 0.15;
        rec.expression_form = state.sim.form || null;  // FBA 查表需表达形式
      }
      state.sim.recipe = rec;
      state.sim.traj = FermentSim.simulate(rec);
      state.sim.baseTraj = (state.sim.baseline && D.baseline && D.baseline[state.sim.host]) ? D.baseline[state.sim.host] : null;
      // 阶段 FBA 耦合：E. coli 专用「半接通 vs 全接通」对比（同一配方跑两种模式）
      if (state.sim.host === "ecoli") {
        var recSemi = Object.assign({}, rec); recSemi.fba_coupled = false;
        var recFull = Object.assign({}, rec); recFull.fba_coupled = true; recFull.expression_form = state.sim.form || null;
        state.sim.trajSemi = FermentSim.simulate(recSemi);
        state.sim.trajFull = FermentSim.simulate(recFull);
        drawFbaCompare();
      }
      // 阶段 FBA 耦合：Pichia 专用「半接通 vs 全接通」对比（同一配方跑两种模式）
      if (state.sim.host === "pichia") {
        var recSemiP = Object.assign({}, rec); recSemiP.fba_coupled = false;
        var recFullP = Object.assign({}, rec); recFullP.fba_coupled = true; recFullP.expression_form = state.sim.form || null;
        state.sim.trajSemi = FermentSim.simulate(recSemiP);
        state.sim.trajFull = FermentSim.simulate(recFullP);
        drawFbaComparePichia();
      }
      if (reset || state.sim.simTime == null) state.sim.simTime = 0;
      if (!state.sim.playing) drawAll(state.sim.traj.n - 1);
      updatePlayBtn();
    }

    function drawAll(idx) {
      drawBioreactor(idx);
      drawChart(idx);
      updateKPI(idx);
      updateProgress(idx);
    }
    function updateProgress(idx) {
      var tr = state.sim.traj; if (!$("sim_progress")) return;
      var t = tr.t[idx], end = tr.t[tr.n - 1];
      $("sim_progress").textContent = "t = " + t.toFixed(1) + " h / " + end.toFixed(0) + " h · 帧 " + (idx + 1) + "/" + tr.n;
    }
    function updateKPI(idx) {
      if (!$("sim_kpi")) return;
      var tr = state.sim.traj, s = tr.summary;
      var cards = [
        ["当前菌浓 X", fmt(tr.X[idx], 2), "g/L", "#1f9e89"],
        ["当前效价 P", fmt(tr.P[idx], 2), "g/L", "#80ed99"],
        ["可溶产物 P_sol", fmt(tr.Psol[idx], 2), "g/L", "#4cc9f0"],
        ["包涵体 P_ib", fmt(tr.Pib[idx], 2), "g/L", "#ffce4d"],
        ["当前 DO", fmt(tr.DO[idx], 1), "%", tr.DO[idx] >= 25 ? "#4cc9f0" : "#ff6b6b"],
        ["终点效价预测", fmt(s.titer, 2), "g/L", "#e6edf3"],
        ["终点菌浓 DCW", fmt(s.dcw, 2), "g/L", "#9fb3c8"],
        ["乙酸峰值", fmt(s.acetatePeak, 2), "g/L", "#ff8c42"]
      ];
      $("sim_kpi").innerHTML = '<div class="kpi-grid">' + cards.map(function (c) {
        return '<div class="kpi-card"><div class="kpi-val" style="color:' + c[3] + '">' + c[1] + ' <small>' + c[2] + '</small></div><div class="kpi-lab">' + c[0] + "</div></div>";
      }).join("") + "</div>";
    }

    function drawBioreactor(idx) {
      var c = $("sim_bioreactor"); if (!c) return;
      var ctx = c.getContext("2d"); var W = c.width, H = c.height;
      ctx.clearRect(0, 0, W, H);
      var tr = state.sim.traj;
      var X = tr.X[idx], DO = tr.DO[idx], air = tr.AIR[idx], rpm = tr.RPM[idx], feed = tr.FEED[idx], induced = tr.INDUCED[idx], t = tr.t[idx];
      var tx = 48, ty = 34, tw = 224, th = 286;
      // 罐体
      ctx.fillStyle = "#0e1518"; roundRect(ctx, tx, ty, tw, th, 14); ctx.fill();
      var Xref = state.sim.host === "ecoli" ? 50 : 30;
      var ratio = clamp(X / Xref, 0, 1);
      var broth = lerpColor([12, 32, 36], [31, 158, 137], ratio); // 清 → 浓
      var fillH = th - 10;
      ctx.save(); roundRect(ctx, tx + 3, ty + 3, tw - 6, th - 6, 12); ctx.clip();
      ctx.fillStyle = "rgb(" + broth + ")";
      ctx.fillRect(tx + 3, ty + th - 3 - fillH, tw - 6, fillH);
      // 气泡（速率 ∝ 通气量）
      var now = (typeof performance !== "undefined" ? performance.now() : Date.now());
      if (state.sim.playing) {
        var spawn = Math.max(0, Math.round((air / 10) * 1.2));
        for (var s2 = 0; s2 < spawn; s2++) state.sim.bubbles.push({ x: tx + 8 + Math.random() * (tw - 16), y: ty + th - 6, r: 1.5 + Math.random() * 2.5, sp: 0.6 + Math.random() * 1.2 });
      }
      var nb = [];
      for (var b = 0; b < state.sim.bubbles.length; b++) {
        var bb = state.sim.bubbles[b]; bb.y -= bb.sp;
        if (bb.y > ty + 6) { nb.push(bb); ctx.strokeStyle = "rgba(220,245,255,0.5)"; ctx.beginPath(); ctx.arc(bb.x, bb.y, bb.r, 0, 6.283); ctx.stroke(); }
      }
      state.sim.bubbles = nb;
      // 搅拌桨（角速度 ∝ 转速）
      var ang = now * 0.001 * (rpm / 60) * 6.283;
      ctx.strokeStyle = "rgba(230,237,243,0.7)"; ctx.lineWidth = 3;
      var cx = tx + tw / 2, cy = ty + th - 28;
      ctx.beginPath(); ctx.moveTo(cx - 34 * Math.cos(ang), cy - 34 * Math.sin(ang)); ctx.lineTo(cx + 34 * Math.cos(ang), cy + 34 * Math.sin(ang)); ctx.stroke();
      ctx.lineWidth = 1;
      ctx.restore();
      ctx.strokeStyle = "#1f9e8933"; roundRect(ctx, tx, ty, tw, th, 14); ctx.stroke();
      // DO 竖条（右）
      var gx = tx + tw + 14, gw = 12;
      ctx.fillStyle = "#16242c"; roundRect(ctx, gx, ty, gw, th, 6); ctx.fill();
      var f = clamp(DO / 100, 0, 1);
      ctx.fillStyle = DO >= 25 ? "#3ddc97" : "#ff6b6b";
      roundRect(ctx, gx, ty + th * (1 - f), gw, th * f, 6); ctx.fill();
      ctx.fillStyle = "#7d93a6"; ctx.font = "10px sans-serif"; ctx.fillText("DO", gx - 2, ty - 6);
      // 文字叠层
      ctx.fillStyle = "#e6edf3"; ctx.font = "bold 13px sans-serif";
      ctx.fillText("t = " + t.toFixed(1) + " h", tx, ty - 12);
      ctx.font = "12px sans-serif"; ctx.fillStyle = "#cdd9e5";
      ctx.fillText("X=" + fmt(X, 1) + "  P=" + fmt(tr.P[idx], 2), tx + 4, ty + th - 8);
      // 补料脉冲（顶部）
      if (feed > 0) {
        var a = 0.4 + 0.6 * Math.abs(Math.sin(now * 0.006));
        ctx.fillStyle = "rgba(255,107,107," + a.toFixed(2) + ")";
        ctx.beginPath(); ctx.arc(tx + 14, ty + 12, 5, 0, 6.283); ctx.fill();
      }
      // 诱导徽标
      if (induced) {
        ctx.fillStyle = "#80ed99"; ctx.font = "bold 11px sans-serif";
        ctx.fillText("● 诱导 ON", tx + tw - 78, ty + 14);
      }
    }

    function drawChart(idx) {
      var c = $("sim_chart"); if (!c) return;
      var ctx = c.getContext("2d"); var W = c.width, H = c.height;
      ctx.clearRect(0, 0, W, H);
      var tr = state.sim.traj, base = state.sim.baseTraj;
      var ser = SERIES[state.sim.host];
      var pad = 34, n = tr.n;
      // 每序列独立量程（与基线共用最大值，便于对照）
      var maxv = {};
      ser.forEach(function (s) {
        var mx = 1e-9;
        for (var i = 0; i < n; i++) mx = Math.max(mx, tr[s.key][i]);
        if (base && base[s.key]) for (var j = 0; j < base[s.key].length; j++) mx = Math.max(mx, base[s.key][j]);
        maxv[s.key] = mx;
      });
      var end = tr.t[n - 1];
      function Xpx(i) { return pad + (i / (n - 1)) * (W - pad - 10); }
      function Ypx(v, m) { return H - pad - (v / m) * (H - pad - 14); }
      // 网格 + 时间轴
      ctx.strokeStyle = "#16242c"; ctx.fillStyle = "#6b7d8f"; ctx.font = "10px sans-serif"; ctx.lineWidth = 1;
      for (var g = 0; g <= 4; g++) {
        var yy = pad + g * (H - pad - 14) / 4;
        ctx.beginPath(); ctx.moveTo(pad, yy); ctx.lineTo(W - 10, yy); ctx.stroke();
        ctx.fillText(((1 - g / 4) * 100).toFixed(0) + "%", 6, yy + 3);
      }
      [0, 0.5, 1].forEach(function (f) {
        var xx = pad + f * (W - pad - 10); ctx.fillText((end * f).toFixed(0) + "h", xx - 6, H - 14);
      });
      // 先画基线（虚线），再画情景（实线至 idx）
      ser.forEach(function (s) {
        if (base && base[s.key]) {
          ctx.strokeStyle = s.color; ctx.globalAlpha = 0.35; ctx.setLineDash([5, 4]); ctx.lineWidth = 1.5;
          ctx.beginPath();
          for (var i = 0; i < base[s.key].length; i++) { var xx = Xpx(i / (base[s.key].length - 1) * (n - 1)); var yy = Ypx(base[s.key][i], maxv[s.key]); i ? ctx.lineTo(xx, yy) : ctx.moveTo(xx, yy); }
          ctx.stroke(); ctx.setLineDash([]); ctx.globalAlpha = 1;
        }
        ctx.strokeStyle = s.color; ctx.lineWidth = 2; ctx.beginPath();
        for (var i2 = 0; i2 <= idx; i2++) { var x2 = Xpx(i2); var y2 = Ypx(tr[s.key][i2], maxv[s.key]); i2 ? ctx.lineTo(x2, y2) : ctx.moveTo(x2, y2); }
        ctx.stroke();
      });
      // 图例（含当前值）
      var lx = pad, ly = 12;
      ser.forEach(function (s) {
        var v = tr[s.key][idx];
        ctx.fillStyle = s.color; ctx.fillRect(lx, ly - 8, 12, 3);
        ctx.fillStyle = "#cdd9e5"; ctx.font = "11px sans-serif";
        ctx.fillText(s.label + " " + fmt(v, s.key === "DO" ? 1 : 2) + s.unit, lx + 16, ly);
        lx += 120;
      });
      ctx.fillStyle = "#6b7d8f"; ctx.fillText(state.sim.baseline ? "虚线=基线(默认配方) · 实线=当前情景" : "实线=当前情景", pad, H - 2);
    }

    // FBA 半接通 vs 全接通 并排对比图（E. coli 专用）：同一配方下跑两种模式，对比 X(t)/μ(t) 与关键 KPI
    function drawFbaCompare() {
      var host = $("sim_fba_compare"); if (!host) return;
      var semi = state.sim.trajSemi, full = state.sim.trajFull;
      if (!semi || !full) { host.innerHTML = ""; return; }
      var rec = state.sim.recipe || {};
      var end = semi.t[semi.n - 1];
      var C_SEMI = "#7d93a6", C_FULL = "#1f9e89";   // 半接通=灰，全接通=青
      var activeFull = !!state.sim.fba_coupled;
      function muPts(tr) {
        var dt = tr.dt, X = tr.X, pts = [];
        for (var i = 0; i < X.length; i++) {
          var mu = (i === 0 || X[i - 1] <= 0) ? 0 : (X[i] - X[i - 1]) / (dt * X[i - 1]);
          pts.push([tr.t[i], mu]);
        }
        return pts;
      }
      function peakMu(tr) { var m = 0, p = muPts(tr); for (var i = 0; i < p.length; i++) if (p[i][1] > m) m = p[i][1]; return m; }
      function phases() {
        return [
          { start: 0, end: (rec.feed_start_h || 7), label: "Batch", color: "#4cc9f0" },
          { start: (rec.feed_start_h || 7), end: (rec.induction_h || 16), label: "Fed-batch", color: "#80ed99" },
          { start: (rec.induction_h || 16), end: end, label: "Induction", color: "#ffce4d" }
        ];
      }
      var cs = state.sim.carbon_source || "glucose";
      var csName = (cs === "glycerol") ? "甘油" : "葡萄糖";
      host.innerHTML =
        '<div class="sub">⚛ FBA 半接通 vs 全接通 · 同一配方（' + csName + '）下对比</div>' +
        '<div class="hint" style="margin-bottom:8px">半接通（默认）：动态层用 <b>手调 μ_max=0.55</b>、现象学乙酸，<b>不继承</b> FBA 的 μ_het。' +
        ' 全接通（勾选启用）：继承 FBA 的 μ_het / 宿主摄取 / 乙酸 yield。当前活动态：<b style="color:' +
        (activeFull ? C_FULL : C_SEMI) + '">' + (activeFull ? "全接通" : "半接通") + "</b>。</div>" +
        '<div style="display:flex;gap:14px;flex-wrap:wrap">' +
          '<div id="fba_cmp_x" style="flex:1;min-width:320px"></div>' +
          '<div id="fba_cmp_mu" style="flex:1;min-width:320px"></div>' +
        "</div>" +
        '<div id="fba_cmp_kpi" style="margin-top:6px"></div>';
      lineChart($("fba_cmp_x"), [
        { name: "半接通 (μ_max=0.55)", color: C_SEMI, points: semi.t.map(function (t, i) { return [t, semi.X[i]]; }) },
        { name: "全接通 (继承 μ_het)", color: C_FULL, points: full.t.map(function (t, i) { return [t, full.X[i]]; }) }
      ], { xtitle: "时间 (h)", ytitle: "菌浓 X (g/L)", phases: phases() });
      lineChart($("fba_cmp_mu"), [
        { name: "半接通 (μ_max=0.55)", color: C_SEMI, points: muPts(semi) },
        { name: "全接通 (继承 μ_het)", color: C_FULL, points: muPts(full) }
      ], { xtitle: "时间 (h)", ytitle: "比生长速率 μ (1/h)", phases: phases(), ymin: 0 });
      // KPI 对比
      var sMu = peakMu(semi), fMu = peakMu(full);
      var sX = semi.summary.dcw, fX = full.summary.dcw;
      var sA = semi.summary.acetatePeak, fA = full.summary.acetatePeak;
      function dlt(a, b) { if (b === 0) return "—"; var p = (a - b) / b * 100; return (p >= 0 ? "+" : "") + p.toFixed(1) + "%"; }
      var rows = [
        ["峰值 μ", sMu, fMu, "1/h", fMu > sMu],
        ["终点菌浓 X", sX, fX, "g/L", fX > sX],
        ["乙酸峰值", sA, fA, "g/L", fA > sA]
      ];
      var kpi = '<div style="display:flex;gap:12px;flex-wrap:wrap">';
      rows.forEach(function (r) {
        var up = r[4];
        var col = up ? C_FULL : C_SEMI;
        kpi += '<div class="predbox" style="flex:1;min-width:160px;padding:10px 14px">' +
          '<div class="pl" style="color:#9fb3c8">' + r[0] + '</div>' +
          '<div style="font-size:12px;color:#7d93a6;margin-top:4px">半接通 <b style="color:' + C_SEMI + '">' + fmt(r[1], 3) + '</b> → 全接通 <b style="color:' + C_FULL + '">' + fmt(r[2], 3) + '</b></div>' +
          '<div class="pv" style="color:' + col + ';font-size:18px;margin-top:2px">' + dlt(r[2], r[1]) + '</div>' +
          '<div class="pl" style="color:#6b7d8f">' + r[3] + '</div></div>';
      });
      kpi += "</div>";
      $("fba_cmp_kpi").innerHTML = kpi;
    }

    // FBA 半接通 vs 全接通 并排对比图（Pichia 专用）：同一配方下跑两种模式，对比 X(t)/μ(t) 与关键 KPI
    function drawFbaComparePichia() {
      var host = $("sim_fba_compare_pichia"); if (!host) return;
      var semi = state.sim.trajSemi, full = state.sim.trajFull;
      if (!semi || !full) { host.innerHTML = ""; return; }
      var rec = state.sim.recipe || {};
      var end = semi.t[semi.n - 1];
      var C_SEMI = "#7d93a6", C_FULL = "#b06ab3";   // 半接通=灰，全接通=紫
      var activeFull = !!state.sim.fba_coupled;
      function muPts(tr) {
        var dt = tr.dt, X = tr.X, pts = [];
        for (var i = 0; i < X.length; i++) {
          var mu = (i === 0 || X[i - 1] <= 0) ? 0 : (X[i] - X[i - 1]) / (dt * X[i - 1]);
          pts.push([tr.t[i], mu]);
        }
        return pts;
      }
      function peakMu(tr) { var m = 0, p = muPts(tr); for (var i = 0; i < p.length; i++) if (p[i][1] > m) m = p[i][1]; return m; }
      function pichiaPhases() {
        var ms = rec.meoh_start_h || 34, gfe = rec.gly_feed_end_h || 32;
        var gEnd = Math.min(gfe, ms);
        return [
          { start: 0, end: gEnd, label: "Glycerol", color: "#4cc9f0" },
          { start: gEnd, end: ms, label: "Starvation", color: "#ffce4d" },
          { start: ms, end: end, label: (rec.meoh_feed_mode === "mixed" ? "Mixed feed" : "MeOH induce"), color: "#80ed99" }
        ];
      }
      var fbaName = (rec.expression_mode === "constitutive") ? "甘油" : "甲醇";
      host.innerHTML =
        '<div class="sub">⚛ FBA 半接通 vs 全接通 · Pichia（' + fbaName + ' 校准）</div>' +
        '<div class="hint" style="margin-bottom:8px">半接通（默认）：动态层用 <b>手调 μ_max</b>、现象学表达负担，<b>不继承</b> FBA 的 μ_het。' +
        ' 全接通（勾选启用）：继承 iPP668 平行 biomass 的 μ_het 降幅 drop_frac（与温度/菌株无关，≡1/(1+f·FORM_FACTOR)）。当前活动态：<b style="color:' +
        (activeFull ? C_FULL : C_SEMI) + '">' + (activeFull ? "全接通" : "半接通") + "</b>。</div>" +
        '<div style="display:flex;gap:14px;flex-wrap:wrap">' +
          '<div id="fba_cmp_p_x" style="flex:1;min-width:320px"></div>' +
          '<div id="fba_cmp_p_mu" style="flex:1;min-width:320px"></div>' +
        "</div>" +
        '<div id="fba_cmp_p_kpi" style="margin-top:6px"></div>';
      lineChart($("fba_cmp_p_x"), [
        { name: "半接通 (手调 μ_max)", color: C_SEMI, points: semi.t.map(function (t, i) { return [t, semi.X[i]]; }) },
        { name: "全接通 (继承 μ_het)", color: C_FULL, points: full.t.map(function (t, i) { return [t, full.X[i]]; }) }
      ], { xtitle: "时间 (h)", ytitle: "菌浓 X (g/L)", phases: pichiaPhases() });
      lineChart($("fba_cmp_p_mu"), [
        { name: "半接通 (手调 μ_max)", color: C_SEMI, points: muPts(semi) },
        { name: "全接通 (继承 μ_het)", color: C_FULL, points: muPts(full) }
      ], { xtitle: "时间 (h)", ytitle: "比生长速率 μ (1/h)", phases: pichiaPhases(), ymin: 0 });
      // KPI 对比
      var sMu = peakMu(semi), fMu = peakMu(full);
      var sX = semi.summary.dcw, fX = full.summary.dcw;
      var sP = semi.summary.titer, fP = full.summary.titer;
      function dlt(a, b) { if (b === 0) return "—"; var p = (a - b) / b * 100; return (p >= 0 ? "+" : "") + p.toFixed(1) + "%"; }
      var rows = [
        ["峰值 μ", sMu, fMu, "1/h", fMu > sMu],
        ["终点菌浓 X", sX, fX, "g/L", fX > sX],
        ["终点效价 P", sP, fP, "g/L", fP > sP]
      ];
      var kpi = '<div style="display:flex;gap:12px;flex-wrap:wrap">';
      rows.forEach(function (r) {
        var up = r[4], col = up ? C_FULL : C_SEMI;
        kpi += '<div class="predbox" style="flex:1;min-width:160px;padding:10px 14px">' +
          '<div class="pl" style="color:#9fb3c8">' + r[0] + '</div>' +
          '<div style="font-size:12px;color:#7d93a6;margin-top:4px">半接通 <b style="color:' + C_SEMI + '">' + fmt(r[1], 3) + '</b> → 全接通 <b style="color:' + C_FULL + '">' + fmt(r[2], 3) + '</b></div>' +
          '<div class="pv" style="color:' + col + ';font-size:18px;margin-top:2px">' + dlt(r[2], r[1]) + '</div>' +
          '<div class="pl" style="color:#6b7d8f">' + r[3] + '</div></div>';
      });
      kpi += "</div>";
      $("fba_cmp_p_kpi").innerHTML = kpi;
    }

    function tick(ts) {
      if (state.page !== "simanim") { state.sim.raf = null; return; }
      if (!state.sim.playing) { state.sim.raf = null; return; }
      if (state.sim.lastTs == null) state.sim.lastTs = ts;
      var dtReal = (ts - state.sim.lastTs) / 1000; state.sim.lastTs = ts;
      var tr = state.sim.traj, totalH = tr.t[tr.n - 1], dt = tr.t[1] - tr.t[0];
      var speedH = totalH / (state.sim.speedMin * 60);
      state.sim.simTime += dtReal * speedH;
      if (state.sim.simTime >= totalH) { state.sim.simTime = totalH; state.sim.playing = false; }
      var idx = Math.min(tr.n - 1, Math.max(0, Math.floor(state.sim.simTime / dt)));
      drawAll(idx);
      if (state.sim.playing) state.sim.raf = requestAnimationFrame(tick);
      else updatePlayBtn();
    }

    function startLoop() {
      if (state.sim.raf) { try { cancelAnimationFrame(state.sim.raf); } catch (e) {} }
      state.sim.lastTs = null; state.sim.playing = true; updatePlayBtn();
      state.sim.raf = requestAnimationFrame(tick);
    }
    function updatePlayBtn() {
      var pb = $("sim_playbtn"); if (pb) pb.textContent = state.sim.playing ? "▶ 播放中…" : "▶ 播放";
    }

    // ---- 骨架 ----
    var ranges = FermentSim.paramRanges[state.sim.host] || [];
    var _defRec = FermentSim.defaultRecipe(state.sim.host, state.sim.medium);  // 含培养基默认的配方
    // 培养基备选项（按宿主）
    var _medGrp = (FermentSim.MEDIA && FermentSim.MEDIA[state.sim.host]) || {};
    var _medOpts = Object.keys(_medGrp).map(function (mk) {
      return '<option value="' + mk + '"' + (state.sim.medium === mk ? " selected" : "") + ">" + esc(_medGrp[mk].label) + "</option>";
    }).join("");
    var _medSel = '<select id="sim_medium" class="fip-select" title="高密表达培养基备选项" onchange="FIP.simMedium(this.value)">' + _medOpts + "</select>";
    var _mi = _medGrp[state.sim.medium] || {};
    var _medInfo = '<b>培养基：</b>' + esc(_mi.label || "") + '<br>' +
      '<b>组分：</b>' + esc(_mi.composition || "") + '<br>' +
      '<b>C:N：</b>' + esc(_mi.cn_ratio || "") + (_mi.note ? ('<br><b>说明：</b>' + esc(_mi.note)) : '');
    var sliderHtml = ranges.map(function (pr) {
      var val = (state.sim.loadedRecipe && state.sim.loadedRecipe[pr.key] != null) ? state.sim.loadedRecipe[pr.key]
        : ((state.sim.recipe && state.sim.recipe[pr.key] != null) ? state.sim.recipe[pr.key]
        : (pr.key in _defRec ? _defRec[pr.key] : pr.def));
      return '<div class="slider-row"><label>' + pr.label + (pr.unit ? " (" + pr.unit + ")" : "") +
        ' · <span class="sv" id="simv_' + pr.key + '">' + val + "</span></label>" +
        '<input type="range" id="sim_' + pr.key + '" min="' + pr.min + '" max="' + pr.max + '" step="' + pr.step + '" value="' + val + '" oninput="FIP.simSet(\'' + pr.key + "', this.value)\"></div>";
    }).join("");

    // 阶段 FBA 耦合面板（E. coli 专用）：碳源 / 菌株因子 / FBA 温度 / 表达水平 / 开关
    var _cs = state.sim.carbon_source || "glucose";
    var _sf = (state.sim.strain_factor != null) ? state.sim.strain_factor : 1.0;
    var _ft = (state.sim.fba_temp != null) ? state.sim.fba_temp : 37.0;
    var _el = (state.sim.expression_level != null) ? state.sim.expression_level : 0.15;
    // Pichia 表达策略面板用的状态量
    var _em = state.sim.expression_mode || "inducible";
    var _mm = state.sim.meoh_feed_mode || "methanol";
    var _mg = (state.sim.mixed_gly_feed_g_l != null) ? state.sim.mixed_gly_feed_g_l : 200.0;
    var _fbaCarbonName = (_em === "constitutive") ? "甘油（组成型全程）" : "甲醇（诱导相）";
    var fbaPanel = (state.sim.host === "ecoli")
      ? ('<div class="card" style="margin-top:12px;background:#10202a;border:1px solid #1f9e8933">' +
          '<div class="sub">⚛ 代谢耦合（FBA 两层接通 · 默认关闭）</div>' +
          '<div class="hint" style="margin-bottom:8px">勾选「启用 FBA 耦合」后，动态层继承 FBA 的 μ_het / 宿主摄取 / 乙酸 yield（浏览器侧由预计算查表驱动）。未勾选时用手调默认值。</div>' +
          '<div class="exp-grid">' +
            '<div class="exp-f"><label>碳源（FIP 二选一）</label><select id="sim_carbon_source" class="fip-select" style="min-width:120px" onchange="FIP.simCarbonSource(this.value)">' +
              '<option value="glucose"' + (_cs === "glucose" ? " selected" : "") + ">葡萄糖</option>" +
              '<option value="glycerol"' + (_cs === "glycerol" ? " selected" : "") + ">甘油</option></select></div>" +
            '<div class="exp-f"><label>菌株改造因子（野生型=1.0）</label><div style="display:flex;gap:6px;align-items:center"><input type="number" id="sim_strain_factor" min="0.5" max="2.0" step="0.05" value="' + _sf + '" style="width:74px" oninput="FIP.simStrain(this.value, \'num\')"><input type="range" id="sim_strain_factor_rng" min="0.5" max="2.0" step="0.05" value="' + _sf + '" style="flex:1;min-width:90px" oninput="FIP.simStrain(this.value, \'rng\')"></div><i id="simv_strain_factor">' + _sf + "</i></div>" +
            '<div class="exp-f"><label>FBA 代表性温度 ℃</label><div style="display:flex;gap:6px;align-items:center"><input type="number" id="sim_fba_temp" min="25" max="42" step="1" value="' + _ft + '" style="width:74px" oninput="FIP.simFbaTemp(this.value, \'num\')"><input type="range" id="sim_fba_temp_rng" min="25" max="42" step="1" value="' + _ft + '" style="flex:1;min-width:90px" oninput="FIP.simFbaTemp(this.value, \'rng\')"></div><i id="simv_fba_temp">' + _ft + "</i></div>" +
            '<div class="exp-f"><label>异源蛋白占 DCW</label><div style="display:flex;gap:6px;align-items:center"><input type="number" id="sim_expression_level" min="0.05" max="0.30" step="0.01" value="' + _el + '" style="width:74px" oninput="FIP.simExprLevel(this.value, \'num\')"><input type="range" id="sim_expression_level_rng" min="0.05" max="0.30" step="0.01" value="' + _el + '" style="flex:1;min-width:90px" oninput="FIP.simExprLevel(this.value, \'rng\')"></div><i id="simv_expression_level">' + _el + "</i></div>" +
          "</div>" +
          '<label style="display:flex;align-items:center;gap:6px;font-size:12px;color:#9fb3c8;margin-top:8px"><input type="checkbox" id="sim_fba_coupled" ' + (state.sim.fba_coupled ? "checked" : "") + ' onchange="FIP.simFbaToggle(this.checked)"> 启用 FBA 耦合（继承 μ_het）</label>' +
        "</div>")
      : "";

    // Pichia 表达策略 + FBA 耦合面板：表达模式 / 诱导型补料模式 / 混合补料甘油浓度 / 菌株因子 / FBA 温度 / 表达水平 / 开关
    var pichiaPanel = (state.sim.host === "pichia")
      ? ('<div class="card" style="margin-top:12px;background:#10202a;border:1px solid #b06ab333">' +
          '<div class="sub">🧬 Pichia 表达策略（组成型 / 诱导型）</div>' +
          '<div class="hint" style="margin-bottom:8px">组成型：碳源全程甘油、无甲醇诱导；诱导型：甘油批 → 甘油补料 → 碳饥饿 → 甲醇诱导（可选甲醇+甘油混合补料共利用）。</div>' +
          '<div class="exp-grid">' +
            '<div class="exp-f"><label>表达模式</label><select id="sim_expr_mode" class="fip-select" onchange="FIP.simExprMode(this.value)">' +
              '<option value="inducible"' + (_em === "inducible" ? " selected" : "") + ">诱导型（甲醇诱导）</option>" +
              '<option value="constitutive"' + (_em === "constitutive" ? " selected" : "") + ">组成型（全程甘油）</option></select></div>" +
            '<div class="exp-f"><label>诱导型补料模式</label><select id="sim_meoh_mode" class="fip-select" onchange="FIP.simMeohMode(this.value)">' +
              '<option value="methanol"' + (_mm === "methanol" ? " selected" : "") + ">纯甲醇</option>" +
              '<option value="mixed"' + (_mm === "mixed" ? " selected" : "") + ">甲醇 + 甘油混合补料</option></select></div>" +
            (_mm === "mixed"
              ? '<div class="exp-f"><label>混合补料甘油浓度</label><div style="display:flex;gap:6px;align-items:center"><input type="number" id="sim_mixed_gly" min="50" max="500" step="10" value="' + _mg + '" style="width:74px" oninput="FIP.simMixedGly(this.value)"><input type="range" id="sim_mixed_gly_rng" min="50" max="500" step="10" value="' + _mg + '" style="flex:1;min-width:90px" oninput="FIP.simMixedGly(this.value)"></div><i id="simv_mixed_gly">' + _mg + " g/L</i></div>"
              : "") +
          "</div>" +
          '<div class="sub" style="margin-top:10px">⚛ 代谢耦合（FBA 两层接通 · 默认关闭）</div>' +
          '<div class="hint" style="margin-bottom:8px">勾选「启用 FBA 耦合」后，动态层继承 iPP668 平行 biomass 的 μ_het 降幅（drop_frac，≡1/(1+f·FORM_FACTOR)，与温度/菌株无关）与甲醇 DO 受限溢出的甲醛 yield。FBA 校准碳源：<b>' + _fbaCarbonName + '</b>（组成型=甘油，诱导型=甲醇，按表达模式自动选择）。</div>' +
          '<div class="exp-grid">' +
            '<div class="exp-f"><label>菌株改造因子（野生型=1.0）</label><div style="display:flex;gap:6px;align-items:center"><input type="number" id="sim_strain_factor" min="0.5" max="2.0" step="0.05" value="' + _sf + '" style="width:74px" oninput="FIP.simStrain(this.value, \'num\')"><input type="range" id="sim_strain_factor_rng" min="0.5" max="2.0" step="0.05" value="' + _sf + '" style="flex:1;min-width:90px" oninput="FIP.simStrain(this.value, \'rng\')"></div><i id="simv_strain_factor">' + _sf + "</i></div>" +
            '<div class="exp-f"><label>FBA 代表性温度 ℃</label><div style="display:flex;gap:6px;align-items:center"><input type="number" id="sim_fba_temp" min="20" max="36" step="1" value="' + _ft + '" style="width:74px" oninput="FIP.simFbaTemp(this.value, \'num\')"><input type="range" id="sim_fba_temp_rng" min="20" max="36" step="1" value="' + _ft + '" style="flex:1;min-width:90px" oninput="FIP.simFbaTemp(this.value, \'rng\')"></div><i id="simv_fba_temp">' + _ft + "</i></div>" +
            '<div class="exp-f"><label>异源蛋白占 DCW</label><div style="display:flex;gap:6px;align-items:center"><input type="number" id="sim_expression_level" min="0.05" max="0.30" step="0.01" value="' + _el + '" style="width:74px" oninput="FIP.simExprLevel(this.value, \'num\')"><input type="range" id="sim_expression_level_rng" min="0.05" max="0.30" step="0.01" value="' + _el + '" style="flex:1;min-width:90px" oninput="FIP.simExprLevel(this.value, \'rng\')"></div><i id="simv_expression_level">' + _el + "</i></div>" +
          "</div>" +
          '<label style="display:flex;align-items:center;gap:6px;font-size:12px;color:#9fb3c8;margin-top:8px"><input type="checkbox" id="sim_fba_coupled" ' + (state.sim.fba_coupled ? "checked" : "") + ' onchange="FIP.simFbaToggle(this.checked)"> 启用 FBA 耦合（继承 μ_het）</label>' +
        "</div>")
      : "";

    body.innerHTML =
      '<div class="card">' +
      '<div class="sub">🎬 What-if 动态仿真 · 机理前向模型（与 Process / Scale-up Twin 共用 van\'t Riet kLa）</div>' +
      '<div style="display:flex;gap:10px;align-items:center;flex-wrap:wrap;margin:6px 0 12px;padding:8px 10px;background:#10202a;border:1px solid #1f9e8933;border-radius:8px">' +
        '<span class="ksub">📥 仿真参数来源</span>' +
        '<select id="sim_batch_src" class="fip-select" onchange="FIP.simLoadBatch(this.value)">' +
          '<option value="">— 选择历史批次 —</option>' +
          DATA.batches.map(function (b) { return '<option value="' + esc(b.batch_id) + '"' + (state.sim.srcBatch === b.batch_id ? " selected" : "") + '>' + esc(b.batch_id) + ' · ' + esc(b.host_organism || "") + ' · ' + fmt(b.working_volume_l, 0) + 'L</option>'; }).join("") +
        '</select>' +
        '<button class="btn2" onclick="FIP.simAdoptExp()">🧪 采用本次实验参数输入</button>' +
        '<button class="btn2" onclick="FIP.simClearSource()">✖ 清除来源</button>' +
        '<span id="sim_src_badge" class="ec-badge">' + simSourceBadge() + '</span>' +
      '</div>' +
      '<div style="display:flex;gap:16px;flex-wrap:wrap">' +
        '<div style="flex:1;min-width:300px">' +
          '<div style="display:flex;gap:8px;align-items:center;flex-wrap:wrap;margin-bottom:8px">' +
            '<select id="sim_host" class="fip-select" onchange="FIP.simHost(this.value)">' +
              '<option value="ecoli"' + (state.sim.host === "ecoli" ? " selected" : "") + ">大肠杆菌 E. coli (24h)</option>" +
              '<option value="pichia"' + (state.sim.host === "pichia" ? " selected" : "") + ">毕赤酵母 Pichia (96h)</option></select>" +
            '<span class="hint" style="margin:0 2px 0 4px">培养基</span>' + _medSel +
            '<span class="hint" style="margin:0 2px 0 4px">表达形式</span>' +
            '<select id="sim_form" class="fip-select" title="蛋白表达形式：决定产物定位与固有表达负担（阶段 H）" onchange="FIP.simForm(this.value)">' +
              '<option value="">— 自定义 —</option>' +
              Object.keys(FermentSim.EXPRESSION_FORMS).map(function (fk) {
                return '<option value="' + fk + '"' + (state.sim.form === fk ? " selected" : "") + ">" + esc(FermentSim.EXPRESSION_FORMS[fk].label) + "</option>";
              }).join("") +
            "</select>" +
            '<button class="btn" id="sim_playbtn" onclick="FIP.simPlay()">▶ 播放</button>' +
            '<button class="btn2" onclick="FIP.simPause()">⏸ 暂停</button>' +
            '<button class="btn2" onclick="FIP.simReset()">⟲ 重置</button>' +
            '<label class="hint" style="margin:0">完成于 <select id="sim_speed" class="fip-select" style="min-width:80px" onchange="FIP.simSpeed(parseFloat(this.value))">' +
              '<option value="3">3 分钟</option><option value="4" selected>4 分钟</option><option value="5">5 分钟</option><option value="1">1 分钟(快)</option><option value="0.0001">即时</option></select></label>' +
            '<label style="display:flex;align-items:center;gap:6px;font-size:12px;color:#9fb3c8"><input type="checkbox" id="sim_base" ' + (state.sim.baseline ? "checked" : "") + ' onchange="FIP.simToggleBase(this.checked)">对比基线</label>' +
          "</div>" +
          '<div id="sim_medium_info" class="method-note" style="margin:6px 0 2px">' + _medInfo + "</div>" +
        '<div id="sim_sliders" class="whatif-grid">' + sliderHtml + "</div>" +
        '<div id="sim_progress" class="hint" style="margin-top:8px">t = 0.0 h</div>' +
        fbaPanel +
        pichiaPanel +
        "</div>" +
        '<div style="width:330px"><canvas id="sim_bioreactor" width="300" height="356" style="background:#0a0f13;border:1px solid #1f9e8933;border-radius:10px"></canvas></div>' +
      "</div>" +
      '<div style="display:flex;gap:16px;flex-wrap:wrap;margin-top:12px">' +
        '<div style="flex:1;min-width:340px"><canvas id="sim_chart" width="620" height="280" style="background:#0a0f13;border:1px solid #1f9e8933;border-radius:10px;width:100%"></canvas></div>' +
        '<div id="sim_kpi" style="width:300px"></div>' +
      "</div>" +
      (state.sim.expRecipe ? '<div class="method-note" style="margin-top:6px">🧪 当前显示：实验设计三阶段方案（分段 pH/温度 + 各阶段碳源浓度）。拖动任意滑块或切换培养基即回到 What-if 滑块模式。' +
        (((state.sim.expRecipe._cn_advisory || []).length)
          ? '<br>⚠ 提示：' + state.sim.expRecipe._cn_advisory.map(esc).join("；") + '（软提示，未进入机理模型）'
          : '') + '</div>' : '') +
      '<div class="hint" id="sim_note" style="margin-top:8px">' + esc(D.note || "") + "</div>" +
      (state.sim.host === "ecoli"
        ? '<div class="card" style="margin-top:14px;background:#10202a;border:1px solid #1f9e8933">' +
          '<div id="sim_fba_compare"></div>' +
        "</div>"
        : (state.sim.host === "pichia"
          ? '<div class="card" style="margin-top:14px;background:#10202a;border:1px solid #b06ab333">' +
            '<div id="sim_fba_compare_pichia"></div>' +
          "</div>"
          : "")) +
      "</div>";

    // 绑定
    FIP.simSet = function (key, val) {
      val = parseFloat(val); var sv = $("simv_" + key); if (sv) sv.textContent = val;
      state.sim.expRecipe = null;
      if (state.sim.recipe) state.sim.recipe[key] = val;
      recompute(false);
    };
    // 阶段 H · 表达形式选择：把预设（ib_frac/burden_sol/burden_ib/tox_k）写入 4 个滑块并重算
    FIP.simForm = function (form) {
      state.sim.form = form || null;
      if (!form || !FermentSim.EXPRESSION_FORMS[form]) return;  // 自定义：保持滑块现状
      var fp = FermentSim.EXPRESSION_FORMS[form];
      [["ib_frac", fp.ib_frac], ["burden_sol", fp.burden_sol], ["burden_ib", fp.burden_ib], ["tox_k", fp.tox_k]].forEach(function (kv) {
        var el = $("sim_" + kv[0]); if (el) { el.value = kv[1]; var sv = $("simv_" + kv[0]); if (sv) sv.textContent = kv[1]; }
      });
      state.sim.expRecipe = null;
      recompute(false);
    };
    FIP.simPlay = function () { if (!state.sim.traj) recompute(true); if (state.sim.simTime >= state.sim.traj.t[state.sim.traj.n - 1]) state.sim.simTime = 0; startLoop(); };
    FIP.simPause = function () { state.sim.playing = false; updatePlayBtn(); };
    FIP.simReset = function () { state.sim.simTime = 0; state.sim.playing = false; state.sim.bubbles = []; if (state.sim.traj) drawAll(0); updatePlayBtn(); };
    FIP.simSpeed = function (m) { state.sim.speedMin = m; };
    // 阶段 FBA 耦合控件（E. coli 专用）
    FIP.simCarbonSource = function (v) { state.sim.carbon_source = v; state.sim.expRecipe = null; recompute(false); };
    FIP.simStrain = function (v, src) {
      v = parseFloat(v); if (isNaN(v)) return;
      var nf = $("sim_strain_factor"), rg = $("sim_strain_factor_rng"), sv = $("simv_strain_factor");
      if (src !== "num" && nf) nf.value = v;
      if (src !== "rng" && rg) rg.value = v;
      if (sv) sv.textContent = v;
      state.sim.strain_factor = v; state.sim.expRecipe = null; recompute(false);
    };
    FIP.simFbaTemp = function (v, src) {
      v = parseFloat(v); if (isNaN(v)) return;
      var nf = $("sim_fba_temp"), rg = $("sim_fba_temp_rng"), sv = $("simv_fba_temp");
      if (src !== "num" && nf) nf.value = v;
      if (src !== "rng" && rg) rg.value = v;
      if (sv) sv.textContent = v;
      state.sim.fba_temp = v; state.sim.expRecipe = null; recompute(false);
    };
    FIP.simExprLevel = function (v, src) {
      v = parseFloat(v); if (isNaN(v)) return;
      var nf = $("sim_expression_level"), rg = $("sim_expression_level_rng"), sv = $("simv_expression_level");
      if (src !== "num" && nf) nf.value = v;
      if (src !== "rng" && rg) rg.value = v;
      if (sv) sv.textContent = v;
      state.sim.expression_level = v; state.sim.expRecipe = null; recompute(false);
    };
    FIP.simFbaToggle = function (on) { state.sim.fba_coupled = !!on; state.sim.expRecipe = null; recompute(false); };
    // Pichia 表达策略控件：改变后重建滑块（混合补料滑块按补料模式显隐）并重算
    FIP.simExprMode = function (v) { state.sim.expression_mode = v; state.sim.expRecipe = null; renderSimAnim(body); };
    FIP.simMeohMode = function (v) { state.sim.meoh_feed_mode = v; state.sim.expRecipe = null; renderSimAnim(body); };
    FIP.simMixedGly = function (v) {
      v = parseFloat(v); if (isNaN(v)) return;
      var nf = $("sim_mixed_gly"), rg = $("sim_mixed_gly_rng"), sv = $("simv_mixed_gly");
      if (nf) nf.value = v; if (rg) rg.value = v; if (sv) sv.textContent = v + " g/L";
      state.sim.mixed_gly_feed_g_l = v; state.sim.expRecipe = null; recompute(false);
    };
    FIP.simToggleBase = function (on) { state.sim.baseline = on; recompute(false); };
    FIP.simMedium = function (mk) {
      state.sim.medium = mk; state.sim.simTime = 0; state.sim.playing = false; state.sim.bubbles = [];
      state.sim.recipe = null; state.sim.expRecipe = null; state.sim.traj = null;
      state.sim.loadedRecipe = null; state.sim.loadedFrom = null; state.sim.srcBatch = null;
      state.sim.form = null;
      renderSimAnim(body);  // 重建滑块（培养基默认）与轨迹
    };
    FIP.simHost = function (h) {
      state.sim.host = h; state.sim.medium = (FermentSim.DEFAULT_MEDIUM && FermentSim.DEFAULT_MEDIUM[h]) || "defined";
      state.sim.simTime = 0; state.sim.playing = false; state.sim.bubbles = [];
      state.sim.recipe = null; state.sim.expRecipe = null; state.sim.traj = null;
      state.sim.loadedRecipe = null; state.sim.loadedFrom = null; state.sim.srcBatch = null;
      state.sim.form = null;
      // 宿主专属 FBA 温度默认（E. coli 37℃ / Pichia 28℃），避免沿用他宿主值
      state.sim.fba_temp = (h === "pichia") ? 28.0 : 37.0;
      state.sim.expression_mode = (h === "pichia") ? (state.sim.expression_mode || "inducible") : state.sim.expression_mode;
      renderSimAnim(body);  // 重建滑块与轨迹
    };

    recompute(true);
    // 首屏：若 speed 选「即时」则直接画满；否则停在 t=0 等待播放
    if (state.sim.speedMin <= 0.001) { state.sim.simTime = state.sim.traj.t[state.sim.traj.n - 1]; drawAll(state.sim.traj.n - 1); }
    else drawAll(0);
  }

  // -------------------------------------------------------------------------
  // 对外接口
  // -------------------------------------------------------------------------
  var FIP = {
    nav: function (p) { state.page = p; renderShell(); },
    _state: function () { return state; },  // 调试/测试用：暴露内部状态
    navToggle: function (k) { if (k === "yeast") state.navFoldYeast = !state.navFoldYeast; else if (k === "et") state.navFoldET = !state.navFoldET; renderShell(); },
    set: function (path, val) {
      var parts = path.split("."); var o = state; for (var i = 0; i < parts.length - 1; i++) o = o[parts[i]]; o[parts[parts.length - 1]] = val;
      renderShell();
    },
    viewBatch: function (bid) { viewBatch(bid); },
    viewProcess: function (bid) { viewProcess(bid); },
    doeSetHost: function (host) {
      var box = $("doe_inner");
      if (box) box.innerHTML = doeHostBlock(host, true, null);
    },
    toggleTags: function () { toggleTags(); },
    toggleCmpBatch: function (bid) {
      var arr = state.process.cmp.batches, i = arr.indexOf(bid);
      if (i >= 0) arr.splice(i, 1); else arr.push(bid);
      renderCompare();
      var btn = $("cmp_batch_btn"); if (btn) btn.innerHTML = "已选 " + arr.length + " 批 ▾";
    },
    cmpParam: function (v) { state.process.cmp.param = v; renderCompare(); },
    cmpToggleDrop: function (e) {
      var d = $("cmp_drop");
      if (d) d.style.display = (d.style.display === "block") ? "none" : "block";
      if (e && e.stopPropagation) e.stopPropagation();
    },
    cmpCloseDrop: function () { var d = $("cmp_drop"); if (d) d.style.display = "none"; },
    renderCompare: function () { renderCompare(); },
    expSet: function (path, val) {
      var parts = path.split("."), o = state;
      for (var i = 0; i < parts.length - 1; i++) o = o[parts[i]];
      o[parts[parts.length - 1]] = val;
      renderExpForm();
    },
    expMediumFill: function (ph, mk) {
      var hk = hostKeyOf(state.process.exp.host_type || "ecoli");
      var grp = (FermentSim.MEDIA && FermentSim.MEDIA[hk]) || {};
      var m = grp[mk]; if (!m) return;
      var p = state.process.exp[ph];
      // 阶段感知碳源：BSM 在生长期为甘油、诱导期为甲醇；避免把「甘油→甲醇」整串盖到单阶段
      var cs = (m.c_source_by_phase && m.c_source_by_phase[ph]) ? m.c_source_by_phase[ph] : (m.c_source || p.c_source);
      p.c_source = cs;
      p.n_source = m.n_source || p.n_source;
      p.cn_ratio = m.cn_ratio || p.cn_ratio;
      renderExpForm();
    },
    expTogglePhase: function (k) { state.process.exp.phases[k] = !state.process.exp.phases[k]; renderExpForm(); },
    expAddStep: function (ph) {
      var arr = state.process.exp[ph].feed_steps; if (!arr) { arr = state.process.exp[ph].feed_steps = []; }
      var last = arr.length ? arr[arr.length - 1].t : 0;
      arr.push({ t: last, rate: 0 });
      renderExpForm();
    },
    expDelStep: function (ph, idx) {
      var arr = state.process.exp[ph].feed_steps; if (arr && arr.length > idx) arr.splice(idx, 1);
      renderExpForm();
    },
    expGenJSON: function () {
      var host = $("p_exp_json"); if (!host) return;
      var o = JSON.parse(JSON.stringify(state.process.exp)); o.generated_at = new Date().toISOString();
      host.style.display = "block"; host.textContent = JSON.stringify(o, null, 2);
    },
    expDownload: function () {
      var o = JSON.parse(JSON.stringify(state.process.exp)); o.generated_at = new Date().toISOString();
      var s = JSON.stringify(o, null, 2);
      try {
        var blob = new Blob([s], { type: "application/json" });
        var a = document.createElement("a"); a.href = URL.createObjectURL(blob);
        a.download = "exp_design_" + state.process.exp.vessel + "_" + state.process.exp.host_type + ".json";
        document.body.appendChild(a); a.click(); document.body.removeChild(a);
      } catch (ex) {}
    },
    ptAdoptExp: function () { ptAdoptExp(); },
    runExpSim: function () { runExpSim(); },
    simLoadBatch: function (bid) { simLoadBatch(bid); },
    simAdoptExp: function () { simAdoptExp(); },
    simClearSource: function () {
      state.sim.loadedRecipe = null; state.sim.loadedFrom = null; state.sim.srcBatch = null;
      if (state.sim._body) renderSimAnim(state.sim._body); else FIP.nav("simanim");
    },
    // 实验参数保存 / 载入（本机 localStorage）
    expSave: function () {
      try {
        fipSet("fip_exp_design_v1", JSON.stringify(state.process.exp));
        var m = $("p_exp_savemsg"); if (m) { m.textContent = "✓ 已保存参数" + (fipStorageMode() === "local" ? " 到本机" : "（当前预览环境仅会话内有效）") + " · " + new Date().toLocaleTimeString(); m.style.color = "#80ed99"; }
      } catch (ex) { var m2 = $("p_exp_savemsg"); if (m2) { m2.textContent = "保存失败: " + ex.message; m2.style.color = "#ff8c8c"; } }
    },
    expLoad: function () {
      try {
        var s = fipGet("fip_exp_design_v1");
        if (!s) { var m = $("p_exp_savemsg"); if (m) { m.textContent = "未找到已保存的参数"; m.style.color = "#ffce4d"; } return; }
        var o = JSON.parse(s);
        // 合并以防缺字段
        ["batch", "fedbatch", "induction"].forEach(function (ph) {
          if (o[ph]) state.process.exp[ph] = Object.assign({}, state.process.exp[ph], o[ph]);
        });
        state.process.exp = Object.assign({}, state.process.exp, o);
        var migrated = expNormalizeCConc(state.process.exp);
        renderExpForm();
        var m = $("p_exp_savemsg"); if (m) {
          m.textContent = "✓ 已载入已保存参数" + (migrated ? "（旧版碳源浓度已按「补料液浓度」语义纠正为 500/792 g/L）" : "");
          m.style.color = "#80ed99";
        }
      } catch (ex) { var m2 = $("p_exp_savemsg"); if (m2) { m2.textContent = "载入失败: " + ex.message; m2.style.color = "#ff8c8c"; } }
    },
    // 培养基配方弹窗（Excel 兼容 .xls 填写 / 导出）
    expOpenRecipe: function (ph) {
      ensureExpRecipeModal();
      window.__expRecipePhase = ph;
      var e = state.process.exp;
      if (!e[ph].recipe) e[ph].recipe = [];
      var label = { batch: "Batch 期", fedbatch: "Fed-batch 期", induction: "Induction 期" }[ph] || ph;
      var title = $("exp_recipe_title"); if (title) title.innerHTML = "📋 培养基配方 · " + label + "（C:N " + esc(e[ph].cn_ratio || "-") + "）";
      FIP.expRecipeRender();
      recipeLibRender();
      var m = $("exp_recipe_modal"); if (m) m.style.display = "flex";
    },
    expRecipeRender: function () {
      ensureExpRecipeModal();
      var ph = window.__expRecipePhase; if (!ph) return;
      var rec = state.process.exp[ph].recipe || [];
      var rows = rec.map(function (s, idx) {
        var cat = s.cat || "";
        var sel = '<select class="fip-select" style="padding:3px 5px;font-size:12px" onchange="FIP.expRecipeSet(' + idx + ', \'cat\', this.value)">' +
          '<option value="C"' + (cat === "C" ? " selected" : "") + '>碳源</option>' +
          '<option value="N"' + (cat === "N" ? " selected" : "") + '>氮源</option>' +
          '<option value="O"' + (cat === "O" ? " selected" : "") + '>其他</option>' +
          '</select>';
        return '<tr><td><input type="text" value="' + esc(s.comp || "") + '" onchange="FIP.expRecipeSet(' + idx + ', \'comp\', this.value)"></td>' +
          '<td><input type="number" step="any" value="' + esc(s.amount || "") + '" onchange="FIP.expRecipeSet(' + idx + ', \'amount\', this.value)"></td>' +
          '<td><input type="text" value="' + esc(s.unit || "") + '" onchange="FIP.expRecipeSet(' + idx + ', \'unit\', this.value)"></td>' +
          '<td>' + sel + '</td>' +
          '<td><button class="btn2" onclick="FIP.expRecipeDel(' + idx + ')">移除</button></td></tr>';
      }).join("");
      var tbl = $("exp_recipe_tbl");
      if (tbl) tbl.innerHTML = '<tr><th>组分 Component</th><th>用量</th><th>单位 / 备注</th><th>类别</th><th></th></tr>' + rows;
    },
    expRecipeAddRow: function () {
      var ph = window.__expRecipePhase; if (!ph) return;
      if (!state.process.exp[ph].recipe) state.process.exp[ph].recipe = [];
      state.process.exp[ph].recipe.push({ comp: "", amount: "", unit: "" });
      FIP.expRecipeRender();
    },
    expRecipeSet: function (idx, field, val) {
      var ph = window.__expRecipePhase; if (!ph) return;
      var rec = state.process.exp[ph].recipe || [];
      if (!rec[idx]) rec[idx] = {};
      rec[idx][field] = val;
    },
    expRecipeDel: function (idx) {
      var ph = window.__expRecipePhase; if (!ph) return;
      var rec = state.process.exp[ph].recipe || [];
      if (rec.length > idx) rec.splice(idx, 1);
      FIP.expRecipeRender();
    },
    expRecipeSave: function () {
      var ph = window.__expRecipePhase;
      if (ph && (state.process.exp[ph].recipe || []).some(function (r) { return r.cat === "C" || r.cat === "N"; })) FIP.expRecipeAutoCN();
      FIP.expRecipeClose(); renderExpForm();
    },
    // 配方库：保存当前阶段配方到本机，或从库载入到当前阶段
    expRecipeLibSave: function () {
      var nameEl = $("exp_recipe_libname"), name = nameEl ? nameEl.value.trim() : "";
      var ph = window.__expRecipePhase; if (!ph) return;
      if (!name) name = "未命名配方";
      var rows = (state.process.exp[ph].recipe || []).map(function (r) { return { comp: r.comp, amount: r.amount, unit: r.unit, cat: r.cat || "" }; });
      var lib = recipeLibGet().filter(function (x) { return x.name !== name; });
      lib.push({ name: name, phase: ph, rows: rows });
      try { fipSet("fip_recipe_lib_v1", JSON.stringify(lib)); } catch (ex) {}
      recipeLibRender();
      recipeLibRender2(); // 同步刷新实验设计卡片上的配方库下拉
      if (nameEl) nameEl.value = "";
    },
    expRecipeLibApply: function (name) {
      if (!name) return;
      var ph = window.__expRecipePhase; if (!ph) return;
      var lib = recipeLibGet(), x = null;
      for (var i = 0; i < lib.length; i++) { if (lib[i].name === name) { x = lib[i]; break; } }
      if (!x) return;
      state.process.exp[ph].recipe = (x.rows || []).map(function (r) { return { comp: r.comp, amount: r.amount, unit: r.unit, cat: r.cat || "" }; });
      FIP.expRecipeAutoCN(); // 还原类别后自动回填 碳源浓度 / C:N（含重渲染弹窗与表单）
    },
    // 从实验设计卡片直接调用已保存配方（无需打开弹窗）：保存到本阶段并自动回填 C:N
    expRecipeLibApply2: function (k, name) {
      if (!name || !k) return;
      var lib = recipeLibGet(), x = null;
      for (var i = 0; i < lib.length; i++) { if (lib[i].name === name) { x = lib[i]; break; } }
      if (!x) return;
      state.process.exp[k].recipe = (x.rows || []).map(function (r) { return { comp: r.comp, amount: r.amount, unit: r.unit, cat: r.cat || "" }; });
      window.__expRecipePhase = k;   // 供 expRecipeAutoCN 定位阶段（重算 C:N）
      FIP.expRecipeAutoCN();
    },
    expRecipeLibDel: function () {
      var sel = $("exp_recipe_libsel"), name = sel ? sel.value : "";
      if (!name) return;
      var lib = recipeLibGet().filter(function (x) { return x.name !== name; });
      try { fipSet("fip_recipe_lib_v1", JSON.stringify(lib)); } catch (ex) {}
      recipeLibRender();
      recipeLibRender2(); // 同步刷新实验设计卡片上的配方库下拉
    },
    expRecipeClose: function () { var m = $("exp_recipe_modal"); if (m) m.style.display = "none"; },
    // 配方上传：JSON（含 recipe 数组或裸数组）或 CSV（comp,amount,unit,cat）
    expRecipeUpload: function (file) {
      if (!file) return;
      var ph = window.__expRecipePhase; if (!ph) return;
      var reader = new FileReader();
      reader.onload = function () {
        try {
          var text = reader.result, rows = [];
          if (/\.json$/i.test(file.name) || /^\s*[[{]/.test(text)) {
            var obj = JSON.parse(text);
            var arr = obj.recipe || (Array.isArray(obj) ? obj : null);
            if (!arr) return;
            rows = arr.map(function (r) { return { comp: r.comp || "", amount: (r.amount != null ? r.amount : ""), unit: r.unit || "", cat: r.cat || "" }; });
          } else {
            rows = text.split(/\r?\n/).filter(function (l) { return l.trim(); }).map(function (l) {
              var p = l.split(/[,;\t]/); return { comp: (p[0] || "").trim(), amount: (p[1] || "").trim(), unit: (p[2] || "").trim(), cat: (p[3] || "").trim() };
            });
          }
          state.process.exp[ph].recipe = rows;
          FIP.expRecipeRender();
        } catch (ex) { alert("上传失败：" + ex.message); }
      };
      reader.readAsText(file);
    },
    // 按配方自动计算 碳源浓度 / C:N 比，回填到阶段字段
    expRecipeAutoCN: function () {
      var ph = window.__expRecipePhase; if (!ph) return;
      expRecipeGuessCat(ph);
      var rec = state.process.exp[ph].recipe || [], cMass = 0, nMass = 0, cConc = 0, cSource = "", nSource = "";
      rec.forEach(function (r) {
        var amt = parseFloat(r.amount); if (isNaN(amt)) return;
        var unit = (r.unit || "").toLowerCase();
        var a = /%|v\/v|w\/v/.test(unit) && !/g\/l|g$/.test(unit) ? amt * 10 : amt; // 粗估 % -> g/L
        var fr = cnFractions(r.comp);
        if (r.cat === "C") { cConc += a; cMass += a * fr.C; if (!cSource) cSource = r.comp; }
        else if (r.cat === "N") { nMass += a * fr.N; if (!nSource) nSource = r.comp; }
      });
      var o = state.process.exp[ph];
      if (cConc > 0) { o.c_conc = Math.round(cConc * 10) / 10; if (cSource) o.c_source = cSource; }
      if (nSource) o.n_source = nSource;
      o.cn_ratio = nMass > 0 ? (fmt(cMass / nMass, 1) + ":1") : (cConc > 0 ? "—" : o.cn_ratio);
      FIP.expRecipeRender(); renderExpForm();
    },
    // 补料策略库：保存 / 选择应用 / 删除 / 上传
    feedStratSave: function (ph) {
      var nameEl = $("feed_strat_name_" + ph), name = nameEl ? nameEl.value.trim() : "";
      if (!name) name = "未命名策略";
      var p = state.process.exp[ph];
      var strat = { feed_mode: p.feed_mode, feed_const: p.feed_const, feed_steps: (p.feed_steps || []).map(function (s) { return { t: s.t, rate: s.rate }; }) };
      var lib = feedStratLibGet().filter(function (x) { return x.name !== name; });
      lib.push({ name: name, phase: ph, strat: strat });
      try { fipSet("fip_feed_lib_v1", JSON.stringify(lib)); } catch (ex) {}
      feedStratLibRender(); if (nameEl) nameEl.value = "";
    },
    feedStratApply: function (ph, name) {
      if (!name) return;
      var lib = feedStratLibGet(), x = null;
      for (var i = 0; i < lib.length; i++) { if (lib[i].name === name) { x = lib[i]; break; } }
      if (!x) return; applyFeedStrat(ph, x.strat); renderExpForm();
    },
    feedStratDel: function (ph) {
      var sel = $("feed_strat_sel_" + ph), name = sel ? sel.value : "";
      if (!name) return;
      var lib = feedStratLibGet().filter(function (x) { return x.name !== name; });
      try { fipSet("fip_feed_lib_v1", JSON.stringify(lib)); } catch (ex) {}
      feedStratLibRender();
    },
    feedStratUpload: function (ph, file) {
      if (!file) return;
      var reader = new FileReader();
      reader.onload = function () {
        try {
          var obj = JSON.parse(reader.result), applied = false;
          if (obj.batch || obj.fedbatch || obj.induction) {
            ["batch", "fedbatch", "induction"].forEach(function (k) { if (obj[k] && obj[k].feed_mode) { applyFeedStrat(k, obj[k]); applied = true; } });
          } else if (obj.feed_mode) { applyFeedStrat(ph, obj); applied = true; }
          else if (Array.isArray(obj)) { applyFeedStrat(ph, { feed_mode: "stepwise", feed_steps: obj.map(function (s) { return { t: s.t, rate: s.rate }; }) }); applied = true; }
          if (applied) renderExpForm();
        } catch (ex) { alert("上传失败：" + ex.message); }
      };
      reader.readAsText(file);
    },
    expRecipeExport: function () {
      var ph = window.__expRecipePhase; if (!ph) return;
      var e = state.process.exp;
      var label = { batch: "Batch", fedbatch: "Fedbatch", induction: "Induction" }[ph] || ph;
      var cells = [
        ["阶段 Phase", label], ["碳源 C-source", e[ph].c_source || ""], ["碳源浓度(g/L)", e[ph].c_conc != null ? e[ph].c_conc : ""], ["氮源 N-source", e[ph].n_source || ""],
        ["C:N 比例", e[ph].cn_ratio || ""], ["发酵罐", e.vessel], ["起始装液量(kg)", e.fill_kg]
      ];
      var xml = '<?xml version="1.0"?>\n<?mso-application progid="Excel.Sheet"?>\n<Workbook xmlns="urn:schemas-microsoft-com:office:spreadsheet" xmlns:ss="urn:schemas-microsoft-com:office:spreadsheet"><Worksheet ss:Name="Recipe"><Table>';
      cells.forEach(function (c) { xml += expRecipeRowXml(c); });
      xml += expRecipeRowXml(["—— 培养基配方 ——", "", "", ""]);
      xml += expRecipeRowXml(["组分 Component", "用量", "单位/备注", "类别"]);
      (e[ph].recipe || []).forEach(function (s) {
        var catTxt = s.cat === "C" ? "碳源" : s.cat === "N" ? "氮源" : s.cat === "O" ? "其他" : "";
        xml += expRecipeRowXml([s.comp || "", s.amount || "", s.unit || "", catTxt]);
      });
      xml += "</Table></Worksheet></Workbook>";
      ioDownload("培养基配方_" + label + ".xls", xml, "application/vnd.ms-excel");
    },
    calcExpr: function () { calcExpr(); },
    calcScale: function () { calcScale(); },
    scaleSetFermenter: function (path, val) {
      var parts = (path || "").split(".");
      if (parts.length < 2) return;
      var tag = parts[0], field = parts[1];
      var f = state.scale.fermenters[tag];
      if (!f) return;
      if (field === "name") { f.name = val; }
      else {
        var n = parseFloat(val);
        if (isNaN(n)) return;
        if (field === "n_imp") {
          n = Math.max(1, Math.min(8, Math.round(n)));
          var arr = (f.imps && f.imps.length) ? f.imps.slice() : [f.imp || "rushton"];
          while (arr.length < n) arr.push("rushton");
          f.imps = arr.slice(0, n);
          f.n_imp = n;
        } else {
          f[field] = n;
        }
      }
      renderFermenters($("s_fermenters"));
      if ($("s_fchart")) renderFermenterChart($("s_fchart"));
      if ($("s_sweep")) renderSweep($("s_sweep"));
      if ($("s_risk")) renderScaleRisk($("s_risk"));
      calcScale();
    },
    scaleSetImp: function (tag, idx, val) {
      var f = state.scale.fermenters[tag];
      if (!f) return;
      if (!f.imps) f.imps = [f.imp || "rushton"];
      f.imps[idx] = val;
      renderFermenters($("s_fermenters"));
      if ($("s_fchart")) renderFermenterChart($("s_fchart"));
      if ($("s_sweep")) renderSweep($("s_sweep"));
      calcScale();
    },
    scaleAddFermenter: function () {
      var fers = state.scale.fermenters;
      var n = 2; while (fers["X" + n]) n++;
      var tag = "X" + n;
      // 以最大尺度为基准放大生成（体积 ×8、几何按比例）
      var tags = Object.keys(fers);
      var base = fers[tags[tags.length - 1]];
      var f = Math.pow(8, 1 / 3);
      var V_L = Math.round(base.V_L * 8);
      var nf = genFermenter(V_L, Math.round(base.N / Math.pow(f, 2 / 3)));
      nf.tag = tag; nf.name = "自定义 " + V_L + " L";
      fers[tag] = nf;
      renderFermenters($("s_fermenters"));
      if ($("s_fchart")) renderFermenterChart($("s_fchart"));
      calcScale();
    },
    scaleDelFermenter: function (tag) {
      var fers = state.scale.fermenters;
      if (Object.keys(fers).length <= 1) return;
      delete fers[tag];
      if (state.scale.src === tag) state.scale.src = Object.keys(fers)[0];
      if (state.scale.tgt === tag) state.scale.tgt = Object.keys(fers)[Object.keys(fers).length > 1 ? 1 : 0];
      var inc = state.scale.chart.include;
      if (inc) state.scale.chart.include = inc.filter(function (t) { return t !== tag; });
      renderFermenters($("s_fermenters"));
      if ($("s_fchart")) renderFermenterChart($("s_fchart"));
      calcScale();
    },
    scaleSetMu: function (v) {
      var n = parseFloat(v); if (isNaN(n) || n <= 0) return;
      state.scale.mu = n;
      renderFermenters($("s_fermenters"));
      if ($("s_fchart")) renderFermenterChart($("s_fchart"));
      if ($("s_sweep")) renderSweep($("s_sweep"));
      calcScale();
    },
    scaleSetProc: function (key, v) {
      var n = parseFloat(v);
      if (isNaN(n)) return;
      state.scale.proc[key] = n;
      if ($("s_risk")) renderScaleRisk($("s_risk"));
      calcScale();
    },
    scaleSetCorr: function (key, v) {
      var n = parseFloat(v);
      if (isNaN(n) || n <= 0) return;
      if (!state.scale.proc.corr) state.scale.proc.corr = { kla: 1, tmix: 1, heat: 1 };
      state.scale.proc.corr[key] = n;
      if ($("s_risk")) renderScaleRisk($("s_risk"));
      calcScale();
    },
    scaleSetHost: function (host) { scaleSetHost(host); },
    scaleLinkProcess: function () { scaleLinkFromProcess(); },
    scaleUnlinkProcess: function () { scaleUnlinkProcess(); },
    scaleSetSweepTag: function (tag) {
      state.scale.sweep.tag = tag;
      if ($("s_sweep")) renderSweep($("s_sweep"));
    },
    scaleSweepCmpToggle: function (tag) {
      var sw = state.scale.sweep;
      var all = Object.keys(state.scale.fermenters).sort(function (a, b) { return state.scale.fermenters[a].V_L - state.scale.fermenters[b].V_L; });
      if (!sw.cmp) sw.cmp = all.slice();
      var i = sw.cmp.indexOf(tag);
      if (i >= 0) { if (sw.cmp.length > 1) sw.cmp = sw.cmp.filter(function (t) { return t !== tag; }); }
      else sw.cmp = sw.cmp.concat([tag]);
      if ($("s_sweep")) renderSweep($("s_sweep"));
    },
    scaleSweepCmpAll: function () {
      state.scale.sweep.cmp = Object.keys(state.scale.fermenters).sort(function (a, b) { return state.scale.fermenters[a].V_L - state.scale.fermenters[b].V_L; });
      if ($("s_sweep")) renderSweep($("s_sweep"));
    },
    scaleSweepCmpMetric: function (k) {
      state.scale.sweep.cmpMetric = k;
      if ($("s_sweep")) renderSweep($("s_sweep"));
    },
    scaleSaveFermenters: function () {
      var ok = scaleSaveFermenters();
      if (!ok) { alert("保存失败：当前环境不支持本地存储（预览沙箱）。请在浏览器中打开本页面后再保存。"); }
      renderFermenters($("s_fermenters"));
      if ($("s_fchart")) renderFermenterChart($("s_fchart"));
      if ($("s_sweep")) renderSweep($("s_sweep"));
      calcScale();
    },
    scaleResetFermenters: function () {
      if (!confirm("确认恢复为 10 档默认发酵罐结构参数？已保存的自定义参数将被清除。")) return;
      scaleResetFermenters();
      renderFermenters($("s_fermenters"));
      if ($("s_fchart")) renderFermenterChart($("s_fchart"));
      if ($("s_sweep")) renderSweep($("s_sweep"));
      calcScale();
    },
    scaleExportFermenters: function () {
      var ok = scaleExportFermenters();
      if (!ok) alert("导出失败：当前环境不支持文件下载。请改用「保存结构参数」(浏览器本地存储)，或在浏览器中直接打开本页面。");
    },
    scaleImportFermenters: function (file) {
      scaleImportFermenters(file);
    },
    scaleSetHcalc: function (field, val) {
      var hc = state.scale.hcalc || (state.scale.hcalc = { D: 0.12, V: 2 });
      var v = parseFloat(val);
      if (!isNaN(v) && v >= 0) hc[field] = v;
      if ($("s_fermenters")) renderFermenters($("s_fermenters"));
    },
    scaleApplyHcalc: function (tag) {
      var f = state.scale.fermenters[tag];
      if (!f) return;
      var hc = state.scale.hcalc || (state.scale.hcalc = { D: 0.12, V: 2 });
      var hcD = +hc.D, hcV = +hc.V;
      if (!(hcD > 0)) return;
      f.H = +(4 * (hcV / 1000) / (Math.PI * hcD * hcD)).toFixed(4);
      if ($("s_fermenters")) renderFermenters($("s_fermenters"));
      if ($("s_fchart")) renderFermenterChart($("s_fchart"));
      if ($("s_sweep")) renderSweep($("s_sweep"));
      calcScale();
    },
    scaleSetChartMetric: function (k) {
      state.scale.chart.metric = k;
      if ($("s_fchart")) renderFermenterChart($("s_fchart"));
    },
    scaleToggleFermenter: function (tag) {
      var ch = state.scale.chart, inc = ch.include || Object.keys(state.scale.fermenters);
      var i = inc.indexOf(tag);
      if (i >= 0) { if (inc.length > 1) inc = inc.filter(function (t) { return t !== tag; }); }
      else inc = inc.concat([tag]);
      ch.include = inc;
      if ($("s_fchart")) renderFermenterChart($("s_fchart"));
    },
    scaleChartAll: function () {
      state.scale.chart.include = null;
      if ($("s_fchart")) renderFermenterChart($("s_fchart"));
    },
    example: function (k) { exampleFill(k); },
    saveFeedback: function () { saveFeedback(); },
    clearFeedback: function () { clearFeedback(); },
    ask: function () {
      var q = $("c_q").value; state.copilot = q;
      $("c_a").textContent = copilotAnswer(q);
    },
    askWith: function (i) {
      var q = COPILOT_EXAMPLES[i]; if (!q) return;
      var ta = $("c_q"); if (ta) ta.value = q;
      state.copilot = q;
      $("c_a").textContent = copilotAnswer(q);
    },
    // ---- Expression Twin / E.coli 模块交互接口 ----
    ecoliAnalyze: function () { ecoliRecompute(); },
    ecoliExample: function () { var ta = $("e2_seq"); if (ta) { ta.value = ECOLI_SAMPLE; state.ecoli.seq = ECOLI_SAMPLE; } ecoliRecompute(); },
    ecoliClear: function () { var ta = $("e2_seq"); if (ta) { ta.value = ""; state.ecoli.seq = ""; } ecoliRecompute(); },
    ecoliSetSp: function (v) { state.ecoli.sp = v; ecoliRecompute(); },
    ecoliSetModel: function (v) { state.ecoli.model = v; ecoliRecompute(); },
    ecoliSaveFb: function () { ecoliSaveFb(); },
    ecoliTrainFb: function () { ecoliTrainFb(); },
    ecoliClearFb: function () { ecoliClearFb(); },
    // ---- Expression Twin / 酵母 L1-L7 层级详情页接口 ----
    yeastRecompute: function () { if (state.page && state.page.indexOf("yL") === 0) yeastLevelRecompute(parseInt(state.page.slice(2), 10) - 1); },
    // ---- 批次数据 导入/导出接口 ----
    exportCSV: function () { ioDownload("fip_batches.csv", ioBatchesToCSV(), "text/csv;charset=utf-8"); var m = $("io_msg"); if (m) m.textContent = "已导出 " + DATA.batches.length + " 批 (CSV)"; },
    exportJSON: function () {
      var payload = { batches: DATA.batches, timeseries: DATA.timeseries, exported_at: new Date().toISOString() };
      ioDownload("fip_data.json", JSON.stringify(payload, null, 2), "application/json");
      var m = $("io_msg"); if (m) m.textContent = "已导出 " + DATA.batches.length + " 批 + 时序 (JSON)";
    },
    downloadTemplate: function () {
      var sample = [BATCH_COLS.join(","), "B2026-099,E.coli_BL21(DE3),S,2,24,12.5,0.85,true,"];
      ioDownload("fip_batch_template.csv", "﻿" + sample.join("\n"), "text/csv;charset=utf-8");
    },
    // 在线位号时序模板：列名即 Tier-1/2 位号，供真实 SCADA/DCS 导出直接对齐。
    // 导入时建议按 1~5 min 采样；缺列不会报错，但对应模块会按下方「为什么需要」降级。
    downloadTsTemplate: function () {
      var T = DATA.tag_tiers || {};
      var cols = ["batch_id", "hours"];
      ["tier1_required", "tier2_recommended", "tier3_optional"].forEach(function (k) {
        Object.keys(T[k] || {}).forEach(function (t) { cols.push(t); });
      });
      var cn = ["批次号", "时间(h)"];
      ["tier1_required", "tier2_recommended", "tier3_optional"].forEach(function (k) {
        var spec = T[k] || {};
        Object.keys(spec).forEach(function (t) { cn.push(spec[t].cn + "(" + spec[t].unit + ")"); });
      });
      var row = ["B2026-099", "0"];
      for (var i = 2; i < cols.length; i++) row.push("");
      var sample = [cols.join(","), cn.join(","), row.join(",")];
      ioDownload("fip_timeseries_template.csv", "﻿" + sample.join("\n"), "text/csv;charset=utf-8");
    },
    importFile: function (file) {
      if (!file) return;
      var name = (file.name || "").toLowerCase(), m = $("io_msg");
      var reader = new FileReader();
      reader.onload = function (ev) {
        try {
          var text = ev.target.result;
          var msg = name.indexOf(".json") >= 0 ? FIP._importJSON(text) : FIP._importCSV(text);
          if (m) m.textContent = msg;
          renderBatches($("body-content"));
        } catch (e) { if (m) m.textContent = "导入失败: " + e.message; }
      };
      reader.onerror = function () { if (m) m.textContent = "文件读取失败"; };
      reader.readAsText(file);
    },
    _importCSV: function (text) {
      var rows = ioParseCSV(text);
      if (!rows.length) throw new Error("CSV 无数据行");
      return FIP._upsert(rows, null);
    },
    _importJSON: function (text) {
      var obj = JSON.parse(text), batches = null, ts = null;
      if (Array.isArray(obj)) batches = obj;
      else { batches = obj.batches || null; ts = obj.timeseries || obj.time_series || null; }
      if (!batches) throw new Error("JSON 缺少 batches 字段");
      return FIP._upsert(batches, ts);
    },
    _upsert: function (batches, ts) {
      var n = 0;
      batches.forEach(function (b) {
        if (b.batch_id === undefined || b.batch_id === null || b.batch_id === "") return;
        ioCoerceBatch(b);
        var i = -1;
        for (var k = 0; k < DATA.batches.length; k++) { if (String(DATA.batches[k].batch_id) === String(b.batch_id)) { i = k; break; } }
        if (i >= 0) DATA.batches[i] = Object.assign({}, DATA.batches[i], b); else DATA.batches.push(b);
        n++;
      });
      var nts = 0, nsoft = 0;
      if (ts) {
        for (var key in ts) {
          if (!ts.hasOwnProperty(key)) continue;
          DATA.timeseries[key] = ts[key]; nts++;
          // 真实数据通常不含软测量列（DCW/WCW 是软测量产物），导入后即时用
          // 与后端同构的氧衡算 ODE 补算，保证 Process Twin / Scale-up 关联可用。
          if ((!ts[key].DCW || !ts[key].WCW) && ts[key].OUR) { if (ptRebuildSoftSensor(key, ts[key])) nsoft++; }
        }
      }
      return "已导入 " + n + " 批" + (ts ? "（含时序 " + nts + " 条" + (nsoft ? "，其中 " + nsoft + " 条已重建软测量 DCW/WCW" : "") + "）" : "");
    },
    _batchesCSV: function () { return ioBatchesToCSV(); },

    // ---- Process Twin 模型 API（纯函数，供自动化测试/外部复用）----
    pt: {
      softMeta: ptSoftMeta,
      dcwFromGas: ptDcwFromGas,
      trapz: ptTrapz,
      scenarioFactors: ptScenarioFactors,
      whatIfTraj: ptWhatIfTraj,
      xmax: ptXmax
    }
  };
  window.FIP = FIP;

  // -------------------------------------------------------------------------
  // 启动
  // -------------------------------------------------------------------------
  document.addEventListener("DOMContentLoaded", function () {
    DATA = window.DASH_DATA || (typeof DASH_DATA !== "undefined" ? DASH_DATA : null);
    // 自动加载已保存的发酵罐结构参数（若用户此前点击过“保存结构参数”）
    var _sf = scaleLoadFermenters();
    if (_sf) state.scale.fermenters = _sf;
    expNormalizeCConc(state.process.exp);   // 碳源浓度语义迁移（发酵液浓度 → 补料液浓度）
    renderShell();
    // 点击页面其它区域时收起批次多选下拉面板
    document.addEventListener("click", function (e) {
      try {
        var w = e.target.closest && e.target.closest(".cmp-drop-wrap");
        if (!w) { var d = document.getElementById("cmp_drop"); if (d) d.style.display = "none"; }
      } catch (ex) {}
    });
  });
})();
