// Sanitising coordinates and red flags (port of R/sanitize_stream.R).

import * as rs from "./rstats.js";
import { asTable, attrs, checkColumns, diffBySegment, filterRows, isFirstOfSegment, numCol, segmentId,
  toLonlat, toMetric, uniqueInOrder, utmEpsg } from "./utils.js";
import { tobler_speed } from "./tobler.js";

const METHODS = ["kalman", "loess", "gaussian", "savgol", "pipeline", "none"];
const ACTIONS = ["flag", "drop", "interpolate"];
const FLAG_NAMES = ["flag_time", "flag_speed", "flag_accel", "flag_jump", "flag_slope", "flag_vert_speed",
  "flag_idle_ele", "flag_gap", "flag_low_res"];

const rowsOf = (seg, s) => {
  const r = [];
  seg.forEach((v, i) => { if (v === s) r.push(i); });
  return r;
};

export function sanitize_stream(data, {
  method = "kalman", action = "flag", smooth_ele = false,
  v_max = 10, accel_max = 6, vert_speed_max = 1, slope_max = 1.19,
  gap = 30, low_res = 20, edge_jump = 30, idle = 0.1,
  tobler = null, tol_factor = 1.1,
  sigma_obs = 5, sigma_acc = 1, rts = true,
  alpha = 0.05, bandwidth = 10,
  sg_window = 9, sg_order = 2,
  window = 25, jump_k = 1.5, local_k = 2,
} = {}) {
  if (!METHODS.includes(method)) throw new Error(`'method' must be one of ${METHODS.join(", ")}`);
  if (!ACTIONS.includes(action)) throw new Error(`'action' must be one of ${ACTIONS.join(", ")}`);
  checkColumns(data, ["lon", "lat"], "sanitize_stream");
  let t = asTable(data);
  let hasTime = "time" in t;
  const hasEle = "ele" in t;

  // ---- 1. flags on raw values -------------------------------------------------
  const lon0 = numCol(t.lon);
  const lat0 = numCol(t.lat);
  const flagCoord = lon0.map((lo, i) => {
    const la = lat0[i];
    return Number.isNaN(lo) || Number.isNaN(la) || Math.abs(la) > 90 || Math.abs(lo) > 180 || (lo === 0 && la === 0);
  });
  const nCoord = flagCoord.filter(Boolean).length;
  t = filterRows(t, flagCoord.map((f) => !f));
  const n = t.lon.length;
  if (n < 3) throw new Error("Fewer than 3 readings with valid coordinates.");

  const lon = numCol(t.lon);
  const lat = numCol(t.lat);
  const epsg = utmEpsg(lon, lat);
  let [x, y] = toMetric(lon, lat, epsg);
  let seg = segmentId(t);
  const first = isFirstOfSegment(seg);

  const ddx = diffBySegment(x, seg);
  const ddy = diffBySegment(y, seg);
  const deltaD = ddx.map((v, i) => Math.sqrt(v * v + ddy[i] * ddy[i]));
  const F = Object.fromEntries(FLAG_NAMES.map((f) => [f, new Array(n).fill(false)]));
  F.flag_low_res = deltaD.map((d) => d > low_res);
  const speed = new Array(n).fill(NaN);

  let ele = hasEle ? numCol(t.ele) : null;
  const deltaE = hasEle ? diffBySegment(ele, seg) : new Array(n).fill(0);
  if (hasEle) {
    const moving = deltaD.map((d) => d >= idle);
    F.flag_slope = deltaE.map((e, i) => moving[i] && Math.abs(e / deltaD[i]) > slope_max);
    F.flag_idle_ele = deltaE.map((e, i) => !moving[i] && !first[i] && Math.abs(e) > 0);
  }
  if (hasTime) {
    const tt = numCol(t.time);
    const dt = diffBySegment(tt, seg);
    F.flag_time = dt.map((d, i) => !first[i] && (Number.isNaN(d) || d <= 0));
    F.flag_gap = dt.map((d, i) => !first[i] && !Number.isNaN(d) && d > gap);
    const okT = dt.map((d, i) => !first[i] && !Number.isNaN(d) && d > 0);
    okT.forEach((ok, i) => { if (ok) speed[i] = deltaD[i] / dt[i]; });
    let speedMax;
    if (tobler === null) speedMax = new Array(n).fill(v_max);
    else {
      const slope = new Array(n).fill(0);
      if (hasEle) okT.forEach((ok, i) => { if (ok) slope[i] = deltaD[i] > idle ? deltaE[i] / deltaD[i] : 0; });
      speedMax = tobler_speed(slope.map((s) => Math.min(Math.max(s, -0.6), 0.6)), tobler).map((v) => v * tol_factor);
    }
    F.flag_speed = speed.map((s, i) => !Number.isNaN(s) && s > speedMax[i]);
    const ds = diffBySegment(speed, seg, NaN);
    const accel = ds.map((v, i) => (okT[i] ? v / dt[i] : NaN));
    F.flag_accel = accel.map((a) => !Number.isNaN(a) && Math.abs(a) > accel_max);
    if (hasEle) {
      F.flag_vert_speed = deltaE.map((e, i) => okT[i] && Math.abs(e / dt[i]) > vert_speed_max);
    }
    for (const s of uniqueInOrder(seg)) {
      const rows = rowsOf(seg, s);
      const k = rows.length;
      if (k < 7) continue;
      for (let i = 1; i <= 3; i++) {
        if (deltaD[rows[i]] > edge_jump) F.flag_jump[rows[i - 1]] = true;
        if (deltaD[rows[k - i]] > edge_jump) F.flag_jump[rows[k - i]] = true;
      }
    }
  }
  if (hasTime) {
    const sp = speed.filter((s, i) => !first[i] && !Number.isNaN(s));
    if (sp.length && sp.filter((s) => s > v_max).length / sp.length > 0.5) {
      console.warn("More than half of the readings exceed the speed limit: time is not physical and the stream is smoothed as a route without time.");
      hasTime = false;
      for (const f of ["flag_speed", "flag_accel", "flag_vert_speed", "flag_gap", "flag_jump"]) F[f].fill(false);
    }
  }
  const positionFlag = F.flag_time.map((_, i) => F.flag_time[i] || F.flag_speed[i] || F.flag_accel[i] || F.flag_jump[i]);
  const eleFlag = F.flag_slope.map((_, i) => F.flag_slope[i] || F.flag_vert_speed[i] || F.flag_idle_ele[i]);

  const counts = {};
  for (const f of FLAG_NAMES) counts[f] = F[f].filter(Boolean).length;
  const idleE = deltaE.filter((_, i) => F.flag_idle_ele[i]);
  const report = {
    n_readings: n + nCoord,
    flag_coord: nCoord,
    counts,
    share_integer_ele: hasEle ? rs.rMean(rs.naOmit(ele).map((e) => (e === Math.round(e) ? 1 : 0))) : NaN,
    share_low_res: rs.rMean(F.flag_low_res.map(Number)),
    idle_gain: hasEle ? rs.rSum(idleE.map((e) => Math.max(e, 0))) : NaN,
    idle_loss: hasEle ? rs.rSum(idleE.map((e) => Math.min(e, 0))) : NaN,
  };
  report.flag_ele_rounding = !Number.isNaN(report.share_integer_ele) && report.share_integer_ele > 0.9;

  // ---- 2. act on flagged readings --------------------------------------------------
  if (action === "drop") {
    if (hasEle) {
      ele = ele.map((e, i) => (eleFlag[i] ? NaN : e));
      t.ele = ele;
    }
    const keep = positionFlag.map((f) => !f);
    t = filterRows(t, keep);
    for (const f of FLAG_NAMES) F[f] = F[f].filter((_, i) => keep[i]);
    x = x.filter((_, i) => keep[i]);
    y = y.filter((_, i) => keep[i]);
    seg = seg.filter((_, i) => keep[i]);
  } else if (action === "interpolate") {
    x = x.map((v, i) => (positionFlag[i] ? NaN : v));
    y = y.map((v, i) => (positionFlag[i] ? NaN : v));
    if (hasEle) ele = ele.map((e, i) => (eleFlag[i] || positionFlag[i] ? NaN : e));
    const tt = hasTime ? numCol(t.time) : null;
    for (const s of uniqueInOrder(seg)) {
      const rows = rowsOf(seg, s);
      const base = hasTime ? rows.map((i) => tt[i]) : rows.map((_, k) => k + 1);
      const fx = fillLinear(base, rows.map((i) => x[i]));
      const fy = fillLinear(base, rows.map((i) => y[i]));
      rows.forEach((r, k) => { x[r] = fx[k]; y[r] = fy[k]; });
      if (hasEle) {
        const fe = fillLinear(base, rows.map((i) => ele[i]));
        rows.forEach((r, k) => { ele[r] = fe[k]; });
      }
    }
    if (hasEle) t.ele = ele;
  }

  // ---- 3. smoothing ------------------------------------------------------------------
  if (method !== "none") {
    const tt = hasTime ? numCol(t.time) : null;
    const eleNow = hasEle ? numCol(t.ele) : null;
    for (const s of uniqueInOrder(seg)) {
      const rows = rowsOf(seg, s);
      if (rows.length < 5) continue;
      const xs = rows.map((i) => x[i]);
      const ys = rows.map((i) => y[i]);
      let tSeg;
      if (hasTime) tSeg = rows.map((i) => tt[i]);
      else {
        tSeg = [0];
        for (let k = 1; k < xs.length; k++) tSeg.push(tSeg[k - 1] + Math.sqrt((xs[k] - xs[k - 1]) ** 2 + (ys[k] - ys[k - 1]) ** 2));
        // cumsum in R accumulates step by step, same as above
      }
      const z = hasEle && smooth_ele ? rows.map((i) => eleNow[i]) : null;
      const out = smoothSegment(xs, ys, z, tSeg, method, hasTime, { sigma_obs, sigma_acc, rts, alpha, bandwidth,
        sg_window, sg_order, window, jump_k, local_k, edge_jump, gap });
      rows.forEach((r, k) => {
        x[r] = out.x[k];
        y[r] = out.y[k];
        if (z) eleNow[r] = out.z[k];
      });
    }
    if (hasEle && smooth_ele) t.ele = eleNow;
  }

  const [lo2, la2] = toLonlat(x, y, epsg);
  t.lon = lo2;
  t.lat = la2;
  for (const f of FLAG_NAMES) t[f] = F[f];
  attrs(t).sanitize_report = report;
  return t;
}

// Linear interpolation of the NA values of y along x (ends take the nearest value)
export function fillLinear(x, y) {
  const ok = y.map((v, i) => !Number.isNaN(v) && !Number.isNaN(x[i]));
  const nOk = ok.filter(Boolean).length;
  if (nOk < 2 || nOk === y.length) return y.slice();
  const xo = x.filter((_, i) => ok[i]);
  const yo = y.filter((_, i) => ok[i]);
  const missing = [];
  ok.forEach((o, i) => { if (!o) missing.push(i); });
  const filled = rs.approx(xo, yo, missing.map((i) => x[i]), { ties: "mean", rule: 2 });
  const out = y.slice();
  missing.forEach((i, k) => { out[i] = filled[k]; });
  return out;
}

const diffPrepend = (v) => v.map((a, i) => (i === 0 ? NaN : a - v[i - 1]));

function smoothSegment(x, y, z, t, method, hasTime, p) {
  const reset = hasTime ? p.gap : Infinity;
  if (method === "kalman") {
    const k = kalmanCv([x, y], t, p.sigma_obs, p.sigma_acc, p.rts, reset);
    return { x: k[0], y: k[1], z: z ? kalmanCv([z], t, p.sigma_obs, p.sigma_acc, p.rts, reset)[0] : null };
  }
  if (method === "loess") {
    const fit = (v) => {
      const ok = v.map((a) => !Number.isNaN(a));
      return rs.loessDirect(t.filter((_, i) => ok[i]), v.filter((_, i) => ok[i]), p.alpha, t);
    };
    return { x: fit(x), y: fit(y), z: z ? fit(z) : null };
  }
  if (method === "gaussian") {
    if (!hasTime) throw new Error("method = \"gaussian\" needs a 'time' column.");
    return { x: gaussianKernel(t, x, p.bandwidth), y: gaussianKernel(t, y, p.bandwidth),
      z: z ? gaussianKernel(t, z, p.bandwidth) : null };
  }
  if (method === "savgol") {
    return { x: savgol(x, p.sg_window, p.sg_order), y: savgol(y, p.sg_window, p.sg_order),
      z: z ? savgol(z, p.sg_window, p.sg_order) : null };
  }
  // method === "pipeline" (Martinez, Sanchez & Diaz, 2026)
  const k = kalmanCv([x, y], t, p.sigma_obs, p.sigma_acc, false, reset);
  x = k[0];
  y = k[1];
  if (z) z = kalmanCv([z], t, p.sigma_obs, p.sigma_acc, false, reset)[0];
  const dt = diffPrepend(t);
  const dd = x.map((v, i) => (i === 0 ? NaN : Math.sqrt((v - x[i - 1]) ** 2 + (y[i] - y[i - 1]) ** 2)));
  const speed = dd.map((d, i) => d / dt[i]);
  const derived = [speed, diffPrepend(speed).map((v, i) => v / dt[i])];
  if (z) derived.push(diffPrepend(z).map((v, i) => v / dd[i]));
  const outlier = new Array(x.length).fill(false);
  for (let v of derived) {
    v = v.map((a) => (Number.isFinite(a) ? a : NaN));
    const q = rs.quantile7(v, [0.25, 0.75]);
    const iqr = q[1] - q[0];
    const dv = diffPrepend(v);
    const sdv = rs.rSd(dv);
    const lm = rs.convFilterCentered(v, p.window);
    const l2 = rs.convFilterCentered(v.map((a) => a * a), p.window);
    v.forEach((a, i) => {
      const c1 = a < q[0] - 1.5 * iqr || a > q[1] + 1.5 * iqr;
      const c2 = Math.abs(dv[i]) > p.jump_k * sdv;
      const lsd = Math.sqrt(Math.max(l2[i] - lm[i] ** 2, 0));
      const c3 = Math.abs(a - lm[i]) > p.local_k * lsd;
      if (c1 + c2 + c3 >= 2) outlier[i] = true;
    });
  }
  x = fillLinear(t, x.map((v, i) => (outlier[i] ? NaN : v)));
  y = fillLinear(t, y.map((v, i) => (outlier[i] ? NaN : v)));
  if (z) z = fillLinear(t, z.map((v, i) => (outlier[i] ? NaN : v)));
  x = savgol(x, p.sg_window, p.sg_order);
  y = savgol(y, p.sg_window, p.sg_order);
  if (z) z = savgol(z, p.sg_window, p.sg_order);
  const n = x.length;
  for (let i = 1; i <= 3; i++) {
    if (n > 7 && Math.sqrt((x[i] - x[i - 1]) ** 2 + (y[i] - y[i - 1]) ** 2) > p.edge_jump) {
      for (let k = 0; k < i; k++) { x[k] = x[i]; y[k] = y[i]; }
    }
    const j = n - i + 1;
    if (n > 7 && Math.sqrt((x[j - 1] - x[j - 2]) ** 2 + (y[j - 1] - y[j - 2]) ** 2) > p.edge_jump) {
      for (let k = j - 1; k < n; k++) { x[k] = x[j - 2]; y[k] = y[j - 2]; }
    }
  }
  return { x, y, z };
}

// Constant-velocity Kalman filter (and RTS smoother) for every axis in obs
export function kalmanCv(obs, t, sigmaObs, sigmaAcc, rts = true, reset = Infinity) {
  const n = t.length;
  const out = obs.map((a) => a.slice());
  let start = 0;
  for (let i = 1; i <= n; i++) {
    const dt = i < n ? t[i] - t[i - 1] : 0;
    if (i === n || Number.isNaN(dt) || dt < 0 || dt > reset) {
      const rows = [];
      for (let k = start; k < i; k++) rows.push(k);
      const block = kalmanBlock(obs.map((a) => rows.map((k) => a[k])), rows.map((k) => t[k]), sigmaObs, sigmaAcc, rts);
      block.forEach((a, ax) => rows.forEach((k, j) => { out[ax][k] = a[j]; }));
      start = i;
    }
  }
  return out;
}

// 2x2 helpers on [[a, b], [c, d]]
const mm = (A, B) => [[A[0][0] * B[0][0] + A[0][1] * B[1][0], A[0][0] * B[0][1] + A[0][1] * B[1][1]],
  [A[1][0] * B[0][0] + A[1][1] * B[1][0], A[1][0] * B[0][1] + A[1][1] * B[1][1]]];
const tr = (A) => [[A[0][0], A[1][0]], [A[0][1], A[1][1]]];
const madd = (A, B, s = 1) => [[A[0][0] + s * B[0][0], A[0][1] + s * B[0][1]], [A[1][0] + s * B[1][0], A[1][1] + s * B[1][1]]];
const mv = (A, v) => [A[0][0] * v[0] + A[0][1] * v[1], A[1][0] * v[0] + A[1][1] * v[1]];

function kalmanBlock(obs, t, sigmaObs, sigmaAcc, rts) {
  const n = t.length;
  const R = sigmaObs ** 2;
  return obs.map((z) => {
    const firstIdx = z.findIndex((v) => !Number.isNaN(v));
    if (firstIdx < 0 || n < 2) return z.slice();
    let m = [z[firstIdx], 0];
    let P = [[R, 0], [0, 100]];
    const mF = [];
    const PF = [];
    const mP = [];
    const PP = [];
    for (let i = 0; i < n; i++) {
      const dt = i === 0 ? 0 : Math.max(t[i] - t[i - 1], 0);
      const Fm = [[1, dt], [0, 1]];
      const s2 = sigmaAcc ** 2;
      const Q = [[s2 * (dt ** 4 / 4), s2 * (dt ** 3 / 2)], [s2 * (dt ** 3 / 2), s2 * dt ** 2]];
      m = mv(Fm, m);
      P = madd(mm(mm(Fm, P), tr(Fm)), Q);
      mP.push(m);
      PP.push(P);
      if (!Number.isNaN(z[i])) {
        const K = [P[0][0] / (P[0][0] + R), P[1][0] / (P[0][0] + R)];
        m = [m[0] + K[0] * (z[i] - m[0]), m[1] + K[1] * (z[i] - m[0])];
        P = [[P[0][0] - K[0] * P[0][0], P[0][1] - K[0] * P[0][1]], [P[1][0] - K[1] * P[0][0], P[1][1] - K[1] * P[0][1]]];
      }
      mF.push(m);
      PF.push(P);
    }
    let ms = mF;
    if (rts) {
      ms = mF.slice();
      const Ps = PF.slice();
      for (let i = n - 2; i >= 0; i--) {
        const dt = Math.max(t[i + 1] - t[i], 0);
        const Fm = [[1, dt], [0, 1]];
        const C = mm(mm(PF[i], tr(Fm)), pseudoInverse(PP[i + 1]));
        const diff = [ms[i + 1][0] - mP[i + 1][0], ms[i + 1][1] - mP[i + 1][1]];
        const corr = mv(C, diff);
        ms[i] = [mF[i][0] + corr[0], mF[i][1] + corr[1]];
        Ps[i] = madd(PF[i], mm(mm(C, madd(Ps[i + 1], PP[i + 1], -1)), tr(C)));
      }
    }
    return ms.map((v) => v[0]);
  });
}

// Moore-Penrose inverse of a symmetric 2x2 matrix (eigen decomposition)
export function pseudoInverse(A) {
  const a = A[0][0];
  const b = (A[0][1] + A[1][0]) / 2;
  const d = A[1][1];
  const mean = (a + d) / 2;
  const disc = Math.sqrt(((a - d) / 2) ** 2 + b * b);
  const l1 = mean + disc;
  const l2 = mean - disc;
  let v1;
  if (b !== 0) v1 = [l1 - d, b];
  else v1 = a >= d ? [1, 0] : [0, 1];
  const nv = Math.hypot(v1[0], v1[1]);
  v1 = [v1[0] / nv, v1[1] / nv];
  const v2 = [-v1[1], v1[0]];
  const sv = [Math.abs(l1), Math.abs(l2)];
  const smax = Math.max(sv[0], sv[1]);
  const tol = 2 * smax * rs.DBL_EPSILON;
  let out = [[0, 0], [0, 0]];
  [[l1, v1], [l2, v2]].forEach(([l, v], k) => {
    if (sv[k] > tol && l !== 0) {
      out = madd(out, [[v[0] * v[0] / l, v[0] * v[1] / l], [v[1] * v[0] / l, v[1] * v[1] / l]]);
    }
  });
  return out;
}

export function gaussianKernel(t, v, bw) {
  const tt = [];
  const vv = [];
  v.forEach((a, i) => { if (!Number.isNaN(a)) { tt.push(t[i]); vv.push(a); } });
  return t.map((ti) => {
    const w = [];
    const wv = [];
    for (let k = 0; k < tt.length; k++) {
      if (Math.abs(tt[k] - ti) <= 4 * bw) {
        const wk = rs.dnorm(tt[k] - ti, bw);
        w.push(wk);
        wv.push(wk * vv[k]);
      }
    }
    return rs.rSum(wv) / rs.rSum(w);
  });
}

// Savitzky-Golay smoothing; edges use the polynomials of the first and last windows
export function savgol(v, window, order) {
  const n = v.length;
  if (window % 2 === 0) window += 1;
  if (n < window) return v.slice();
  const half = (window - 1) / 2;
  // projection matrix X (X'X)^-1 X'
  const X = [];
  for (let i = -half; i <= half; i++) {
    const r = [];
    for (let j = 0; j <= order; j++) r.push(i ** j);
    X.push(r);
  }
  const proj = X.map(() => new Array(window).fill(0));
  for (let col = 0; col < window; col++) {
    const e = new Array(window).fill(0);
    e[col] = 1;
    const beta = rs.lstsq(X, e);
    for (let r = 0; r < window; r++) proj[r][col] = X[r].reduce((s, a, j) => s + a * beta[j], 0);
  }
  const out = v.slice();
  const coef = proj[half];
  for (let i = half; i < n - half; i++) {
    let s = 0;
    for (let k = 0; k < window; k++) s += coef[k] * v[i - half + k];
    out[i] = s;
  }
  for (let r = 0; r < half; r++) out[r] = proj[r].reduce((s, a, k) => s + a * v[k], 0);
  for (let r = half + 1; r < window; r++) out[n - window + r] = proj[r].reduce((s, a, k) => s + a * v[n - window + k], 0);
  return out;
}
