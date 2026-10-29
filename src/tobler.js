// Tobler hiking function (port of R/tobler_stream.R).

import * as rs from "./rstats.js";
import { checkColumns } from "./utils.js";
import { differential_stream } from "./measures.js";

const PARAM_NAMES = ["v1", "m1", "c_plus", "c_minus"];

// Tobler parameters from tobler_stream() output, an object or an array
export function toblerParams(params) {
  if (Array.isArray(params)) return params.slice(0, 4).map(Number);
  const get = (k) => {
    const v = params[k];
    return Array.isArray(v) ? Number(v[0]) : Number(v);
  };
  const p = PARAM_NAMES.map(get);
  if (p.every((v) => !Number.isNaN(v))) return p;
  return Object.values(params).slice(0, 4).map((v) => (Array.isArray(v) ? Number(v[0]) : Number(v)));
}

// Speed predicted by the Tobler function at slope (rise over run)
export function tobler_speed(slope, params) {
  const p = toblerParams(params);
  const f = (s) => (s >= p[1] ? p[0] * Math.exp(-p[2] * (s - p[1])) : p[0] * Math.exp(-p[3] * (s - p[1])));
  return Array.isArray(slope) ? slope.map(f) : f(slope);
}

// Fit the Tobler hiking function as the tau-quantile frontier of speed against slope
export function tobler_stream(data, { tau = 0.9, method = "pinball", slope_range = [-0.6, 0.6], min_speed = 0.1,
  bin_width = 0.02, min_bin = 10 } = {}) {
  if (!["pinball", "binned"].includes(method)) throw new Error("'method' must be 'pinball' or 'binned'");
  checkColumns(data, ["lon", "lat", "ele", "time"], "tobler_stream");
  const d = differential_stream(data);
  const v = [];
  const m = [];
  d.speed.forEach((s, i) => {
    const sl = d.slope[i];
    if (Number.isFinite(s) && Number.isFinite(sl) && s > min_speed && sl >= slope_range[0] && sl <= slope_range[1]) {
      v.push(s);
      m.push(sl);
    }
  });
  if (v.length < 30) throw new Error("Fewer than 30 usable readings to fit the Tobler function.");

  const breaks = rs.seqBy(Math.floor(Math.min(...m) / bin_width) * bin_width,
    Math.ceil(Math.max(...m) / bin_width) * bin_width + bin_width, bin_width);
  const code = rs.bincode(m, breaks);
  const bm = [];
  const bv = [];
  const bn = [];
  for (let b = 1; b < breaks.length; b++) {
    const idx = [];
    code.forEach((c, i) => { if (c === b) idx.push(i); });
    if (idx.length >= min_bin) {
      bm.push(rs.rMedian(idx.map((i) => m[i])));
      bv.push(rs.quantile7(idx.map((i) => v[i]), tau));
      bn.push(idx.length);
    }
  }
  if (bm.length < 4) throw new Error("Not enough populated slope bins; lower 'min_bin' or widen 'bin_width'.");
  const logv = bv.map(Math.log);
  const fitAt = (m1) => rs.wls(bm.map((x) => [1, Math.max(x - m1, 0), Math.min(x - m1, 0)]), logv, bn);
  const rss = (m1) => rs.rSum(fitAt(m1)[1].map((r) => r * r));
  const inner = rs.quantile7(bm, [0.1, 0.9]);
  const m1 = rs.optimizeBrent(rss, inner[0], inner[1]);
  const coef = fitAt(m1)[0];
  let theta = [Math.exp(coef[0]), m1, -coef[1], -coef[2]];

  const pinball = (r) => rs.rSum(r.map((x) => x * (tau - (x < 0 ? 1 : 0))));
  const resid = (p) => {
    const pred = tobler_speed(m, p);
    return v.map((x, i) => x - pred[i]);
  };
  if (method === "pinball") {
    const vCap = 2 * Math.max(...v);
    const loss = (p) => {
      if (p[0] <= 0 || p[0] > vCap || p[1] < slope_range[0] || p[1] > slope_range[1]) return Infinity;
      return pinball(resid(p));
    };
    const starts = [theta, theta.map((x, i) => x * [1.1, 1, 1.2, 1.2][i]), theta.map((x, i) => x * [0.9, 1, 0.8, 0.8][i]),
      theta.map((x, i) => x + [0, 0.03, 0, 0][i]), theta.map((x, i) => x - [0, 0.03, 0, 0][i])];
    let best = null;
    for (const st of starts) {
      const opt = rs.nelderMead(loss, st, { maxit: 4000, reltol: 1e-10 });
      if (best === null || opt.value < best.value) best = opt;
    }
    theta = best.par;
  }
  if (theta[0] > 10) {
    console.warn(`v1 = ${Math.round(theta[0] * 10) / 10} m/s is implausible for human locomotion; check the time of the stream.`);
  }
  const q = rs.quantile7(v, tau);
  const r1 = 1 - pinball(resid(theta)) / pinball(v.map((x) => x - q));
  return { v1: [theta[0]], m1: [theta[1]], c_plus: [theta[2]], c_minus: [theta[3]], tau: [tau], n: [v.length], r1: [r1], method: [method] };
}
