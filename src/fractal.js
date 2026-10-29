// Fractal dimension by the divider method (port of R/fractal_stream.R).

import * as rs from "./rstats.js";
import { checkColumns, numCol, segmentId, takeRows, toMetric, utmEpsg, asTable } from "./utils.js";

export function divider_walk(x, y, r) {
  const nPts = x.length;
  let cx = x[0];
  let cy = y[0];
  let i = 0;
  let u0 = 0;
  let n = 0;
  const r2 = r * r * (1 - 1e-9);
  for (;;) {
    let k = i;
    while (k < nPts - 1 && (x[k + 1] - cx) ** 2 + (y[k + 1] - cy) ** 2 < r2) k++;
    if (k >= nPts - 1) return { n, length: n * r + Math.sqrt((x[nPts - 1] - cx) ** 2 + (y[nPts - 1] - cy) ** 2) };
    const dx = x[k + 1] - x[k];
    const dy = y[k + 1] - y[k];
    const fx = x[k] - cx;
    const fy = y[k] - cy;
    const a = dx * dx + dy * dy;
    const b = 2 * (fx * dx + fy * dy);
    const c = fx * fx + fy * fy - r * r;
    let u = (-b + Math.sqrt(Math.max(b * b - 4 * a * c, 0))) / (2 * a);
    const lower = k === i ? u0 : 0;
    u = Math.min(Math.max(u, lower), 1);
    cx = x[k] + u * dx;
    cy = y[k] + u * dy;
    i = k;
    u0 = u;
    n++;
  }
}

// Divider (Richardson) estimate of the fractal dimension of the 2D path
export function fractal_stream(data, { scales = null, n_scales = 10 } = {}) {
  checkColumns(data, ["lon", "lat"], "fractal_stream");
  let t = asTable(data);
  const lon0 = numCol(t.lon);
  const lat0 = numCol(t.lat);
  const keep = [];
  lon0.forEach((v, i) => { if (!Number.isNaN(v) && !Number.isNaN(lat0[i])) keep.push(i); });
  t = takeRows(t, keep);
  const lon = numCol(t.lon);
  const lat = numCol(t.lat);
  const [x, y] = toMetric(lon, lat, utmEpsg(lon, lat));
  const seg = segmentId(t);

  if (scales === null) {
    const step = [];
    for (let i = 1; i < x.length; i++) {
      if (seg[i] === seg[i - 1]) step.push(Math.sqrt((x[i] - x[i - 1]) ** 2 + (y[i] - y[i - 1]) ** 2));
    }
    const total = rs.rSum(step);
    const moving = step.filter((s) => s >= 0.1);
    if (moving.length < 3) throw new Error("The stream does not move.");
    const rMin = 2 * rs.rMedian(moving);
    const rMax = total / 10;
    if (rMax <= rMin) throw new Error("The path is too short for its resolution to estimate a fractal dimension.");
    scales = rs.seqLength(Math.log(rMin), Math.log(rMax), n_scales).map(Math.exp);
  }
  const nSteps = scales.map(() => 0);
  const lengths = scales.map(() => 0);
  const segs = [...new Set(seg)];
  scales.forEach((r, k) => {
    for (const s of segs) {
      const rows = [];
      seg.forEach((v, i) => { if (v === s) rows.push(i); });
      if (rows.length < 2) continue;
      const w = divider_walk(rows.map((i) => x[i]), rows.map((i) => y[i]), r);
      nSteps[k] += w.n;
      lengths[k] += w.length;
    }
  });
  const lx = scales.map(Math.log);
  const ly = lengths.map(Math.log);
  const [beta] = rs.wls(lx.map((v) => [1, v]), ly);
  return {
    dimension: 1 - beta[1],
    r_squared: rs.rCor(lx, ly) ** 2,
    scaling: { r: scales, n: nSteps, length: lengths },
  };
}
