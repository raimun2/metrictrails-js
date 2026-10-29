// Spatial resolution and resampling (port of R/resample_stream.R).

import * as rs from "./rstats.js";
import { asTable, attrs, checkColumns, columnKind, isFirstOfSegment, names, numCol, pathDistance, reorder,
  segmentId, toLonlat, toMetric, uniqueInOrder, utmEpsg } from "./utils.js";

// Distribution of the distance between consecutive readings and the Nyquist
// resolution. `data` can be a stream, an array of streams or an object of named streams.
export function resolution_stream(data, { q = 0.9, idle = 0.1 } = {}) {
  const isStream = (d) => d && !Array.isArray(d) && "lon" in d && Array.isArray(d.lon);
  if (!isStream(data)) {
    const entries = Array.isArray(data) ? data.map((d, i) => [i + 1, d]) : Object.entries(data);
    const out = { stream: [] };
    for (const [key, d] of entries) {
      const row = resolution_stream(d, { q, idle });
      out.stream.push(key);
      for (const k of names(row)) (out[k] ||= []).push(row[k][0]);
    }
    return out;
  }
  checkColumns(data, ["lon", "lat"], "resolution_stream");
  const lon = numCol(data.lon);
  const lat = numCol(data.lat);
  const [x, y] = toMetric(lon, lat, utmEpsg(lon, lat));
  const seg = segmentId(data);
  const first = isFirstOfSegment(seg);
  const deltaD = [];
  for (let i = 1; i < x.length; i++) {
    if (!first[i]) deltaD.push(Math.sqrt((x[i] - x[i - 1]) ** 2 + (y[i] - y[i - 1]) ** 2));
  }
  const moving = deltaD.filter((d) => d >= idle && !Number.isNaN(d));
  if (moving.length === 0) throw new Error("The stream does not move.");
  const qs = rs.quantile7(moving, [0.05, 0.25, 0.5, 0.75, 0.95]);
  const valid = deltaD.filter((d) => !Number.isNaN(d));
  return {
    n: [lon.length],
    pct_idle: [100 * rs.rMean(valid.map((d) => (d < idle ? 1 : 0)))],
    mean: [rs.rMean(moving)],
    sd: [rs.rSd(moving)],
    q05: [qs[0]], q25: [qs[1]], q50: [qs[2]], q75: [qs[3]], q95: [qs[4]],
    max: [Math.max(...moving)],
    q: [q],
    res: [rs.quantile7(moving, q)],
  };
}

// indices (0-based) kept by GEOS DouglasPeuckerLineSimplifier
export function douglasPeucker(x, y, tol) {
  const n = x.length;
  const use = new Array(n).fill(true);
  const stack = [[0, n - 1]];
  while (stack.length) {
    const [i, j] = stack.pop();
    if (i + 1 >= j) continue;
    const ax = x[i], ay = y[i], bx = x[j], by = y[j];
    let maxD = -1;
    let maxK = i;
    for (let k = i + 1; k < j; k++) {
      const px = x[k], py = y[k];
      let d;
      if (ax === bx && ay === by) d = Math.sqrt((px - ax) ** 2 + (py - ay) ** 2);
      else {
        const len2 = (bx - ax) * (bx - ax) + (by - ay) * (by - ay);
        const r = ((px - ax) * (bx - ax) + (py - ay) * (by - ay)) / len2;
        if (r <= 0) d = Math.sqrt((px - ax) ** 2 + (py - ay) ** 2);
        else if (r >= 1) d = Math.sqrt((px - bx) ** 2 + (py - by) ** 2);
        else d = Math.abs(((ay - py) * (bx - ax) - (ax - px) * (by - ay)) / len2) * Math.sqrt(len2);
      }
      if (d > maxD) { maxD = d; maxK = k; }
    }
    if (maxD <= tol) for (let k = i + 1; k < j; k++) use[k] = false;
    else {
      stack.push([maxK, j]);
      stack.push([i, maxK]);
    }
  }
  const out = [];
  use.forEach((u, k) => { if (u) out.push(k); });
  return out;
}

// Resample each segment of a stream to a uniform spatial resolution
export function resample_stream(data, { res = "nyquist", q = 0.9, method = "arclength", idle = 0.1 } = {}) {
  if (!["arclength", "simplify"].includes(method)) throw new Error("'method' must be 'arclength' or 'simplify'");
  checkColumns(data, ["lon", "lat"], "resample_stream");
  const t = asTable(data);
  if (res === "nyquist") res = resolution_stream(t, { q, idle }).res[0];
  if (typeof res !== "number" || !(res > 0)) throw new Error("'res' must be a positive number or \"nyquist\".");

  const lon = numCol(t.lon);
  const lat = numCol(t.lat);
  const epsg = utmEpsg(lon, lat);
  const [X, Y] = toMetric(lon, lat, epsg);
  const seg = segmentId(t);
  const cols = names(t);
  const kinds = Object.fromEntries(cols.map((c) => [c, columnKind(t[c])]));
  const numericCols = cols.filter((c) => kinds[c] === "number" && !["lon", "lat", "segment"].includes(c));
  const timeCols = cols.filter((c) => kinds[c] === "time");
  const otherCols = cols.filter((c) => !["lon", "lat"].includes(c) && !numericCols.includes(c) && !timeCols.includes(c));

  const out = { x: [], y: [] };
  for (const c of [...numericCols, ...timeCols, ...otherCols]) out[c] = [];
  let pieces = 0;
  for (const s of uniqueInOrder(seg)) {
    const rows = [];
    seg.forEach((v, i) => { if (v === s) rows.push(i); });
    const x = rows.map((i) => X[i]);
    const y = rows.map((i) => Y[i]);
    const sOrig = pathDistance(x, y);
    const total = sOrig[sOrig.length - 1];
    if (rows.length < 2 || total < res) continue;
    let sNew;
    if (method === "arclength") {
      sNew = rs.seqLength(0, total, Math.floor(total / res) + 1);
    } else {
      const sDense = rs.seqLength(0, total, Math.ceil(total / 0.2) + 1);
      const xd = rs.approx(sOrig, x, sDense, { ties: "mean" });
      const yd = rs.approx(sOrig, y, sDense, { ties: "mean" });
      const kept = douglasPeucker(xd, yd, res);
      let pointer = 0;
      sNew = kept.map((kk) => {
        while (pointer < xd.length - 1 && (xd[pointer] !== xd[kk] || yd[pointer] !== yd[kk])) pointer++;
        return sDense[pointer];
      });
    }
    out.x.push(...rs.approx(sOrig, x, sNew, { ties: "mean" }));
    out.y.push(...rs.approx(sOrig, y, sNew, { ties: "mean" }));
    for (const c of numericCols) {
      const v = rows.map((i) => t[c][i]);
      const ok = v.map((a) => !Number.isNaN(a));
      const nOk = ok.filter(Boolean).length;
      out[c].push(...(nOk >= 2 ? rs.approx(sOrig.filter((_, i) => ok[i]), v.filter((_, i) => ok[i]), sNew, { ties: "mean", rule: 2 })
        : sNew.map(() => NaN)));
    }
    for (const c of timeCols) {
      const v = rows.map((i) => (t[c][i] === null ? NaN : t[c][i].getTime() / 1000));
      const ok = v.map((a) => !Number.isNaN(a));
      const vals = ok.filter(Boolean).length >= 2
        ? rs.approx(sOrig.filter((_, i) => ok[i]), v.filter((_, i) => ok[i]), sNew, { ties: "mean", rule: 2 })
        : sNew.map(() => NaN);
      out[c].push(...vals.map((secs) => (Number.isNaN(secs) ? null : new Date(Math.round(secs * 1000)))));
    }
    const before = rs.findInterval(sNew, sOrig, true).map((b) => Math.max(b, 1) - 1);
    for (const c of otherCols) {
      if (c === "segment") out[c].push(...sNew.map(() => t.segment[rows[0]]));
      else out[c].push(...before.map((b) => t[c][rows[b]]));
    }
    pieces++;
  }
  if (pieces === 0) throw new Error("No segment is longer than 'res'.");
  const [lo, la] = toLonlat(out.x, out.y, epsg);
  const result = { lon: lo, lat: la };
  for (const c of names(out)) if (c !== "x" && c !== "y") result[c] = out[c];
  const final = reorder(result, ["lon", "lat", "ele", "timestamp", "time"]);
  attrs(final).res = res;
  return final;
}
