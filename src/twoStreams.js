// Comparing streams, linear referencing and map matching (port of
// R/compare_stream.R, R/detect_stream.R, R/project_stream.R, R/match_stream.R).

import * as rs from "./rstats.js";
import { asTable, attrs, checkColumns, nearestNeighbour, numCol, pathDistance, toLonlat, toMetric, utmEpsg,
  isNA } from "./utils.js";
import { SegmentIndex } from "./strtree.js";
import { tobler_speed } from "./tobler.js";

const rowsXY = (x, y, z = null) => x.map((v, i) => (z ? [v, y[i], z[i]] : [v, y[i]]));

export function compare_stream(a, b, { threshold = 1, dims = 2 } = {}) {
  const cols = dims === 3 ? ["lon", "lat", "ele"] : ["lon", "lat"];
  checkColumns(a, cols, "compare_stream");
  checkColumns(b, cols, "compare_stream");
  const epsg = utmEpsg(numCol(a.lon), numCol(a.lat));
  const [ax, ay] = toMetric(numCol(a.lon), numCol(a.lat), epsg);
  const [bx, by] = toMetric(numCol(b.lon), numCol(b.lat), epsg);
  let pa = rowsXY(ax, ay, dims === 3 ? numCol(a.ele) : null).filter((r) => !r.some(isNA));
  let pb = rowsXY(bx, by, dims === 3 ? numCol(b.ele) : null).filter((r) => !r.some(isNA));
  const dab = nearestNeighbour(pa, pb).dist;
  const dba = nearestNeighbour(pb, pa).dist;
  return {
    dims: [dims],
    hausdorff: [Math.max(Math.max(...dab), Math.max(...dba))],
    hausdorff_ab: [Math.max(...dab)],
    hausdorff_ba: [Math.max(...dba)],
    separation_ab: [rs.rMedian(dab)],
    separation_ba: [rs.rMedian(dba)],
    metric_score_ab: [100 * rs.rMean(dab.map((d) => (d <= threshold ? 1 : 0)))],
    metric_score_ba: [100 * rs.rMean(dba.map((d) => (d <= threshold ? 1 : 0)))],
  };
}

function emptyPasses(diagnosis) {
  const out = { pass: [], start_row: [], end_row: [], direction: [], coverage: [], status: [] };
  attrs(out).diagnosis = diagnosis;
  return out;
}

// R prints round(x, 1) with up to 15 significant digits
const rRound1 = (x) => String(Math.round(x * 10) / 10);

// Detect the passes of a stream along a route (rows are 1-based, as in R)
export function detect_stream(data, route, { tolerance = 50, min_completion = 1, min_coverage = 0.97, min_partial = 0.25 } = {}) {
  checkColumns(data, ["lon", "lat"], "detect_stream");
  checkColumns(route, ["lon", "lat"], "detect_stream");
  const epsg = utmEpsg(numCol(route.lon), numCol(route.lat));
  const [rx0, ry0] = toMetric(numCol(route.lon), numCol(route.lat), epsg);
  const [dx, dy] = toMetric(numCol(data.lon), numCol(data.lat), epsg);
  const rxy = rowsXY(rx0, ry0).filter((r) => !r.some(isNA));
  const dxy = rowsXY(dx, dy);
  const n = dxy.length;
  const m = rxy.length;
  const dData = nearestNeighbour(dxy, rxy);
  const dRoute = nearestNeighbour(rxy, dxy).dist;
  if (Math.min(...rs.naOmit(dData.dist)) > tolerance) return emptyPasses("streams do not intersect");
  const completion = rs.rMean(dRoute.map((d) => (d <= tolerance ? 1 : 0)));
  if (completion < min_completion) {
    return emptyPasses(`streams intersect but the stream covers only ${rRound1(100 * completion)}% of the route`);
  }
  const routeLength = rs.rSum(rxy.slice(1).map((r, i) => Math.sqrt((r[0] - rxy[i][0]) ** 2 + (r[1] - rxy[i][1]) ** 2)));
  const step = dxy.map((r, i) => {
    if (i === 0) return 0;
    const s = Math.sqrt((r[0] - dxy[i - 1][0]) ** 2 + (r[1] - dxy[i - 1][1]) ** 2);
    return Number.isNaN(s) ? 0 : s;
  });
  let acc = 0;
  const cum = step.map((s) => (acc += s));
  const nearest = dData.index;
  const ends = [rxy[0], rxy[m - 1]];
  const nearEnd = ends.map((e) => dxy.map((r) => Math.sqrt((r[0] - e[0]) ** 2 + (r[1] - e[1]) ** 2) <= tolerance));

  let passes = [];
  for (let e = 0; e < 2; e++) {
    const close = nearEnd[e];
    const target = nearEnd[1 - e];
    const starts = [];
    close.forEach((c, i) => { if (c && !(i > 0 && close[i - 1])) starts.push(i + 1); });
    for (const start of starts) {
      let end = -1;
      for (let i = start; i < n; i++) {
        if (target[i] && cum[i] - cum[start - 1] >= 0.5 * routeLength) { end = i + 1; break; }
      }
      if (end < 0) {
        for (let i = 0; i < n; i++) if (cum[i] <= cum[start - 1] + routeLength + tolerance) end = i + 1;
      }
      if (end - start < 2) continue;
      const sub = dxy.slice(start - 1, end);
      const covered = nearestNeighbour(rxy, sub).dist.map((d) => (d <= tolerance ? 1 : 0));
      const ranks = [];
      for (let k = 1; k <= end - start + 1; k++) ranks.push(k);
      const rho = rs.corSpearman(ranks, nearest.slice(start - 1, end));
      passes.push({ start_row: start, end_row: end, direction: !Number.isNaN(rho) && rho < 0 ? "reverse" : "forward",
        coverage: rs.rMean(covered) });
    }
  }
  passes = passes.filter((p) => p.coverage >= min_partial);
  passes.sort((a, b) => b.coverage - a.coverage || a.start_row - b.start_row);
  const keep = passes.map(() => true);
  passes.forEach((p, i) => {
    if (!keep[i]) return;
    passes.forEach((q, j) => {
      if (j > i && q.start_row <= p.end_row && q.end_row >= p.start_row) keep[j] = false;
    });
  });
  passes = passes.filter((_, i) => keep[i]).sort((a, b) => a.start_row - b.start_row);
  if (passes.length === 0) return emptyPasses("the stream never passes through the start or end of the route");
  const out = {
    pass: passes.map((_, i) => i + 1),
    start_row: passes.map((p) => p.start_row),
    end_row: passes.map((p) => p.end_row),
    direction: passes.map((p) => p.direction),
    coverage: passes.map((p) => p.coverage),
    status: passes.map((p) => (p.coverage >= min_coverage ? "match" : "partial")),
  };
  attrs(out).diagnosis = out.status.includes("match") ? "match" : "partial match only";
  return out;
}

export function project_stream(data, { epsg = null } = {}) {
  checkColumns(data, ["lon", "lat"], "project_stream");
  const t = asTable(data);
  const lon = numCol(t.lon);
  const lat = numCol(t.lat);
  if (epsg === null) epsg = utmEpsg(lon, lat);
  const [x, y] = toMetric(lon, lat, epsg);
  t.x = x;
  t.y = y;
  attrs(t).epsg = epsg;
  return t;
}

// Orthogonal projection on a polyline: distance along it and signed offset
export function projectOnLine(px, py, rx, ry, rsAlong) {
  const index = new SegmentIndex(rx, ry);
  const s = [];
  const cross = [];
  px.forEach((x, k) => {
    const y = py[k];
    const a = index.nearest(x, y);
    const dx = rx[a + 1] - rx[a];
    const dy = ry[a + 1] - ry[a];
    const len2 = dx * dx + dy * dy;
    let u = len2 > 0 ? ((x - rx[a]) * dx + (y - ry[a]) * dy) / len2 : 0;
    u = Math.min(Math.max(u, 0), 1);
    const qx = rx[a] + u * dx;
    const qy = ry[a] + u * dy;
    const d = Math.sqrt((x - qx) ** 2 + (y - qy) ** 2);
    const side = Math.sign(dx * (y - qy) - dy * (x - qx));
    s.push(rsAlong[a] + u * Math.sqrt(len2));
    cross.push(side * d);
  });
  return { s, cross };
}

export function locate_stream(data, reference) {
  checkColumns(data, ["lon", "lat"], "locate_stream");
  checkColumns(reference, ["lon", "lat"], "locate_stream");
  const t = asTable(data);
  const epsg = utmEpsg(numCol(reference.lon), numCol(reference.lat));
  const [rx, ry] = toMetric(numCol(reference.lon), numCol(reference.lat), epsg);
  if (rx.some(isNA) || ry.some(isNA)) throw new Error("The reference has missing coordinates.");
  if (rx.length < 2) throw new Error("The reference needs at least 2 readings.");
  const rsAlong = pathDistance(rx, ry);
  const [ox, oy] = toMetric(numCol(t.lon), numCol(t.lat), epsg);
  const ok = ox.map((v, i) => !isNA(v) && !isNA(oy[i]));
  const p = projectOnLine(ox.filter((_, i) => ok[i]), oy.filter((_, i) => ok[i]), rx, ry, rsAlong);
  const s = [];
  const cross = [];
  let k = 0;
  ok.forEach((o) => {
    if (o) { s.push(p.s[k]); cross.push(p.cross[k]); k++; } else { s.push(NaN); cross.push(NaN); }
  });
  const px = rs.approx(rsAlong, rx, s, { ties: "ordered" });
  const py = rs.approx(rsAlong, ry, s, { ties: "ordered" });
  const [lo, la] = toLonlat(px, py, epsg);
  t.s_ref = s;
  t.offset_cross = cross;
  t.lon_ref = lo;
  t.lat_ref = la;
  return t;
}

// Align a stream along a reference trail by dynamic programming (Viterbi)
export function match_stream(data, reference, params, { direction = "auto", by = null, dims = 2, lambda = 1,
  tol_factor = 1.1, slope_clip = 0.6, max_cells = 5e7 } = {}) {
  if (!["auto", "forward", "reverse"].includes(direction)) throw new Error("'direction' must be 'auto', 'forward' or 'reverse'");
  checkColumns(data, ["lon", "lat", "time", ...(dims === 3 ? ["ele"] : []), ...(by ? [by] : [])], "match_stream");
  checkColumns(reference, ["lon", "lat"], "match_stream");
  const t = asTable(data);
  const refLon = numCol(reference.lon);
  const refLat = numCol(reference.lat);
  const epsg = utmEpsg(refLon, refLat);
  const [rx, ry] = toMetric(refLon, refLat, epsg);
  const rz = "ele" in reference ? numCol(reference.ele) : rx.map(() => 0);
  if (rx.some(isNA) || ry.some(isNA) || rz.some(isNA)) throw new Error("The reference has missing coordinates or elevations.");
  const rsAlong = pathDistance(rx, ry);
  const m = rx.length;
  const n0 = t.lon.length;

  const groupKeys = new Map();
  const group = by === null ? new Array(n0).fill(1) : t[by].map((v) => {
    const k = String(v);
    if (!groupKeys.has(k)) groupKeys.set(k, groupKeys.size + 1);
    return groupKeys.get(k);
  });

  const out = asTable(t);
  out.lon_obs = t.lon.slice();
  out.lat_obs = t.lat.slice();
  if ("ele" in t) out.ele_obs = t.ele.slice();
  const newCols = ["ref_index", "s_ref", "offset_d", "offset_cross", "offset_along", "speed_h", "speed_v", "speed_total", "slope_ref"];
  for (const c of newCols) out[c] = new Array(n0).fill(NaN);
  out.ele = new Array(n0).fill(NaN);
  out.lon = numCol(out.lon);
  out.lat = numCol(out.lat);
  const lonAll = numCol(t.lon);
  const latAll = numCol(t.lat);
  const timeAll = numCol(t.time);
  const eleAll = "ele" in t ? numCol(t.ele) : null;
  const report = [];

  for (const g of [...new Set(group)]) {
    const rows = [];
    group.forEach((v, i) => { if (v === g) rows.push(i); });
    const [ox, oy] = toMetric(rows.map((i) => lonAll[i]), rows.map((i) => latAll[i]), epsg);
    const tt = rows.map((i) => timeAll[i]);
    const n = rows.length;
    if (n * m > max_cells) {
      throw new Error(`Round with ${n} readings and a reference of ${m} vertices exceeds 'max_cells'; split the stream into rounds or use a coarser reference.`);
    }
    const proj = projectOnLine(ox, oy, rx, ry, rsAlong);
    let dirG = direction;
    if (dirG === "auto") {
      const lo = Math.ceil(0.1 * n);
      const hi = Math.max(lo, Math.floor(0.9 * n));
      const inner = [];
      for (let k = lo; k <= hi; k++) inner.push(k);
      const rho = rs.corSpearman(inner.map((_, k) => k + 1), inner.map((k) => proj.s[k - 1]));
      dirG = !Number.isNaN(rho) && rho < 0 ? "reverse" : "forward";
    }
    const ord = [];
    for (let k = 0; k < m; k++) ord.push(dirG === "forward" ? k : m - 1 - k);
    const gx = ord.map((k) => rx[k]);
    const gy = ord.map((k) => ry[k]);
    const gz = ord.map((k) => rz[k]);
    const gs = dirG === "forward" ? rsAlong.slice() : ord.map((k) => rsAlong[m - 1] - rsAlong[k]);
    const oz = dims === 3 ? rows.map((i) => eleAll[i]) : null;
    const path = viterbiMatch(ox, oy, oz, tt, gx, gy, gz, gs, params, lambda, tol_factor, slope_clip);
    const j = path.j;
    const jRef = j.map((k) => ord[k]);

    const ds = j.map((k, i) => (i === 0 ? NaN : gs[k] - gs[j[i - 1]]));
    const dz = j.map((k, i) => (i === 0 ? NaN : gz[k] - gz[j[i - 1]]));
    const dt = tt.map((v, i) => (i === 0 ? NaN : v - tt[i - 1]));
    const ok = dt.map((v) => !Number.isNaN(v) && v > 0);
    rows.forEach((r, i) => {
      const jr = jRef[i];
      out.lon[r] = refLon[jr];
      out.lat[r] = refLat[jr];
      out.ele[r] = rz[jr];
      out.ref_index[r] = jr + 1;
      out.s_ref[r] = rsAlong[jr];
      out.offset_d[r] = Math.sqrt((ox[i] - rx[jr]) ** 2 + (oy[i] - ry[jr]) ** 2);
      const pj = Math.max(j[i] - 1, 0);
      const nj = Math.min(j[i] + 1, m - 1);
      let tx = gx[nj] - gx[pj];
      let ty = gy[nj] - gy[pj];
      let norm = Math.sqrt(tx * tx + ty * ty);
      if (norm === 0) norm = NaN;
      tx /= norm;
      ty /= norm;
      const ex = ox[i] - gx[j[i]];
      const ey = oy[i] - gy[j[i]];
      out.offset_along[r] = ex * tx + ey * ty;
      out.offset_cross[r] = tx * ey - ty * ex;
      out.speed_h[r] = ok[i] ? ds[i] / dt[i] : NaN;
      out.speed_v[r] = ok[i] ? Math.abs(dz[i]) / dt[i] : NaN;
      out.speed_total[r] = ok[i] ? Math.sqrt(ds[i] ** 2 + dz[i] ** 2) / dt[i] : NaN;
      out.slope_ref[r] = !Number.isNaN(ds[i]) && ds[i] > 0 ? dz[i] / ds[i] : NaN;
    });
    const load = new Array(m).fill(0);
    for (const k of j) load[k]++;
    const cap = Math.max(2, Math.ceil(n / m));
    const slopeC = ds.map((d, i) => Math.min(Math.max(!Number.isNaN(d) && d > 0 ? dz[i] / d : 0, -slope_clip), slope_clip));
    const ceilH = tobler_speed(slopeC, params).map((v) => v * tol_factor);
    const atCeiling = ds.map((d, i) => ok[i] && d > 0 && d / dt[i] >= 0.98 * ceilH[i]);
    const offs = rows.map((r) => out.offset_d[r]);
    report.push({
      group: by === null ? null : t[by][rows[0]],
      direction: dirG,
      n,
      vertices_used: load.filter((v) => v > 0).length,
      max_per_vertex: Math.max(...load),
      vertices_over_cap: load.filter((v) => v > cap).length,
      pct_at_ceiling: 100 * rs.rMean(atCeiling.filter((_, i) => ok[i]).map(Number)),
      objective: path.objective,
      offset_median: rs.rMedian(offs),
      offset_q95: rs.quantile7(offs, 0.95),
    });
  }
  attrs(out).match_report = report;
  return out;
}

export function viterbiMatch(ox, oy, oz, t, gx, gy, gz, gs, params, lambda, tolFactor, slopeClip) {
  const n = ox.length;
  const m = gx.length;
  const slopes = rs.seqBy(-slopeClip, slopeClip, 0.01);
  const vt = tobler_speed(slopes, params);
  const vTop = Math.max(...slopes.map((s, i) => vt[i] * tolFactor * Math.sqrt(1 + s * s)));
  const unary = (i) => gx.map((x, k) => {
    let u = (ox[i] - x) ** 2 + (oy[i] - gy[k]) ** 2;
    if (oz) u += (oz[i] - gz[k]) ** 2;
    return u;
  });
  let D = unary(0);
  const back = [new Int32Array(m)];
  for (let i = 1; i < n; i++) {
    const dt = t[i] - t[i - 1];
    const best = D.slice();
    const arg = new Int32Array(m);
    for (let k = 0; k < m; k++) arg[k] = k;
    if (!Number.isNaN(dt) && dt > 0) {
      const reach = rs.findInterval(gs.map((s) => s + vTop * dt), gs);
      let kMax = -Infinity;
      reach.forEach((r, k) => { if (r - (k + 1) > kMax) kMax = r - (k + 1); });
      for (let k = 1; k <= kMax; k++) {
        for (let from = 0; from < m - k; from++) {
          const to = from + k;
          const ds = gs[to] - gs[from];
          const dz = gz[to] - gz[from];
          const slope = ds > 0 ? dz / ds : 0;
          const slopeC = Math.min(Math.max(slope, -slopeClip), slopeClip);
          const vhC = tobler_speed(slopeC, params) * tolFactor;
          const allowed = ds === 0 || Math.sqrt(ds * ds + dz * dz) / dt <= vhC * Math.sqrt(1 + slopeC * slopeC) + 1e-9;
          if (!allowed) continue;
          let pen = 0;
          if (ds !== 0) {
            const vh = tobler_speed(slope, params) * tolFactor;
            const excess = Math.max(ds / dt - vh, 0);
            pen = (excess ** 2 * (1 + slope ** 2 + (1 + slope ** 2))) / lambda ** 2;
          }
          const cand = D[from] + pen;
          if (cand < best[to]) {
            best[to] = cand;
            arg[to] = from;
          }
        }
      }
    }
    const u = unary(i);
    D = best.map((b, k) => b + u[k]);
    back.push(arg);
  }
  const j = new Array(n);
  let bestK = 0;
  for (let k = 1; k < m; k++) if (D[k] < D[bestK]) bestK = k;
  j[n - 1] = bestK;
  for (let i = n - 1; i > 0; i--) j[i - 1] = back[i][j[i]];
  return { j, objective: D[bestK] };
}
