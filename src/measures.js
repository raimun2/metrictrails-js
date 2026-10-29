// Derived measures (port of R/differential_stream.R and R/summarise_stream.R).

import * as rs from "./rstats.js";
import { asTable, checkColumns, diffBySegment, isFirstOfSegment, names, numCol, segmentId, takeRows,
  toMetric, utmEpsg } from "./utils.js";
import { fractal_stream } from "./fractal.js";

// R's %% for doubles
const rmod = (a, b) => a - Math.floor(a / b) * b;

export function differential_stream(data, { distance = "utm", idle = 0.1 } = {}) {
  if (!["utm", "haversine"].includes(distance)) throw new Error("'distance' must be 'utm' or 'haversine'");
  checkColumns(data, ["lon", "lat"], "differential_stream");
  const t = asTable(data);
  const seg = segmentId(t);
  const n = seg.length;
  const first = isFirstOfSegment(seg);
  const lon = numCol(t.lon);
  const lat = numCol(t.lat);
  const [x, y] = toMetric(lon, lat, utmEpsg(lon, lat));
  const dx = diffBySegment(x, seg);
  const dy = diffBySegment(y, seg);

  let deltaD;
  if (distance === "utm") deltaD = dx.map((v, i) => Math.sqrt(v * v + dy[i] * dy[i]));
  else {
    const rad = Math.PI / 180;
    deltaD = lat.map((la, i) => {
      if (first[i] || i === 0) return 0;
      const lat1 = lat[i - 1] * rad;
      const lat2 = la * rad;
      const dlat = lat2 - lat1;
      const dlon = (lon[i] - lon[i - 1]) * rad;
      const a = Math.sin(dlat / 2) ** 2 + Math.cos(lat1) * Math.cos(lat2) * Math.sin(dlon / 2) ** 2;
      return 2 * 6371008.8 * Math.asin(Math.min(1, Math.sqrt(a)));
    });
  }
  t.delta_d = deltaD;
  let acc = 0;
  t.distance = deltaD.map((v) => (acc += v));

  const moving = deltaD.map((v) => v >= idle);
  const heading = dx.map((v, i) => (moving[i] && !first[i] ? Math.atan2(dy[i], v) : NaN));
  t.turn_angle = heading.map((h, i) => {
    if (first[i] || i === 0) return NaN;
    const turn = rmod(h - heading[i - 1] + Math.PI, 2 * Math.PI) - Math.PI;
    return (Math.abs(turn) * 180) / Math.PI;
  });

  const hasEle = "ele" in t;
  let dEle;
  if (hasEle) {
    dEle = diffBySegment(numCol(t.ele), seg);
    t.delta_ele = dEle;
    t.dplus = dEle.map((v) => (Number.isNaN(v) ? NaN : Math.max(v, 0)));
    t.dminus = dEle.map((v) => (Number.isNaN(v) ? NaN : Math.max(-v, 0)));
    t.delta_d3d = dEle.map((v, i) => Math.sqrt(deltaD[i] ** 2 + v * v));
    t.slope = dEle.map((v, i) => (moving[i] && !first[i] ? v / deltaD[i] : NaN));
  }
  if ("time" in t) {
    const dt = diffBySegment(numCol(t.time), seg);
    t.delta_time = dt;
    const ok = dt.map((v, i) => !first[i] && v > 0);
    const speed = deltaD.map((v, i) => (ok[i] ? v / dt[i] : NaN));
    t.speed = speed;
    const ds = diffBySegment(speed, seg, NaN);
    t.accel = ds.map((v, i) => (ok[i] ? v / dt[i] : NaN));
    t.pace = speed.map((s) => (s > 0 ? 1000 / 60 / s : NaN));
    if (hasEle) {
      const vert = dEle.map((v, i) => (ok[i] ? v / dt[i] : NaN));
      t.vert_speed = vert;
      const dv = diffBySegment(vert, seg, NaN);
      t.vert_accel = dv.map((v, i) => (ok[i] ? v / dt[i] : NaN));
      t.speed_3d = t.delta_d3d.map((v, i) => (ok[i] ? v / dt[i] : NaN));
    }
  }
  return t;
}

export function summarise_stream(data, { by = null, idle = 0.1, fractal = false } = {}) {
  const byCols = by === null ? [] : Array.isArray(by) ? by : [by];
  checkColumns(data, ["lon", "lat", ...byCols], "summarise_stream");
  const t = asTable(data);
  const n = t.lon.length;
  const keys = new Map();
  const group = [];
  for (let i = 0; i < n; i++) {
    const k = byCols.map((b) => String(t[b][i])).join("\r");
    if (!keys.has(k)) keys.set(k, keys.size + 1);
    group.push(keys.get(k));
  }
  const rows = [];
  for (const g of [...new Set(group)]) {
    const idx = [];
    group.forEach((v, i) => { if (v === g) idx.push(i); });
    const sub = takeRows(t, idx);
    const d = differential_stream(sub, { idle });
    const isIdle = d.delta_d.map((v, i) => i > 0 && v < idle);
    const lon = numCol(sub.lon);
    const lat = numCol(sub.lat);
    const [x, y] = toMetric(lon, lat, utmEpsg(lon, lat));
    const k = x.length - 1;
    const straight = Math.sqrt((x[k] - x[0]) ** 2 + (y[k] - y[0]) ** 2);
    const dist2 = rs.rSum(d.delta_d);
    const row = {};
    for (const b of byCols) row[b] = sub[b][0];
    Object.assign(row, {
      n: idx.length,
      distance_2d: dist2,
      resolution: rs.rMean(d.delta_d.slice(1)),
      sinuosity: dist2 / straight,
      mean_turn: rs.rMean(rs.naOmit(d.turn_angle)),
      pct_idle: 100 * rs.rMean(isIdle.slice(1).map(Number)),
    });
    if ("ele" in sub) {
      const ele = rs.naOmit(numCol(sub.ele));
      row.distance_3d = rs.rSum(rs.naOmit(d.delta_d3d));
      row.elevation_gain = rs.rSum(rs.naOmit(d.dplus));
      row.elevation_loss = rs.rSum(rs.naOmit(d.dminus));
      row.ele_min = Math.min(...ele);
      row.ele_max = Math.max(...ele);
      row.mean_slope = 100 * rs.rMean(rs.naOmit(d.slope));
      row.iqr_slope = 100 * rs.rIqr(d.slope);
      row.km_effort = row.distance_2d / 1000 + row.elevation_gain / 100;
      row.idle_gain = rs.rSum(rs.naOmit(d.dplus.filter((_, i) => isIdle[i])));
      row.idle_loss = rs.rSum(rs.naOmit(d.dminus.filter((_, i) => isIdle[i])));
    }
    if ("time" in sub) {
      const tt = rs.naOmit(numCol(sub.time));
      row.duration = Math.max(...tt) - Math.min(...tt);
    }
    if (fractal) {
      try {
        row.fractal_dim = fractal_stream(sub).dimension;
      } catch (e) {
        row.fractal_dim = NaN;
      }
    }
    rows.push(row);
  }
  const out = {};
  for (const r of rows) for (const k of Object.keys(r)) if (!(k in out)) out[k] = [];
  for (const k of Object.keys(out)) out[k] = rows.map((r) => r[k]);
  return out;
}

export { names };
