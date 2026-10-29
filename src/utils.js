// Internal helpers (port of R/utils.R and _utils.py).
//
// A stream is a plain object whose keys are column names (in order) and whose
// values are arrays of equal length: numbers (NaN for NA), booleans, strings
// or null, and Date objects for `timestamp`. Attributes that R stores with
// attr() (sanitize_report, res, epsg, ...) are kept in a non-enumerable
// `attrs` property.

import proj4 from "proj4";
import { naOmit } from "./rstats.js";
import { PointIndex } from "./strtree.js";

export const NAMING_HINT = ". See the documentation of uniform_stream for the naming convention.";

// ---- tables -------------------------------------------------------------------------

export function names(t) {
  return Object.keys(t);
}

export function nrow(t) {
  const n = names(t);
  return n.length ? t[n[0]].length : 0;
}

export function attrs(t) {
  if (!Object.prototype.hasOwnProperty.call(t, "attrs")) {
    Object.defineProperty(t, "attrs", { value: {}, enumerable: false, writable: true });
  }
  return t.attrs;
}

export function copyTable(t) {
  const out = {};
  for (const k of names(t)) out[k] = t[k].slice();
  return out;
}

export function takeRows(t, idx) {
  const out = {};
  for (const k of names(t)) out[k] = idx.map((i) => t[k][i]);
  return out;
}

export function filterRows(t, keep) {
  const idx = [];
  keep.forEach((k, i) => { if (k) idx.push(i); });
  return takeRows(t, idx);
}

// rebuild a table with `first` columns in front, keeping the rest in order
export function reorder(t, first) {
  const out = {};
  const present = first.filter((c) => c in t);
  for (const c of present) out[c] = t[c];
  for (const c of names(t)) if (!present.includes(c)) out[c] = t[c];
  return out;
}

// table from an array of row objects (or a table already)
export function asTable(data) {
  if (Array.isArray(data)) {
    const out = {};
    for (const row of data) for (const k of Object.keys(row)) if (!(k in out)) out[k] = [];
    for (const k of Object.keys(out)) out[k] = data.map((r) => (r[k] === undefined ? null : r[k]));
    return out;
  }
  return copyTable(data);
}

export function toRows(t) {
  const n = nrow(t);
  const cols = names(t);
  const out = [];
  for (let i = 0; i < n; i++) {
    const r = {};
    for (const c of cols) r[c] = t[c][i];
    out.push(r);
  }
  return out;
}

export const isNA = (v) => v === null || v === undefined || (typeof v === "number" && Number.isNaN(v));

export function toNumber(v) {
  if (v === null || v === undefined) return NaN;
  if (typeof v === "number") return v;
  if (typeof v === "boolean") return v ? 1 : 0;
  if (v instanceof Date) return v.getTime() / 1000;
  const s = String(v).trim();
  if (s === "" || s === "NA") return NaN;
  const x = Number(s);
  return Number.isNaN(x) ? NaN : x;
}

export const numCol = (v) => Array.from(v, toNumber);

export function columnKind(v) {
  let kind = null;
  for (const x of v) {
    if (x === null || x === undefined || (typeof x === "number" && Number.isNaN(x))) continue;
    const k = x instanceof Date ? "time" : typeof x;
    if (kind === null) kind = k;
    else if (kind !== k) return "other";
  }
  if (kind === null) return "number";
  if (kind === "number") return "number";
  if (kind === "time") return "time";
  if (kind === "boolean") return "logical";
  return "other";
}

// ---- projections ----------------------------------------------------------------------

function crsDef(code) {
  const m = /^EPSG:(\d+)$/i.exec(String(code));
  if (!m) return code; // a proj string
  const epsg = Number(m[1]);
  if (epsg === 4326) return "+proj=longlat +datum=WGS84 +no_defs";
  if (epsg === 3857) return "+proj=merc +a=6378137 +b=6378137 +lat_ts=0 +lon_0=0 +x_0=0 +y_0=0 +k=1 +units=m +nadgrids=@null +no_defs";
  if (epsg >= 32601 && epsg <= 32660) return `+proj=utm +zone=${epsg - 32600} +datum=WGS84 +units=m +no_defs`;
  if (epsg >= 32701 && epsg <= 32760) return `+proj=utm +zone=${epsg - 32700} +south +datum=WGS84 +units=m +no_defs`;
  throw new Error(`CRS ${code} is not supported; use EPSG:4326, EPSG:3857 or a WGS 84 / UTM zone.`);
}

const converters = new Map();
function converter(from, to) {
  const key = `${from}|${to}`;
  if (!converters.has(key)) converters.set(key, proj4(crsDef(from), crsDef(to)));
  return converters.get(key);
}

// project coordinates; returns [xs, ys]; NA stays NA
export function project(a, b, from, to) {
  const conv = converter(from, to);
  const xs = new Array(a.length);
  const ys = new Array(a.length);
  for (let i = 0; i < a.length; i++) {
    if (isNA(a[i]) || isNA(b[i])) {
      xs[i] = NaN;
      ys[i] = NaN;
    } else {
      const p = conv.forward([a[i], b[i]]);
      xs[i] = p[0];
      ys[i] = p[1];
    }
  }
  return [xs, ys];
}

export function utmEpsg(lon, lat) {
  const lo = naOmit(lon);
  const la = naOmit(lat);
  const lon0 = (Math.min(...lo) + Math.max(...lo)) / 2;
  const lat0 = (Math.min(...la) + Math.max(...la)) / 2;
  let zone = Math.floor((lon0 + 180) / 6) + 1;
  zone = Math.min(Math.max(zone, 1), 60);
  return lat0 < 0 ? 32700 + zone : 32600 + zone;
}

export const toMetric = (lon, lat, epsg) => project(lon, lat, "EPSG:4326", `EPSG:${epsg}`);
export const toLonlat = (x, y, epsg) => project(x, y, `EPSG:${epsg}`, "EPSG:4326");

// ---- segments ----------------------------------------------------------------------------

export function segmentId(t) {
  const n = nrow(t);
  if (!("segment" in t)) return new Array(n).fill(1);
  const map = new Map();
  return t.segment.map((v) => {
    const k = isNA(v) ? "NA" : String(v);
    if (!map.has(k)) map.set(k, map.size + 1);
    return map.get(k);
  });
}

export function checkColumns(t, columns, fun) {
  const missing = columns.filter((c) => !(c in t));
  if (missing.length) throw new Error(`${fun}() needs column(s): ${missing.join(", ")}${NAMING_HINT}`);
  return true;
}

export function isFirstOfSegment(seg) {
  return seg.map((s, i) => i === 0 || s !== seg[i - 1]);
}

export function diffBySegment(x, seg, first = 0) {
  return x.map((v, i) => (i === 0 || seg[i] !== seg[i - 1] ? first : v - x[i - 1]));
}

export function uniqueInOrder(v) {
  return [...new Set(v)];
}

export function cumsum(v) {
  let s = 0;
  return v.map((x) => (s += x));
}

export const pathDistance = (x, y) => {
  const out = [0];
  for (let i = 1; i < x.length; i++) out.push(Math.sqrt((x[i] - x[i - 1]) ** 2 + (y[i] - y[i - 1]) ** 2));
  return cumsum(out);
};

// Nearest row of target for every row of query: 1-based index into the complete
// rows of target, and distance. 2D reproduces sf::st_nearest_feature; 3D is a
// brute-force search keeping the first of equidistant rows.
export function nearestNeighbour(query, target) {
  const tgt = target.filter((r) => !r.some(isNA));
  const index = new Array(query.length).fill(-1);
  const dist = new Array(query.length).fill(NaN);
  const dims = query.length ? query[0].length : 2;
  const tree = dims === 2 ? new PointIndex(tgt.map((r) => r[0]), tgt.map((r) => r[1])) : null;
  query.forEach((q, i) => {
    if (q.some(isNA)) return;
    let best;
    if (tree) best = tree.nearest(q[0], q[1]);
    else {
      let bd = Infinity;
      for (let j = 0; j < tgt.length; j++) {
        let d2 = 0;
        for (let k = 0; k < dims; k++) d2 += (q[k] - tgt[j][k]) ** 2;
        if (d2 < bd) { bd = d2; best = j; }
      }
    }
    let d2 = 0;
    for (let k = 0; k < dims; k++) d2 += (q[k] - tgt[best][k]) ** 2;
    index[i] = best + 1;
    dist[i] = Math.sqrt(d2);
  });
  return { index, dist };
}
