// Reading and writing GPX (port of R/read_stream.R and R/write_stream.R).
// The R package reads GPX with GDAL; this port parses the text directly and
// returns what the GDAL GPX driver returns for track and route points. Point
// extensions use the same regular expressions as the R code. Everything works
// on text so it runs in the browser; readStream/writeStream add file access
// in Node.

import { attrs, checkColumns, numCol, segmentId, isNA } from "./utils.js";

const ATTR = /([A-Za-z_:][A-Za-z0-9_.:-]*)\s*=\s*("([^"]*)"|'([^']*)')/g;
const RESERVED = new Set(["lon", "lat", "ele", "timestamp", "time", "segment", "lng", "longitude",
  "latitude", "elevation", "altitude", "alt", "datetime"]);

const local = (tag) => tag.split(":").pop();

function parseAttrs(startTag) {
  const out = {};
  for (const m of startTag.matchAll(ATTR)) out[local(m[1])] = m[3] !== undefined ? m[3] : m[4];
  return out;
}

function childText(block, name) {
  const re = new RegExp(`<(?:[A-Za-z0-9_]+:)?${name}(?:\\s[^>]*)?>([\\s\\S]*?)</(?:[A-Za-z0-9_]+:)?${name}\\s*>`);
  const m = re.exec(block);
  return m ? m[1].trim() : null;
}

function num(s) {
  if (s === null || s === undefined) return NaN;
  const t = String(s).trim();
  if (t === "") return NaN;
  const x = Number(t);
  return Number.isNaN(x) ? NaN : x;
}

// ISO 8601 date-time to seconds since the epoch (UTC); null if invalid
export function parseIsoTime(text) {
  if (text === null || text === undefined || text === "") return null;
  const m = /^(\d{4})-(\d{2})-(\d{2})[T ](\d{2}):(\d{2}):(\d{2})(\.\d+)?(Z|[+-]\d{2}:?\d{2})?$/.exec(text.trim());
  if (!m) return null;
  let secs = Date.UTC(+m[1], +m[2] - 1, +m[3], +m[4], +m[5], +m[6]) / 1000 + (m[7] ? Number(m[7]) : 0);
  const tz = m[8];
  if (tz && tz !== "Z") {
    const sign = tz[0] === "+" ? 1 : -1;
    const d = tz.slice(1).replace(":", "");
    secs -= sign * (Number(d.slice(0, 2)) * 3600 + Number(d.slice(2)) * 60);
  }
  return secs;
}

function points(text, container, segTag, pointTag) {
  const rows = [];
  const contRe = new RegExp(`<(?:[A-Za-z0-9_]+:)?${container}[\\s>][\\s\\S]*?</(?:[A-Za-z0-9_]+:)?${container}\\s*>`, "g");
  const ptRe = new RegExp(`<(?:[A-Za-z0-9_]+:)?${pointTag}(\\s[^>]*?)?(/>|>([\\s\\S]*?)</(?:[A-Za-z0-9_]+:)?${pointTag}\\s*>)`, "g");
  let fid = 0;
  for (const cm of text.matchAll(contRe)) {
    const block = cm[0];
    let groups;
    if (segTag === null) groups = [block];
    else {
      const segRe = new RegExp(`<(?:[A-Za-z0-9_]+:)?${segTag}[\\s>][\\s\\S]*?</(?:[A-Za-z0-9_]+:)?${segTag}\\s*>`, "g");
      groups = [...block.matchAll(segRe)].map((m) => m[0]);
    }
    groups.forEach((g, sid) => {
      for (const pm of g.matchAll(ptRe)) {
        const a = parseAttrs(pm[1] || "");
        const core = (pm[3] || "").replace(/<extensions>[\s\S]*?<\/extensions>/g, "");
        rows.push({ lon: num(a.lon), lat: num(a.lat), ele: num(childText(core, "ele")),
          time: parseIsoTime(childText(core, "time")), key: `${fid} ${sid}` });
      }
    });
    fid++;
  }
  return rows;
}

// read_stream() on the text of a GPX file
export function readStreamText(raw, { type = "auto", name = "input" } = {}) {
  if (!["auto", "track", "route"].includes(type)) throw new Error("'type' must be one of 'auto', 'track', 'route'");
  const text = raw.split(/\r\n|\r|\n/).join("\n");
  const clean = text.replace(/<!--[\s\S]*?-->/g, "");
  const track = points(clean, "trk", "trkseg", "trkpt");
  const route = points(clean, "rte", null, "rtept");
  const layer = type === "auto" ? (track.length > 0 ? "track" : "route") : type;
  const rows = layer === "track" ? track : route;
  if (rows.length === 0) throw new Error(`No ${layer} points found in ${name} (waypoint-only files are not streams).`);

  const stream = { lon: rows.map((r) => r.lon), lat: rows.map((r) => r.lat) };
  if (rows.some((r) => !Number.isNaN(r.ele))) stream.ele = rows.map((r) => r.ele);
  if (rows.some((r) => r.time !== null)) stream.timestamp = rows.map((r) => (r.time === null ? null : new Date(Math.round(r.time * 1000))));
  const keys = new Map();
  stream.segment = rows.map((r) => {
    if (!keys.has(r.key)) keys.set(r.key, keys.size + 1);
    return keys.get(r.key);
  });

  if (text.includes("<extensions>")) {
    const tag = layer === "track" ? "trkpt" : "rtept";
    const blocks = [...text.matchAll(new RegExp(`<${tag}[ >][\\s\\S]*?(?:</${tag}>|/>)`, "g"))].map((m) => m[0]);
    if (blocks.length === rows.length) {
      const ext = blocks.map((b) => (b.includes("<extensions>") ? b.replace(/^[\s\S]*?<extensions>([\s\S]*)<\/extensions>[\s\S]*$/, "$1") : ""));
      const leafRe = /<([A-Za-z0-9_]+:)?[A-Za-z0-9_]+>[^<]*<\//g;
      const leafNames = [];
      for (const e of ext) {
        for (const m of e.matchAll(leafRe)) {
          const leaf = m[0].replace(/^<([A-Za-z0-9_]+:)?([A-Za-z0-9_]+)>[\s\S]*$/, "$2");
          if (!leafNames.includes(leaf)) leafNames.push(leaf);
        }
      }
      for (const leaf of leafNames) {
        const re = new RegExp(`<([A-Za-z0-9_]+:)?${leaf}>([^<]*)</`);
        const value = ext.map((e) => {
          const m = re.exec(e);
          return m ? num(m[2]) : NaN;
        });
        let col = leaf;
        if (RESERVED.has(col.toLowerCase()) || col in stream) col = `ext_${leaf}`;
        if (value.some((v) => !Number.isNaN(v)) && !(col in stream)) stream[col] = value;
      }
    } else {
      console.warn("Point extensions could not be aligned with the points and were skipped.");
    }
  }
  return stream;
}

// read a GPX file (Node only)
export async function readStream(file, { type = "auto" } = {}) {
  const fs = await import("node:fs");
  const path = await import("node:path");
  if (!fs.existsSync(file)) throw new Error(`File not found: ${file}`);
  if (!/\.gpx$/i.test(file)) throw new Error("read_stream() only reads GPX files.");
  return readStreamText(fs.readFileSync(file, "utf8"), { type, name: path.basename(file) });
}

// C printf("%.*f"): round the exact binary value, ties to even
export function formatFixed(x, digits) {
  if (Number.isNaN(x)) return "NA";
  const neg = x < 0 || Object.is(x, -0);
  const exact = Math.abs(x).toFixed(100);
  const [ip, fp] = exact.split(".");
  const keep = fp.slice(0, digits);
  const rest = fp.slice(digits);
  let digitsStr = ip + keep;
  const firstRest = rest.charCodeAt(0) - 48;
  const tail = rest.slice(1).replace(/0+$/, "");
  const last = digitsStr.charCodeAt(digitsStr.length - 1) - 48;
  const roundUp = firstRest > 5 || (firstRest === 5 && (tail.length > 0 || last % 2 === 1));
  if (roundUp) {
    const arr = digitsStr.split("").map(Number);
    let i = arr.length - 1;
    while (i >= 0) {
      if (arr[i] === 9) { arr[i] = 0; i--; } else { arr[i]++; break; }
    }
    if (i < 0) arr.unshift(1);
    digitsStr = arr.join("");
  }
  const intPart = digitsStr.slice(0, digitsStr.length - digits) || "0";
  const out = digits > 0 ? `${intPart}.${digitsStr.slice(digitsStr.length - digits)}` : intPart;
  const isZero = /^[0.]+$/.test(out);
  return (neg && !(isZero && x === 0) ? "-" : "") + out;
}

function formatTime(d) {
  const p = (v, n = 2) => String(v).padStart(n, "0");
  const t = new Date(Math.floor(d.getTime() / 1000) * 1000);
  return `${p(t.getUTCFullYear(), 4)}-${p(t.getUTCMonth() + 1)}-${p(t.getUTCDate())}T${p(t.getUTCHours())}:${p(t.getUTCMinutes())}:${p(t.getUTCSeconds())}Z`;
}

// GPX text that write_stream() writes
export function writeStreamText(data, { name = null } = {}) {
  checkColumns(data, ["lon", "lat"], "write_stream");
  const seg = segmentId(data);
  const lat = numCol(data.lat);
  const lon = numCol(data.lon);
  const hasEle = "ele" in data;
  const hasTime = "timestamp" in data;
  const ele = hasEle ? numCol(data.ele) : null;
  const pts = lat.map((_, i) => {
    let p = `<trkpt lat="${formatFixed(lat[i], 8)}" lon="${formatFixed(lon[i], 8)}">`;
    if (hasEle && !Number.isNaN(ele[i])) p += `<ele>${formatFixed(ele[i], 2)}</ele>`;
    if (hasTime && !isNA(data.timestamp[i])) p += `<time>${formatTime(data.timestamp[i])}</time>`;
    return `      ${p}</trkpt>`;
  });
  const out = ['<?xml version="1.0" encoding="UTF-8"?>',
    '<gpx version="1.1" creator="metrictrails" xmlns="http://www.topografix.com/GPX/1/1">', "  <trk>"];
  if (name !== null) out.push(`    <name>${name}</name>`);
  for (const s of [...new Set(seg)]) {
    out.push("    <trkseg>");
    pts.forEach((p, i) => { if (seg[i] === s) out.push(p); });
    out.push("    </trkseg>");
  }
  out.push("  </trk>", "</gpx>");
  return `${out.join("\n")}\n`;
}

export async function writeStream(data, file, { name = null } = {}) {
  const fs = await import("node:fs");
  fs.writeFileSync(file, writeStreamText(data, { name }), "utf8");
  return file;
}

export { attrs };
