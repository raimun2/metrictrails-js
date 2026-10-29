// Naming convention (port of R/uniform_stream.R).

import { asTable, checkColumns, columnKind, isNA, names, numCol, reorder } from "./utils.js";

const DICTIONARY = {
  lon: ["lon", "lng", "long", "longitude", "position_long", "position_lon"],
  lat: ["lat", "latitude", "position_lat"],
  ele: ["ele", "elevation", "altitude", "alt", "enhanced_altitude", "elev"],
  timestamp: ["timestamp", "datetime", "date_time"],
};

// formats tried by as.POSIXct() on text, in order; the first that parses the
// first non-missing value is used for all values; trailing text is ignored
const FORMATS = [
  [/^(\d{4})-(\d{1,2})-(\d{1,2}) (\d{1,2}):(\d{1,2}):(\d{1,2}(?:\.\d*)?)/, 6],
  [/^(\d{4})\/(\d{1,2})\/(\d{1,2}) (\d{1,2}):(\d{1,2}):(\d{1,2}(?:\.\d*)?)/, 6],
  [/^(\d{4})-(\d{1,2})-(\d{1,2}) (\d{1,2}):(\d{1,2})/, 5],
  [/^(\d{4})\/(\d{1,2})\/(\d{1,2}) (\d{1,2}):(\d{1,2})/, 5],
  [/^(\d{4})-(\d{1,2})-(\d{1,2})/, 3],
  [/^(\d{4})\/(\d{1,2})\/(\d{1,2})/, 3],
];

function parseWith([re, parts], s) {
  const m = re.exec(s);
  if (!m) return null;
  const y = +m[1], mo = +m[2], d = +m[3];
  const h = parts >= 5 ? +m[4] : 0;
  const mi = parts >= 5 ? +m[5] : 0;
  const sec = parts === 6 ? Number(m[6]) : 0;
  if (mo < 1 || mo > 12 || d < 1 || d > 31 || h > 24 || mi > 59) return null;
  return Date.UTC(y, mo - 1, d, h, mi, 0) / 1000 + sec;
}

function textToSeconds(values) {
  const texts = values.map((v) => (isNA(v) ? null : String(v).replace(/Z$/, "").replaceAll("T", " ")));
  const first = texts.find((t) => t !== null);
  if (first === undefined) return texts.map(() => NaN);
  const fmt = FORMATS.find((f) => parseWith(f, first) !== null);
  if (!fmt) throw new Error("character string is not in a standard unambiguous format");
  return texts.map((t) => {
    if (t === null) return NaN;
    const s = parseWith(fmt, t);
    return s === null ? NaN : s;
  });
}

function toSeconds(v) {
  if (columnKind(v) === "time") return v.map((d) => (isNA(d) ? NaN : d.getTime() / 1000));
  return textToSeconds(v);
}

// Rename the columns of a stream to the metrictrails naming convention
export function uniform_stream(data) {
  let t = asTable(data);
  let cols = names(t);
  const rename = (from, to) => {
    const out = {};
    for (const c of names(t)) out[c === from ? to : c] = t[c];
    t = out;
    cols = names(t);
  };
  for (const [target, options] of Object.entries(DICTIONARY)) {
    const hits = cols.filter((c) => options.includes(c.toLowerCase()));
    if (hits.length > 1) {
      throw new Error(`Several columns could be '${target}': ${hits.join(", ")}. Remove or rename all but one.`);
    }
    if (hits.length === 1 && hits[0] !== target) rename(hits[0], target);
  }
  const hit = cols.filter((c) => c.toLowerCase() === "time");
  if (hit.length === 1) {
    if (hit[0] !== "time") rename(hit[0], "time");
    const kind = columnKind(t.time);
    if (kind === "time" || kind === "other") {
      if ("timestamp" in t) throw new Error("Both 'timestamp' and a date-time 'time' column are present.");
      rename("time", "timestamp");
    }
  }

  checkColumns(t, ["lon", "lat"], "uniform_stream");
  t.lon = numCol(t.lon);
  t.lat = numCol(t.lat);
  if ("ele" in t) t.ele = numCol(t.ele);

  if ("timestamp" in t) {
    const secs = toSeconds(t.timestamp);
    const ok = secs.filter((s) => !Number.isNaN(s));
    const min = ok.length ? Math.min(...ok) : NaN;
    t.timestamp = secs.map((s) => (Number.isNaN(s) ? null : new Date(Math.round(s * 1000))));
    t.time = secs.map((s) => s - min);
  } else if ("time" in t) {
    t.time = numCol(t.time);
  }

  if ("time" in t && t.time.length > 2) {
    let inc = 0;
    let tot = 0;
    for (let i = 1; i < t.time.length; i++) {
      const d = t.time[i] - t.time[i - 1];
      if (!Number.isNaN(d)) { tot++; if (d > 0) inc++; }
    }
    if (tot === 0 || inc / tot < 0.5) {
      console.warn("Time increases in less than half of the readings; it is removed and the stream is treated as a route without time.");
      delete t.time;
      delete t.timestamp;
    }
  }
  return reorder(t, ["lon", "lat", "ele", "timestamp", "time"]);
}
