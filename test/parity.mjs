// Parity of the JavaScript port with the R package: every golden case
// (parity/make_golden.R, run on the CRAN release) is rerun here and compared.

import fs from "node:fs";
import path from "node:path";
import { fileURLToPath } from "node:url";
import * as mt from "../src/index.js";
import { attrs } from "../src/utils.js";

const ATOL = 1e-6;
const RTOL = 1e-9;
const here = path.dirname(fileURLToPath(import.meta.url));
const extdata = (f) => path.join(here, "..", "extdata", f);

function parseCsv(text) {
  const rows = [];
  let row = [];
  let field = "";
  let quoted = false;
  let wasQuoted = false;
  for (let i = 0; i < text.length; i++) {
    const ch = text[i];
    if (quoted) {
      if (ch === '"') {
        if (text[i + 1] === '"') { field += '"'; i++; } else quoted = false;
      } else field += ch;
    } else if (ch === '"') { quoted = true; wasQuoted = true; }
    else if (ch === ",") { row.push(wasQuoted ? field : field === "NA" ? null : field); field = ""; wasQuoted = false; }
    else if (ch === "\n" || ch === "\r") {
      if (ch === "\r" && text[i + 1] === "\n") i++;
      row.push(wasQuoted ? field : field === "NA" ? null : field);
      rows.push(row);
      row = []; field = ""; wasQuoted = false;
    } else field += ch;
  }
  if (field !== "" || row.length) { row.push(wasQuoted ? field : field === "NA" ? null : field); rows.push(row); }
  return rows;
}

function readTable(file) {
  const types = JSON.parse(fs.readFileSync(file.replace(/\.csv$/, ".types.json"), "utf8"));
  const rows = parseCsv(fs.readFileSync(file, "utf8"));
  const header = rows[0];
  const out = {};
  header.forEach((h, j) => {
    const kind = types[h] || "double";
    out[h] = rows.slice(1).map((r) => {
      const v = r[j];
      if (kind === "double" || kind === "integer") return v === null ? NaN : Number(v);
      if (kind === "timestamp") return v === null ? null : new Date(Math.round(Number(v) * 1000));
      if (kind === "logical") return v === null ? null : v === "TRUE";
      return v;
    });
  });
  return { table: out, types };
}

const asNum = (v) => (v === null || v === undefined ? NaN : v instanceof Date ? v.getTime() / 1000
  : typeof v === "boolean" ? Number(v) : Number(v));

function compareValues(name, got, want, kind, atol = ATOL) {
  if (got.length !== want.length) return { ok: false, worst: Infinity, msg: `${name}: length ${got.length} vs ${want.length}` };
  if (kind === "character" || kind === "logical") {
    for (let i = 0; i < want.length; i++) {
      const g = got[i] === undefined ? null : got[i];
      const w = want[i];
      const gn = g === null || (typeof g === "number" && Number.isNaN(g)) ? null : kind === "logical" ? Boolean(g) : String(g);
      if (gn !== w) return { ok: false, worst: 0, msg: `${name}: mismatch at ${i}: ${g} vs ${w}` };
    }
    return { ok: true, worst: 0 };
  }
  let worst = 0;
  if (kind === "timestamp") atol = Math.max(atol, 1e-3); // Date has millisecond resolution
  for (let i = 0; i < want.length; i++) {
    const g = asNum(got[i]);
    const w = asNum(want[i]);
    if (Number.isNaN(g) !== Number.isNaN(w)) return { ok: false, worst: Infinity, msg: `${name}: NA pattern differs at ${i}: ${g} vs ${w}` };
    if (Number.isNaN(w)) continue;
    if (!Number.isFinite(w) || !Number.isFinite(g)) {
      if (g !== w) return { ok: false, worst: Infinity, msg: `${name}: ${g} vs ${w} at ${i}` };
      continue;
    }
    const e = Math.abs(g - w);
    worst = Math.max(worst, e);
    if (e > atol + RTOL * Math.abs(w)) return { ok: false, worst: e, msg: `${name}: error ${e.toPrecision(3)} at ${i} (got ${g}, want ${w})` };
  }
  return { ok: true, worst };
}

function compareTables(got, want, types) {
  const msgs = [];
  let worst = 0;
  const gc = Object.keys(got);
  const wc = Object.keys(want);
  if (gc.join("|") !== wc.join("|")) msgs.push(`columns differ:\n   got  ${gc}\n   want ${wc}`);
  for (const c of wc) {
    if (!(c in got)) continue;
    const r = compareValues(c, Array.from(got[c]), want[c], types[c] || "double");
    worst = Math.max(worst, r.worst);
    if (!r.ok) msgs.push(r.msg);
  }
  return { worst, msgs };
}

function compareJson(got, want, p = "attrs") {
  if (want === null) return got === null || got === undefined || (typeof got === "number" && Number.isNaN(got)) ? [] : [`${p}: ${got} vs NA`];
  if (Array.isArray(want)) {
    let g = got;
    if (g && !Array.isArray(g) && typeof g === "object") g = Object.values(g);
    if (!Array.isArray(g) || g.length !== want.length) return [`${p}: length differs`];
    return want.flatMap((w, i) => compareJson(g[i], w, `${p}[${i}]`));
  }
  if (typeof want === "object") {
    return Object.entries(want).flatMap(([k, w]) => (got && k in got ? compareJson(got[k], w, `${p}.${k}`) : [`${p}.${k} missing`]));
  }
  if (typeof want === "string" || typeof want === "boolean") return got === want ? [] : [`${p}: ${got} vs ${want}`];
  const g = Number(got);
  return Number.isNaN(g) || Math.abs(g - want) > ATOL + RTOL * Math.abs(want) ? [`${p}: ${got} vs ${want}`] : [];
}

async function runCase(name, spec, golden) {
  const d = path.join(golden, name);
  const inputs = {};
  const spec_inputs = spec.inputs && !Array.isArray(spec.inputs) ? spec.inputs : {};
  for (const [nm, rel] of Object.entries(spec_inputs)) inputs[nm] = readTable(path.join(golden, rel)).table;
  const args = spec.args && !Array.isArray(spec.args) ? { ...spec.args } : {};
  const fun = spec.fun;
  const origWarn = console.warn;
  console.warn = () => {};
  let result;
  let extra = {};
  try {
    if (fun === "read_stream") {
      result = mt.read_stream_text(fs.readFileSync(extdata(args.file), "utf8"));
    } else if (fun === "write_stream") {
      const text = mt.write_stream_text(inputs.data, args);
      const want = fs.readFileSync(path.join(d, "output.gpx"), "utf8").replace(/\r\n/g, "\n");
      return { ok: text === want, worst: 0, msgs: text === want ? [] : ["GPX text differs"] };
    } else if (fun === "resample_dem") {
      const dem = await mt.read_raster(extdata(args.dem));
      const r = mt.resample_dem(dem, args.res, { method: args.method });
      const want = JSON.parse(fs.readFileSync(path.join(d, "output.json"), "utf8"));
      const msgs = [];
      for (const k of ["xmin", "xmax", "ymin", "ymax", "nrow", "ncol"]) if (Math.abs(r[k] - want[k]) > 1e-6) msgs.push(`${k}: ${r[k]} vs ${want[k]}`);
      if (msgs.length) return { ok: false, worst: Infinity, msgs };
      const c = compareValues("values", Array.from(r.values), want.values.map((v) => (v === null ? NaN : v)), "double",
        args.method === "bilinear" ? 1e-3 : ATOL);
      return { ok: c.ok, worst: c.worst, msgs: c.ok ? [] : [c.msg] };
    } else if (fun === "ele_correction") {
      const dem = await mt.read_raster(extdata(args.dem));
      result = await mt.ele_correction(inputs.data, { ...args, dem });
    } else if (fun === "tobler_speed") {
      result = { speed: mt.tobler_speed(inputs.data.slope, args.params) };
    } else if (fun === "compare_stream") result = mt.compare_stream(inputs.a, inputs.b, args);
    else if (fun === "detect_stream") result = mt.detect_stream(inputs.data, inputs.route, args);
    else if (fun === "locate_stream") result = mt.locate_stream(inputs.data, inputs.reference);
    else if (fun === "match_stream") result = mt.match_stream(inputs.data, inputs.reference, inputs.params, args);
    else if (fun === "fractal_stream") {
      const fd = mt.fractal_stream(inputs.data, args);
      result = fd.scaling;
      extra = { dimension: fd.dimension, r_squared: fd.r_squared };
    } else result = mt[fun](inputs.data, args);
  } finally {
    console.warn = origWarn;
  }
  const { table: want, types } = readTable(path.join(d, "output.csv"));
  const { worst, msgs } = compareTables(result, want, types);
  const attrsWant = JSON.parse(fs.readFileSync(path.join(d, "attrs.json"), "utf8"));
  if (attrsWant && Object.keys(attrsWant).length) {
    let got = { ...(result.attrs || {}), ...extra };
    if (fun === "sanitize_stream") got = got.sanitize_report;
    if (fun === "match_stream") got = { match_report: got.match_report };
    msgs.push(...compareJson(got, attrsWant));
  }
  return { ok: msgs.length === 0, worst, msgs };
}

const golden = process.argv[2] || path.join(here, "..", "..", "parity", "golden");
const only = process.argv.slice(3);
const manifest = JSON.parse(fs.readFileSync(path.join(golden, "manifest.json"), "utf8"));
let failures = 0;
let count = 0;
for (const [name, spec] of Object.entries(manifest.cases)) {
  if (only.length && !only.some((o) => name.includes(o))) continue;
  count++;
  let r;
  try {
    r = await runCase(name, spec, golden);
  } catch (e) {
    r = { ok: false, worst: Infinity, msgs: [`${e.name}: ${e.message}\n${e.stack.split("\n").slice(1, 3).join("\n")}`] };
  }
  if (!r.ok) failures++;
  console.log(`${r.ok ? "PASS" : "FAIL"}  ${name.padEnd(28)} max|err|=${r.worst.toExponential(2)}`);
  for (const m of r.msgs.slice(0, 5)) console.log(`       ${m}`);
}
console.log(`\n${count} cases, ${failures} failing`);
process.exit(failures ? 1 : 0);
