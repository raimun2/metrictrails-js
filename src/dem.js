// Digital Elevation Models (port of R/dem.R and dem.py).
//
// A DEM is a plain Raster object: { values (row-major Float64Array, row 0 at
// the top), nrow, ncol, xmin, ymax, xres, yres, crs: "EPSG:xxxx" }. The
// extraction, refinement and aggregation arithmetic is written out so that it
// matches terra. GeoTIFF files are read with geotiff.js.
//
// get_dem() differs from R by design: R reprojects the Terrain Tiles mosaic to
// UTM with GDAL (a second bilinear interpolation) before extraction. Without
// GDAL, this port keeps the mosaic in Web Mercator (EPSG:3857), the tiles'
// own grid, and ele_correction() projects the readings to it. Elevations are
// therefore close to, but not identical with, R's.

import * as rs from "./rstats.js";
import { asTable, checkColumns, isNA, numCol, pathDistance, project, segmentId, toLonlat, toMetric,
  uniqueInOrder, utmEpsg } from "./utils.js";

export const TILE_URL = "https://s3.amazonaws.com/elevation-tiles-prod/geotiff/{z}/{x}/{y}.tif";

export function makeRaster(values, nrow, ncol, xmin, ymax, xres, yres, crs) {
  return { values: Float64Array.from(values), nrow, ncol, xmin, ymax, xres, yres, crs,
    get xmax() { return this.xmin + this.ncol * this.xres; },
    get ymin() { return this.ymax - this.nrow * this.yres; } };
}

// Read the first band of a GeoTIFF from an ArrayBuffer (or a path in Node)
export async function read_raster(source) {
  const GeoTIFF = await import("geotiff");
  let buffer = source;
  if (typeof source === "string") {
    const fs = await import("node:fs");
    const b = fs.readFileSync(source);
    buffer = b.buffer.slice(b.byteOffset, b.byteOffset + b.byteLength);
  }
  const tiff = await GeoTIFF.fromArrayBuffer(buffer);
  const image = await tiff.getImage();
  const [band] = await image.readRasters();
  const [xmin, , , ymax] = image.getBoundingBox();
  const [xres, yres] = image.getResolution();
  const nodata = image.getGDALNoData();
  const keys = image.getGeoKeys() || {};
  const epsg = keys.ProjectedCSTypeGeoKey || keys.GeographicTypeGeoKey || 4326;
  const values = Float64Array.from(band, (v) => (nodata !== null && v === nodata ? NaN : v));
  return makeRaster(values, image.getHeight(), image.getWidth(), xmin, ymax, Math.abs(xres), Math.abs(yres), `EPSG:${epsg}`);
}

const isLonlat = (dem) => dem.crs === "EPSG:4326";

// ---- extraction ---------------------------------------------------------------------

export function extractSimple(dem, x, y) {
  return x.map((xi, k) => {
    const yi = y[k];
    let col = Math.floor((xi - dem.xmin) / dem.xres);
    let row = Math.floor((dem.ymax - yi) / dem.yres);
    if (xi === dem.xmax) col = dem.ncol - 1;
    if (yi === dem.ymin) row = dem.nrow - 1;
    if (!(col >= 0 && col < dem.ncol && row >= 0 && row < dem.nrow)) return NaN;
    return dem.values[row * dem.ncol + col];
  });
}

const lerp = (a, b, w) => (Number.isNaN(a) ? b : Number.isNaN(b) ? a : a + w * (b - a));

// Bilinear interpolation between the four nearest cell centres (terra 'bilinear')
export function extractBilinear(dem, x, y) {
  const { ncol, nrow, values: v } = dem;
  return x.map((xi, k) => {
    const yi = y[k];
    if (!(xi >= dem.xmin && xi <= dem.xmax && yi >= dem.ymin && yi <= dem.ymax)) return NaN;
    const cf = Math.min(Math.max((xi - dem.xmin) / dem.xres - 0.5, 0), ncol - 1);
    const rf = Math.min(Math.max((dem.ymax - yi) / dem.yres - 0.5, 0), nrow - 1);
    const c0 = Math.min(Math.floor(cf), Math.max(ncol - 2, 0));
    const r0 = Math.min(Math.floor(rf), Math.max(nrow - 2, 0));
    const c1 = Math.min(c0 + 1, ncol - 1);
    const r1 = Math.min(r0 + 1, nrow - 1);
    const top = lerp(v[r0 * ncol + c0], v[r0 * ncol + c1], cf - c0);
    const bottom = lerp(v[r1 * ncol + c0], v[r1 * ncol + c1], cf - c0);
    return lerp(top, bottom, rf - r0);
  });
}

// ---- resampling ------------------------------------------------------------------------

// GDAL-style warp (bilinear or nearest) of a raster to a grid in the same CRS
function warpSameCrs(dem, xmin, ymax, res, nrow, ncol, method) {
  const out = new Float64Array(nrow * ncol).fill(NaN);
  const v = dem.values;
  const xscale = Math.min(dem.xres / res, 1);
  const yscale = Math.min(dem.yres / res, 1);
  const rx = 1 / xscale;
  const ry = 1 / yscale;
  const nx = Math.ceil(2 * rx) + 1;
  const ny = Math.ceil(2 * ry) + 1;
  for (let r = 0; r < nrow; r++) {
    const py = (dem.ymax - (ymax - (r + 0.5) * res)) / dem.yres;
    for (let c = 0; c < ncol; c++) {
      const px = (xmin + (c + 0.5) * res - dem.xmin) / dem.xres;
      if (!(px >= 0 && px < dem.ncol && py >= 0 && py < dem.nrow)) continue;
      if (method === "near") {
        const ix = Math.floor(px + 1e-10);
        const iy = Math.floor(py + 1e-10);
        if (ix >= 0 && ix < dem.ncol && iy >= 0 && iy < dem.nrow) out[r * ncol + c] = v[iy * dem.ncol + ix];
        continue;
      }
      const x0 = Math.floor(px - 0.5 - rx + 1);
      const y0 = Math.floor(py - 0.5 - ry + 1);
      let acc = 0;
      let wsum = 0;
      for (let dy = 0; dy < ny; dy++) {
        const iy = y0 + dy;
        const wy = 1 - Math.abs((iy + 0.5 - py) * yscale);
        if (wy <= 0 || iy < 0 || iy >= dem.nrow) continue;
        for (let dx = 0; dx < nx; dx++) {
          const ix = x0 + dx;
          const wx = 1 - Math.abs((ix + 0.5 - px) * xscale);
          if (wx <= 0 || ix < 0 || ix >= dem.ncol) continue;
          const val = v[iy * dem.ncol + ix];
          if (Number.isNaN(val)) continue;
          acc += wx * wy * val;
          wsum += wx * wy;
        }
      }
      if (wsum >= 1e-6) out[r * ncol + c] = acc / wsum;
    }
  }
  return out;
}

function aggregate(dem, fact, fun) {
  const nr = Math.ceil(dem.nrow / fact);
  const nc = Math.ceil(dem.ncol / fact);
  const f = { mean: rs.rMean, median: rs.rMedian, min: (a) => Math.min(...a), max: (a) => Math.max(...a) }[fun];
  const out = new Float64Array(nr * nc).fill(NaN);
  for (let i = 0; i < nr; i++) {
    for (let j = 0; j < nc; j++) {
      const block = [];
      for (let r = i * fact; r < Math.min((i + 1) * fact, dem.nrow); r++) {
        for (let c = j * fact; c < Math.min((j + 1) * fact, dem.ncol); c++) {
          const val = dem.values[r * dem.ncol + c];
          if (!Number.isNaN(val)) block.push(val);
        }
      }
      if (block.length) out[i * nc + j] = f(block);
    }
  }
  return makeRaster(out, nr, nc, dem.xmin, dem.ymax, dem.xres * fact, dem.yres * fact, dem.crs);
}

// Refine (bilinear, near) or coarsen (mean, median, min, max) a DEM in a metric CRS
export function resample_dem(dem, res, { method = "bilinear" } = {}) {
  if (!["bilinear", "near", "mean", "median", "min", "max"].includes(method)) {
    throw new Error("'method' must be one of bilinear, near, mean, median, min, max");
  }
  if (isLonlat(dem)) throw new Error("The DEM must be in a metric CRS; project it first.");
  if (method === "bilinear" || method === "near") {
    const nc = Math.ceil((dem.xmax - dem.xmin) / res);
    const nr = Math.ceil((dem.ymax - dem.ymin) / res);
    return makeRaster(warpSameCrs(dem, dem.xmin, dem.ymax, res, nr, nc, method), nr, nc, dem.xmin, dem.ymax, res, res, dem.crs);
  }
  const fact = Math.round(res / dem.xres);
  if (fact < 2) {
    throw new Error(`Aggregation needs 'res' to be at least twice the current resolution (${Number(dem.xres.toPrecision(3))} m).`);
  }
  return aggregate(dem, fact, method);
}

// ---- download -----------------------------------------------------------------------------

const tileX = (lon, z) => Math.floor(((lon + 180) / 360) * 2 ** z);
const tileY = (lat, z) => {
  const r = (lat * Math.PI) / 180;
  return Math.floor(((1 - Math.log(Math.tan(r) + 1 / Math.cos(r)) / Math.PI) / 2) * 2 ** z);
};

// Slippy-map tiles that intersect a corridor of `buffer` metres around the path
export function dem_tiles(data, { z = 14, buffer = 500 } = {}) {
  const lon = numCol(data.lon);
  const lat = numCol(data.lat);
  const epsg = utmEpsg(lon, lat);
  const [x0, y0] = toMetric(lon, lat, epsg);
  const seg0 = segmentId(data);
  const keep = x0.map((v, i) => !isNA(v) && !isNA(y0[i]));
  const x = x0.filter((_, i) => keep[i]);
  const y = y0.filter((_, i) => keep[i]);
  const seg = seg0.filter((_, i) => keep[i]);
  const maxLat = Math.max(...rs.naOmit(lat).map(Math.abs));
  const tileM = (40075016.7 * Math.cos((maxLat * Math.PI) / 180)) / 2 ** z;
  const spacing = Math.min(buffer, tileM / 2);
  const px = [];
  const py = [];
  for (const s of uniqueInOrder(seg)) {
    const rows = [];
    seg.forEach((v, i) => { if (v === s) rows.push(i); });
    const xs = rows.map((i) => x[i]);
    const ys = rows.map((i) => y[i]);
    if (rows.length === 1) { px.push(xs[0]); py.push(ys[0]); continue; }
    const along = pathDistance(xs, ys);
    const total = Math.max(...along);
    const sNew = [...new Set([...rs.seqBy(0, total, Math.max(spacing, 1)), total])];
    px.push(...rs.approx(along, xs, sNew, { ties: "mean" }));
    py.push(...rs.approx(along, ys, sNew, { ties: "mean" }));
  }
  const offsets = rs.seqLength(-buffer, buffer, Math.max(2, Math.ceil((2 * buffer) / (tileM / 2)) + 1));
  const tiles = [];
  const seen = new Set();
  for (const dy of offsets) {
    for (const dx of offsets) {
      const [lo, la] = toLonlat(px.map((v) => v + dx), py.map((v) => v + dy), epsg);
      lo.forEach((l, i) => {
        const key = `${tileX(l, z)},${tileY(la[i], z)}`;
        if (!seen.has(key)) { seen.add(key); tiles.push([tileX(l, z), tileY(la[i], z)]); }
      });
    }
  }
  return { tiles, epsg, x, y };
}

// Download Terrain Tiles around the stream and mosaic them (EPSG:3857, see above)
export async function get_dem(data, { z = 14, buffer = 500, max_tiles = Infinity, fetchFn = globalThis.fetch } = {}) {
  checkColumns(data, ["lon", "lat"], "get_dem");
  if (z < 1 || z > 14) throw new Error("'z' must be between 1 and 14.");
  const { tiles } = dem_tiles(data, { z, buffer });
  if (tiles.length > max_tiles) throw new Error(`${tiles.length} tiles needed, more than 'max_tiles'.`);
  const rasters = [];
  for (const [tx, ty] of tiles) {
    const url = TILE_URL.replace("{z}", z).replace("{x}", tx).replace("{y}", ty);
    let raster = null;
    for (let attempt = 0; attempt < 3 && raster === null; attempt++) {
      try {
        const resp = await fetchFn(url);
        if (resp.ok) raster = await read_raster(await resp.arrayBuffer());
      } catch (e) {
        raster = null;
      }
    }
    if (raster === null) throw new Error(`Could not download a valid tile from ${url}`);
    rasters.push(raster);
  }
  // mosaic on the common Web Mercator grid of the tiles
  const res = rasters[0].xres;
  const xmin = Math.min(...rasters.map((r) => r.xmin));
  const ymax = Math.max(...rasters.map((r) => r.ymax));
  const xmax = Math.max(...rasters.map((r) => r.xmax));
  const ymin = Math.min(...rasters.map((r) => r.ymin));
  const ncol = Math.round((xmax - xmin) / res);
  const nrow = Math.round((ymax - ymin) / res);
  const values = new Float64Array(nrow * ncol).fill(NaN);
  for (const r of rasters) {
    const c0 = Math.round((r.xmin - xmin) / res);
    const r0 = Math.round((ymax - r.ymax) / res);
    for (let i = 0; i < r.nrow; i++) {
      for (let j = 0; j < r.ncol; j++) values[(r0 + i) * ncol + c0 + j] = r.values[i * r.ncol + j];
    }
  }
  return makeRaster(values, nrow, ncol, xmin, ymax, res, res, "EPSG:3857");
}

// Replace ele with DEM values (or add them as ele_dem)
export async function ele_correction(data, { dem = null, z = 14, extract = "bilinear", replace = true } = {}) {
  if (!["bilinear", "simple"].includes(extract)) throw new Error("'extract' must be 'bilinear' or 'simple'");
  checkColumns(data, ["lon", "lat"], "ele_correction");
  if (dem === null) dem = await get_dem(data, { z });
  else if (typeof dem === "string" || dem instanceof ArrayBuffer) dem = await read_raster(dem);
  const [px, py] = project(numCol(data.lon), numCol(data.lat), "EPSG:4326", dem.crs);
  const ele = (extract === "bilinear" ? extractBilinear : extractSimple)(dem, px, py);
  const nNA = ele.filter((v) => Number.isNaN(v)).length;
  if (nNA) console.warn(`${nNA} readings fall outside the DEM and get NA elevation.`);
  const t = asTable(data);
  if (replace) t.ele = ele;
  else t.ele_dem = ele;
  return t;
}
