// metrictrails: processing of GNSS streams from sports watches.
// JavaScript port of the R package metrictrails (CRAN), the source of truth.
// Function and argument names follow the R package; optional arguments go in
// a final options object, e.g. sanitize_stream(data, { method: "loess" }).

export const VERSION = "1.0.0";
export const R_VERSION = "1.0.0";

export { readStream as read_stream, readStreamText as read_stream_text, writeStream as write_stream,
  writeStreamText as write_stream_text } from "./io.js";
export { uniform_stream } from "./uniform.js";
export { sanitize_stream } from "./sanitize.js";
export { resolution_stream, resample_stream } from "./resample.js";
export { differential_stream, summarise_stream } from "./measures.js";
export { fractal_stream } from "./fractal.js";
export { tobler_stream, tobler_speed } from "./tobler.js";
export { compare_stream, detect_stream, project_stream, locate_stream, match_stream } from "./twoStreams.js";
export { read_raster, resample_dem, get_dem, ele_correction, makeRaster } from "./dem.js";
export { attrs, toRows, asTable } from "./utils.js";
