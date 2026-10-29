# metrictrails (JavaScript)

JavaScript port of the R package [metrictrails](https://github.com/raimun2/metrictrails) (CRAN), from the [Metric Trails](https://metrictrails.com) lab: tools to read, clean and standardise GNSS streams recorded by sports watches, with a focus on trail and mountain running. It runs in Node and in the browser (everything works on GPX text; file access is optional).

The R package is the source of truth. This port follows its CRAN releases (`R_VERSION`) and is tested against R on every function: the parity workflow runs the CRAN release on golden cases and requires this port to reproduce them (57 cases; differences are at floating-point level, see below).

```js
import * as mt from "metrictrails";

let stream = mt.uniform_stream(mt.read_stream_text(gpxText));
stream = mt.sanitize_stream(stream, { method: "kalman" });
stream = mt.resample_stream(stream);
mt.attrs(stream).res;           // resolution used, as attr(x, "res") in R
mt.summarise_stream(stream);
```

Function and argument names are those of the R package; optional arguments go in a final options object. A stream is a plain object of equal-length column arrays (numbers with `NaN` for NA, `Date` for `timestamp`); what R stores with `attr()` is under the non-enumerable `attrs` property. Row numbers returned by `detect_stream()` and `match_stream()` are 1-based, as in R.

## Differences with R

* Results match R within 1e-6 (coordinates, distances, speeds) and usually within 1e-12. Timestamps have millisecond resolution (`Date`).
* `resample_dem(method = "bilinear")`: terra computes in 32-bit floats; differences up to 0.1 mm.
* `read_stream_text()` parses GPX directly instead of going through GDAL; it has been checked on the package's example files only.
* `get_dem()` differs by design: R reprojects the Terrain Tiles mosaic to UTM with GDAL before extracting elevation; this port keeps the mosaic in Web Mercator and `ele_correction()` projects the readings to it. Elevations are close to R's but not identical. Pass a DEM made by R (GeoTIFF) for identical results.
* Projections use proj4js, which agrees with PROJ to about 1e-10 m in UTM.
