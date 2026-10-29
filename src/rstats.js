// Replicas of the R routines used by metrictrails (port of _rstats.py).
// The R algorithms are reimplemented from the R sources so that results
// match R: quantile type 7, approx, findInterval, stats::filter, loess with
// surface = "direct", optimize (Brent) and optim Nelder-Mead.

export const DBL_EPSILON = 2.220446049250313e-16;

const isNA = (v) => v === null || v === undefined || Number.isNaN(v);

// ---- summaries -----------------------------------------------------------------

// Neumaier compensated sum (close to R's long double accumulation)
export function rSum(x) {
  let s = 0;
  let c = 0;
  for (const v of x) {
    const t = s + v;
    if (Math.abs(s) >= Math.abs(v)) c += (s - t) + v;
    else c += (v - t) + s;
    s = t;
  }
  return s + c;
}

export function rMean(x) {
  const n = x.length;
  if (n === 0) return NaN;
  for (const v of x) if (Number.isNaN(v)) return NaN;
  const s = rSum(x) / n;
  if (!Number.isFinite(s)) return s;
  return s + rSum(x.map((v) => v - s)) / n;
}

export const naOmit = (x) => Array.from(x).filter((v) => !isNA(v));

export function quantile7(x, probs) {
  const s = naOmit(x).sort((a, b) => a - b);
  const n = s.length;
  const one = !Array.isArray(probs);
  const ps = one ? [probs] : probs;
  const out = ps.map((p) => {
    if (n === 0) return NaN;
    const index = 1 + Math.max(n - 1, 0) * p;
    const lo = Math.floor(index);
    const hi = Math.ceil(index);
    let qs = s[lo - 1];
    const h = index - lo;
    if (index > lo && s[hi - 1] !== qs) qs = (1 - h) * qs + h * s[hi - 1];
    return qs;
  });
  return one ? out[0] : out;
}

export function rMedian(x) {
  const s = naOmit(x).sort((a, b) => a - b);
  const n = s.length;
  if (n === 0) return NaN;
  const half = Math.floor((n + 1) / 2);
  if (n % 2 === 1) return s[half - 1];
  return rMean([s[half - 1], s[half]]);
}

export function rSd(x) {
  const v = naOmit(x);
  if (v.length < 2) return NaN;
  const m = rMean(v);
  return Math.sqrt(rSum(v.map((a) => (a - m) ** 2)) / (v.length - 1));
}

export function rIqr(x) {
  const q = quantile7(x, [0.25, 0.75]);
  return q[1] - q[0];
}

export function rankAverage(x) {
  const order = x.map((_, i) => i).sort((a, b) => x[a] - x[b] || a - b);
  const ranks = new Array(x.length);
  let i = 0;
  while (i < x.length) {
    let j = i;
    while (j + 1 < x.length && x[order[j + 1]] === x[order[i]]) j++;
    for (let k = i; k <= j; k++) ranks[order[k]] = (i + j) / 2 + 1;
    i = j + 1;
  }
  return ranks;
}

export function rCor(x, y) {
  if (x.length < 2 || x.some(isNA) || y.some(isNA)) return NaN;
  const mx = rMean(x);
  const my = rMean(y);
  const sxx = rSum(x.map((v) => (v - mx) ** 2));
  const syy = rSum(y.map((v) => (v - my) ** 2));
  if (sxx === 0 || syy === 0) return NaN;
  return rSum(x.map((v, i) => (v - mx) * (y[i] - my))) / Math.sqrt(sxx * syy);
}

export function corSpearman(x, y) {
  if (x.some(isNA) || y.some(isNA)) return NaN;
  return rCor(rankAverage(x), rankAverage(y));
}

// ---- sequences and intervals --------------------------------------------------

export function seqBy(from, to, by) {
  if (to < from) throw new Error("wrong sign in 'by' argument");
  const del = to - from;
  if (del === 0 && to === 0) return [to];
  const dd = Math.abs(del) / Math.max(Math.abs(to), Math.abs(from));
  if (dd < 100 * DBL_EPSILON) return [from];
  const n = Math.trunc(del / by + 1e-10);
  const out = new Array(n + 1);
  for (let i = 0; i <= n; i++) out[i] = Math.min(from + i * by, to);
  return out;
}

export function seqLength(from, to, n) {
  if (n === 1) return [from];
  const by = (to - from) / (n - 1);
  const out = new Array(n);
  for (let i = 0; i < n; i++) out[i] = from + i * by;
  out[n - 1] = to;
  return out;
}

// number of values of sorted vec that are <= v
function bisectRight(vec, v) {
  let lo = 0;
  let hi = vec.length;
  while (lo < hi) {
    const mid = (lo + hi) >> 1;
    if (vec[mid] <= v) lo = mid + 1;
    else hi = mid;
  }
  return lo;
}

// findInterval(x, vec, rightmost.closed): 1-based, 0 below the first value
export function findInterval(x, vec, rightmostClosed = false) {
  const n = vec.length;
  return x.map((v) => {
    let i = bisectRight(vec, v);
    if (rightmostClosed && n && i === n && v === vec[n - 1]) i = n - 1;
    return i;
  });
}

// .bincode(x, breaks, right = TRUE, include.lowest = TRUE): 1-based, 0 = NA
export function bincode(x, breaks) {
  const nb = breaks.length;
  return x.map((v) => {
    if (isNA(v)) return 0;
    let lo = 0;
    let hi = nb - 1;
    if (v < breaks[lo] || breaks[hi] < v) return 0;
    while (hi - lo >= 2) {
      const mid = (hi + lo) >> 1;
      if (v > breaks[mid]) lo = mid;
      else hi = mid;
    }
    return lo + 1;
  });
}

// ---- interpolation ----------------------------------------------------------------

function regularize(x, y, ties) {
  let xs = [];
  let ys = [];
  for (let i = 0; i < x.length; i++) {
    if (!isNA(x[i]) && !isNA(y[i])) {
      xs.push(x[i]);
      ys.push(y[i]);
    }
  }
  if (ties === "ordered") return [xs, ys];
  let sorted = true;
  for (let i = 1; i < xs.length; i++) if (xs[i] < xs[i - 1]) { sorted = false; break; }
  if (!sorted) {
    const o = xs.map((_, i) => i).sort((a, b) => xs[a] - xs[b] || a - b);
    xs = o.map((i) => xs[i]);
    ys = o.map((i) => ys[i]);
  }
  let dup = false;
  for (let i = 1; i < xs.length; i++) if (xs[i] === xs[i - 1]) { dup = true; break; }
  if (!dup) return [xs, ys];
  const ux = [];
  const uy = [];
  let i = 0;
  while (i < xs.length) {
    let j = i;
    while (j + 1 < xs.length && xs[j + 1] === xs[i]) j++;
    ux.push(xs[i]);
    uy.push(rMean(ys.slice(i, j + 1)));
    i = j + 1;
  }
  return [ux, uy];
}

// stats::approx(x, y, xout, method = "linear", ties, rule)$y
export function approx(x, y, xout, { ties = "mean", rule = 1 } = {}) {
  const [xs, ys] = regularize(x, y, ties);
  const n = xs.length;
  if (n < 2) throw new Error("need at least two non-NA values to interpolate");
  const ylow = rule === 2 ? ys[0] : NaN;
  const yhigh = rule === 2 ? ys[n - 1] : NaN;
  return xout.map((v) => {
    if (isNA(v)) return NaN;
    let i = 0;
    let j = n - 1;
    if (v < xs[i]) return ylow;
    if (v > xs[j]) return yhigh;
    while (i < j - 1) {
      const ij = (i + j) >> 1;
      if (v < xs[ij]) j = ij;
      else i = ij;
    }
    if (v === xs[j]) return ys[j];
    if (v === xs[i]) return ys[i];
    return ys[i] + (ys[j] - ys[i]) * ((v - xs[i]) / (xs[j] - xs[i]));
  });
}

// ---- filters and smoothers --------------------------------------------------------

// stats::filter(x, rep(1 / window, window), sides = 2)
export function convFilterCentered(x, window) {
  const n = x.length;
  const f = 1 / window;
  const o = Math.floor(window / 2);
  const out = new Array(n).fill(NaN);
  for (let i = 0; i < n; i++) {
    if (i + o - (window - 1) < 0 || i + o >= n) continue;
    let z = 0;
    let bad = false;
    for (let j = 0; j < window; j++) {
      const tmp = x[i + o - j];
      if (isNA(tmp)) { bad = true; break; }
      z += f * tmp;
    }
    if (!bad) out[i] = z;
  }
  return out;
}

// Least squares by Householder QR: coefficients of min ||X b - y||
export function lstsq(X, y) {
  const n = X.length;
  const p = X[0].length;
  const A = X.map((r) => r.slice());
  const b = y.slice();
  const R = [];
  for (let k = 0; k < p; k++) {
    let norm = 0;
    for (let i = k; i < n; i++) norm += A[i][k] * A[i][k];
    norm = Math.sqrt(norm);
    if (norm === 0) continue;
    const alpha = A[k][k] > 0 ? -norm : norm;
    const v = new Array(n).fill(0);
    for (let i = k; i < n; i++) v[i] = A[i][k];
    v[k] -= alpha;
    let vv = 0;
    for (let i = k; i < n; i++) vv += v[i] * v[i];
    if (vv === 0) continue;
    for (let j = k; j < p; j++) {
      let s = 0;
      for (let i = k; i < n; i++) s += v[i] * A[i][j];
      s = (2 * s) / vv;
      for (let i = k; i < n; i++) A[i][j] -= s * v[i];
    }
    let s = 0;
    for (let i = k; i < n; i++) s += v[i] * b[i];
    s = (2 * s) / vv;
    for (let i = k; i < n; i++) b[i] -= s * v[i];
  }
  for (let k = 0; k < p; k++) R.push(A[k].slice());
  const beta = new Array(p).fill(0);
  const tol = 1e-12 * Math.max(...R.map((r, i) => Math.abs(r[i])));
  for (let k = p - 1; k >= 0; k--) {
    let s = b[k];
    for (let j = k + 1; j < p; j++) s -= R[k][j] * beta[j];
    beta[k] = Math.abs(R[k][k]) > tol ? s / R[k][k] : 0;
  }
  return beta;
}

// Weighted least squares; returns [coefficients, weighted residuals]
export function wls(X, y, w = null) {
  const sw = y.map((_, i) => (w ? Math.sqrt(w[i]) : 1));
  const Xw = X.map((r, i) => r.map((v) => v * sw[i]));
  const yw = y.map((v, i) => v * sw[i]);
  const beta = lstsq(Xw, yw);
  const resid = X.map((r, i) => (y[i] - r.reduce((s, v, j) => s + v * beta[j], 0)) * sw[i]);
  return [beta, resid];
}

// k-th smallest value (0-based) of an array, not modifying it
function kthSmallest(a, k) {
  const s = Float64Array.from(a).sort();
  return s[k];
}

// loess(v ~ t, span, degree = 2, family = "gaussian", surface = "direct")
export function loessDirect(t, v, span, newdata) {
  const n = t.length;
  const q = Math.floor(n * span + 1e-5);
  if (q <= 0) throw new Error("span is too small");
  const dist = new Float64Array(n);
  return newdata.map((x0) => {
    if (isNA(x0)) return NaN;
    for (let i = 0; i < n; i++) dist[i] = Math.abs(t[i] - x0);
    const h = kthSmallest(dist, q - 1) * Math.max(1, span);
    if (h <= 0) return NaN;
    const X = [];
    const y = [];
    for (let i = 0; i < n; i++) {
      if (dist[i] < h) {
        const w = (1 - (dist[i] / h) ** 3) ** 3;
        if (w > 0) {
          const sw = Math.sqrt(w);
          const d = t[i] - x0;
          X.push([sw, sw * d, sw * d * d]);
          y.push(sw * v[i]);
        }
      }
    }
    return lstsq(X, y)[0];
  });
}

export const dnorm = (x, sd) => Math.exp(-0.5 * (x / sd) ** 2) / (Math.sqrt(2 * Math.PI) * sd);

// ---- optimisers -----------------------------------------------------------------------

// stats::optimize(f, c(lower, upper))$minimum (Brent's fmin, as in R)
export function optimizeBrent(f, lower, upper, tol = DBL_EPSILON ** 0.25) {
  const c = (3 - Math.sqrt(5)) * 0.5;
  const eps = Math.sqrt(DBL_EPSILON);
  let a = lower;
  let b = upper;
  let v = a + c * (b - a);
  let w = v;
  let x = v;
  let d = 0;
  let e = 0;
  let fx = f(x);
  let fv = fx;
  let fw = fx;
  const tol3 = tol / 3;
  for (;;) {
    const xm = (a + b) * 0.5;
    const tol1 = eps * Math.abs(x) + tol3;
    const t2 = tol1 * 2;
    if (Math.abs(x - xm) <= t2 - (b - a) * 0.5) break;
    let p = 0;
    let q = 0;
    let r = 0;
    if (Math.abs(e) > tol1) {
      r = (x - w) * (fx - fv);
      q = (x - v) * (fx - fw);
      p = (x - v) * q - (x - w) * r;
      q = (q - r) * 2;
      if (q > 0) p = -p;
      else q = -q;
      r = e;
      e = d;
    }
    let u;
    if (Math.abs(p) >= Math.abs(q * 0.5 * r) || p <= q * (a - x) || p >= q * (b - x)) {
      e = x < xm ? b - x : a - x;
      d = c * e;
    } else {
      d = p / q;
      u = x + d;
      if (u - a < t2 || b - u < t2) {
        d = tol1;
        if (x >= xm) d = -d;
      }
    }
    if (Math.abs(d) >= tol1) u = x + d;
    else if (d > 0) u = x + tol1;
    else u = x - tol1;
    const fu = f(u);
    if (fu <= fx) {
      if (u < x) b = x;
      else a = x;
      v = w; w = x; x = u;
      fv = fw; fw = fx; fx = fu;
    } else {
      if (u < x) a = u;
      else b = u;
      if (fu <= fw || w === x) {
        v = w; fv = fw;
        w = u; fw = fu;
      } else if (fu <= fv || v === x || v === w) {
        v = u; fv = fu;
      }
    }
  }
  return x;
}

// stats::optim(par, fn, method = "Nelder-Mead") (R's nmmin). Returns {par, value, fail}
export function nelderMead(fn, par, { maxit = 500, reltol = Math.sqrt(DBL_EPSILON), abstol = -Infinity } = {}) {
  const alpha = 1;
  const beta = 0.5;
  const gamma = 2;
  const big = 1.0e35;
  const n = par.length;
  let bvec = par.slice();
  const n1 = n + 1;
  const C = n + 2;
  const P = Array.from({ length: n1 }, () => new Array(C).fill(0));
  let f = fn(bvec.slice());
  if (!Number.isFinite(f)) throw new Error("function cannot be evaluated at initial parameters");
  let funcount = 1;
  const convtol = reltol * (Math.abs(f) + reltol);
  P[n1 - 1][0] = f;
  for (let i = 0; i < n; i++) P[i][0] = bvec[i];
  let L = 1;
  let size = 0;
  let step = 0;
  for (let i = 0; i < n; i++) if (0.1 * Math.abs(bvec[i]) > step) step = 0.1 * Math.abs(bvec[i]);
  if (step === 0) step = 0.1;
  for (let j = 2; j <= n1; j++) {
    for (let i = 0; i < n; i++) P[i][j - 1] = bvec[i];
    let trystep = step;
    while (P[j - 2][j - 1] === bvec[j - 2]) {
      P[j - 2][j - 1] = bvec[j - 2] + trystep;
      trystep *= 10;
    }
    size += trystep;
  }
  let oldsize = size;
  let calcvert = true;
  let fail = 0;
  const evalf = (b) => {
    let v = fn(b.slice());
    if (!Number.isFinite(v)) v = big;
    funcount++;
    return v;
  };
  for (;;) {
    if (calcvert) {
      for (let j = 0; j < n1; j++) {
        if (j + 1 !== L) {
          for (let i = 0; i < n; i++) bvec[i] = P[i][j];
          P[n1 - 1][j] = evalf(bvec);
        }
      }
      calcvert = false;
    }
    let VL = P[n1 - 1][L - 1];
    let VH = VL;
    let H = L;
    for (let j = 1; j <= n1; j++) {
      if (j !== L) {
        const fj = P[n1 - 1][j - 1];
        if (fj < VL) { L = j; VL = fj; }
        if (fj > VH) { H = j; VH = fj; }
      }
    }
    if (VH <= VL + convtol || VL <= abstol) break;
    for (let i = 0; i < n; i++) {
      let temp = -P[i][H - 1];
      for (let j = 0; j < n1; j++) temp += P[i][j];
      P[i][C - 1] = temp / n;
    }
    for (let i = 0; i < n; i++) bvec[i] = (1 + alpha) * P[i][C - 1] - alpha * P[i][H - 1];
    f = evalf(bvec);
    const VR = f;
    if (VR < VL) {
      P[n1 - 1][C - 1] = f;
      for (let i = 0; i < n; i++) {
        const g = gamma * bvec[i] + (1 - gamma) * P[i][C - 1];
        P[i][C - 1] = bvec[i];
        bvec[i] = g;
      }
      f = evalf(bvec);
      if (f < VR) {
        for (let i = 0; i < n; i++) P[i][H - 1] = bvec[i];
        P[n1 - 1][H - 1] = f;
      } else {
        for (let i = 0; i < n; i++) P[i][H - 1] = P[i][C - 1];
        P[n1 - 1][H - 1] = VR;
      }
    } else {
      if (VR < VH) {
        for (let i = 0; i < n; i++) P[i][H - 1] = bvec[i];
        P[n1 - 1][H - 1] = VR;
      }
      for (let i = 0; i < n; i++) bvec[i] = (1 - beta) * P[i][H - 1] + beta * P[i][C - 1];
      f = evalf(bvec);
      if (f < P[n1 - 1][H - 1]) {
        for (let i = 0; i < n; i++) P[i][H - 1] = bvec[i];
        P[n1 - 1][H - 1] = f;
      } else if (VR >= VH) {
        calcvert = true;
        size = 0;
        for (let j = 0; j < n1; j++) {
          if (j + 1 !== L) {
            for (let i = 0; i < n; i++) {
              P[i][j] = beta * (P[i][j] - P[i][L - 1]) + P[i][L - 1];
              size += Math.abs(P[i][j] - P[i][L - 1]);
            }
          }
        }
        if (size < oldsize) oldsize = size;
        else { fail = 10; break; }
      }
    }
    if (funcount > maxit) break;
  }
  if (funcount > maxit) fail = 1;
  const out = [];
  for (let i = 0; i < n; i++) out.push(P[i][L - 1]);
  return { par: out, value: P[n1 - 1][L - 1], fail };
}
