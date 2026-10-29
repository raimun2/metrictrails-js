// Nearest feature as returned by sf::st_nearest_feature (GEOS 3.12 STRtree),
// port of _strtree.py. Exact ties (a reading equidistant to two segments, or
// duplicated points) are broken as GEOS does: sort-tile-recursive packing with
// node capacity 10 (libstdc++ std::sort) and a libstdc++ priority queue.

import { stdSort } from "./libstdcpp.js";

const NODE_CAPACITY = 10;

function parentOf(children) {
  let minx = Infinity, miny = Infinity, maxx = -Infinity, maxy = -Infinity;
  for (const c of children) {
    if (c.minx < minx) minx = c.minx;
    if (c.miny < miny) miny = c.miny;
    if (c.maxx > maxx) maxx = c.maxx;
    if (c.maxy > maxy) maxy = c.maxy;
  }
  return { minx, miny, maxx, maxy, children, item: -1 };
}

const lessX = (a, b) => a.minx + a.maxx < b.minx + b.maxx;
const lessY = (a, b) => a.miny + a.maxy < b.miny + b.maxy;

function build(leaves) {
  let level = leaves.slice();
  while (level.length > 1) {
    const n = level.length;
    const nSlices = Math.ceil(Math.sqrt(Math.ceil(n / NODE_CAPACITY)));
    const perSlice = Math.ceil(n / nSlices);
    stdSort(level, 0, n, lessX);
    const parents = [];
    let start = 0;
    for (let s = 0; s < nSlices; s++) {
      const end = start + Math.min(n - start, perSlice);
      stdSort(level, start, end, lessY);
      for (let c = start; c < end; c += NODE_CAPACITY) {
        parents.push(parentOf(level.slice(c, Math.min(c + NODE_CAPACITY, end))));
      }
      start = end;
    }
    level = parents;
  }
  return level[0];
}

function envDistance(a, px, py) {
  if (a.minx <= px && px <= a.maxx && a.miny <= py && py <= a.maxy) return 0;
  let dx = 0;
  if (a.maxx < px) dx = px - a.maxx;
  else if (a.minx > px) dx = a.minx - px;
  let dy = 0;
  if (a.maxy < py) dy = py - a.maxy;
  else if (a.miny > py) dy = a.miny - py;
  if (dx === 0) return dy;
  if (dy === 0) return dx;
  return Math.sqrt(dx * dx + dy * dy);
}

// GEOS Distance::pointToSegment
export function pointSegmentDistance(px, py, ax, ay, bx, by) {
  if (ax === bx && ay === by) return Math.sqrt((px - ax) * (px - ax) + (py - ay) * (py - ay));
  const len2 = (bx - ax) * (bx - ax) + (by - ay) * (by - ay);
  const r = ((px - ax) * (bx - ax) + (py - ay) * (by - ay)) / len2;
  if (r <= 0) return Math.sqrt((px - ax) * (px - ax) + (py - ay) * (py - ay));
  if (r >= 1) return Math.sqrt((px - bx) * (px - bx) + (py - by) * (py - by));
  const s = ((ay - py) * (bx - ax) - (ax - px) * (by - ay)) / len2;
  return Math.abs(s) * Math.sqrt(len2);
}

// libstdc++ priority_queue with comp(a, b) = a.d > b.d (top = smallest)
function heapPush(h, v) {
  h.push(v);
  let hole = h.length - 1;
  let parent = Math.floor((hole - 1) / 2);
  while (hole > 0 && h[parent].d > v.d) {
    h[hole] = h[parent];
    hole = parent;
    parent = Math.floor((hole - 1) / 2);
  }
  h[hole] = v;
}

function heapPop(h) {
  const last = h.length - 1;
  if (last > 0) {
    const value = h[last];
    h[last] = h[0];
    const length = last;
    let hole = 0;
    let second = 0;
    while (second < Math.floor((length - 1) / 2)) {
      second = 2 * (second + 1);
      if (h[second].d > h[second - 1].d) second--;
      h[hole] = h[second];
      hole = second;
    }
    if ((length & 1) === 0 && second === Math.floor((length - 2) / 2)) {
      second = 2 * (second + 1);
      h[hole] = h[second - 1];
      hole = second - 1;
    }
    let parent = Math.floor((hole - 1) / 2);
    while (hole > 0 && h[parent].d > value.d) {
      h[hole] = h[parent];
      hole = parent;
      parent = Math.floor((hole - 1) / 2);
    }
    h[hole] = value;
  }
  h.pop();
}

function search(root, leafDistance, px, py) {
  const dist = (node) => (node.children ? envDistance(node, px, py) : leafDistance(node.item));
  const heap = [];
  heapPush(heap, { d: dist(root), node: root });
  let lower = Infinity;
  let best = null;
  while (heap.length && lower > 0) {
    const { d, node } = heap[0];
    if (best !== null && d >= lower) break;
    heapPop(heap);
    if (!node.children) {
      lower = d;
      best = node;
    } else {
      for (const child of node.children) {
        const cd = dist(child);
        if (cd < lower) heapPush(heap, { d: cd, node: child });
      }
    }
  }
  return best.item;
}

export class SegmentIndex {
  constructor(rx, ry) {
    this.rx = rx;
    this.ry = ry;
    const leaves = [];
    for (let i = 0; i < rx.length - 1; i++) {
      leaves.push({ minx: Math.min(rx[i], rx[i + 1]), miny: Math.min(ry[i], ry[i + 1]),
        maxx: Math.max(rx[i], rx[i + 1]), maxy: Math.max(ry[i], ry[i + 1]), children: null, item: i });
    }
    this.root = build(leaves);
  }

  // 0-based index of the nearest segment
  nearest(px, py) {
    const { rx, ry } = this;
    return search(this.root, (i) => pointSegmentDistance(px, py, rx[i], ry[i], rx[i + 1], ry[i + 1]), px, py);
  }
}

export class PointIndex {
  constructor(x, y) {
    this.x = x;
    this.y = y;
    this.root = build(x.map((_, i) => ({ minx: x[i], miny: y[i], maxx: x[i], maxy: y[i], children: null, item: i })));
  }

  nearest(px, py) {
    const { x, y } = this;
    return search(this.root, (i) => {
      const dx = px - x[i];
      const dy = py - y[i];
      return Math.sqrt(dx * dx + dy * dy);
    }, px, py);
  }
}
