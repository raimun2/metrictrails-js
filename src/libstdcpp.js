// std::sort as implemented in libstdc++ (GCC), port of _libstdcpp.py.
// GEOS packs its STRtree with std::sort, which is not stable; reproducing the
// exact sequence of swaps makes ties break as in R (sf/GEOS on Linux).

const THRESHOLD = 16;
const lg = (n) => 31 - Math.clz32(n);

function adjustHeap(a, first, hole, length, value, less) {
  const top = hole;
  let second = hole;
  while (second < Math.floor((length - 1) / 2)) {
    second = 2 * (second + 1);
    if (less(a[first + second], a[first + second - 1])) second--;
    a[first + hole] = a[first + second];
    hole = second;
  }
  if ((length & 1) === 0 && second === Math.floor((length - 2) / 2)) {
    second = 2 * (second + 1);
    a[first + hole] = a[first + second - 1];
    hole = second - 1;
  }
  let parent = Math.floor((hole - 1) / 2);
  while (hole > top && less(a[first + parent], value)) {
    a[first + hole] = a[first + parent];
    hole = parent;
    parent = Math.floor((hole - 1) / 2);
  }
  a[first + hole] = value;
}

function makeHeap(a, first, last, less) {
  const length = last - first;
  if (length < 2) return;
  let parent = Math.floor((length - 2) / 2);
  for (;;) {
    adjustHeap(a, first, parent, length, a[first + parent], less);
    if (parent === 0) return;
    parent--;
  }
}

function popHeap(a, first, last, result, less) {
  const value = a[result];
  a[result] = a[first];
  adjustHeap(a, first, 0, last - first, value, less);
}

function heapSelect(a, first, middle, last, less) {
  makeHeap(a, first, middle, less);
  for (let i = middle; i < last; i++) if (less(a[i], a[first])) popHeap(a, first, middle, i, less);
}

function sortHeap(a, first, last, less) {
  while (last - first > 1) {
    last--;
    popHeap(a, first, last, last, less);
  }
}

function swap(a, i, j) {
  const t = a[i];
  a[i] = a[j];
  a[j] = t;
}

function moveMedianToFirst(a, result, x, y, z, less) {
  if (less(a[x], a[y])) {
    if (less(a[y], a[z])) swap(a, result, y);
    else if (less(a[x], a[z])) swap(a, result, z);
    else swap(a, result, x);
  } else if (less(a[x], a[z])) swap(a, result, x);
  else if (less(a[y], a[z])) swap(a, result, z);
  else swap(a, result, y);
}

function unguardedPartition(a, first, last, pivot, less) {
  for (;;) {
    while (less(a[first], a[pivot])) first++;
    last--;
    while (less(a[pivot], a[last])) last--;
    if (!(first < last)) return first;
    swap(a, first, last);
    first++;
  }
}

function unguardedPartitionPivot(a, first, last, less) {
  const mid = first + Math.floor((last - first) / 2);
  moveMedianToFirst(a, first, first + 1, mid, last - 1, less);
  return unguardedPartition(a, first + 1, last, first, less);
}

function unguardedLinearInsert(a, last, less) {
  const val = a[last];
  let next = last - 1;
  while (less(val, a[next])) {
    a[last] = a[next];
    last = next;
    next--;
  }
  a[last] = val;
}

function insertionSort(a, first, last, less) {
  if (first === last) return;
  for (let i = first + 1; i < last; i++) {
    if (less(a[i], a[first])) {
      const val = a[i];
      for (let k = i; k > first; k--) a[k] = a[k - 1];
      a[first] = val;
    } else unguardedLinearInsert(a, i, less);
  }
}

function introsortLoop(a, first, last, depth, less) {
  while (last - first > THRESHOLD) {
    if (depth === 0) {
      heapSelect(a, first, last, last, less);
      sortHeap(a, first, last, less);
      return;
    }
    depth--;
    const cut = unguardedPartitionPivot(a, first, last, less);
    introsortLoop(a, cut, last, depth, less);
    last = cut;
  }
}

// std::sort(a + first, a + last, less), in place
export function stdSort(a, first, last, less) {
  if (first === last) return;
  introsortLoop(a, first, last, lg(last - first) * 2, less);
  if (last - first > THRESHOLD) {
    insertionSort(a, first, first + THRESHOLD, less);
    for (let i = first + THRESHOLD; i < last; i++) unguardedLinearInsert(a, i, less);
  } else insertionSort(a, first, last, less);
}
