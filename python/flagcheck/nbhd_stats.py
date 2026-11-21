"""Three-vertex statistics of (a,b)-graphs (no K_a, no independent b-set; edges are red),
used as hypotheses by the neighbourhood-structure constraints ("nstats" in ramsey_run).

For a graph on s vertices let e = #edges, p = #3-sets spanning exactly two edges (induced
paths) and t = #triangles. TRIPLE_POINTS[(a, b)][s] is the set of distinct (e, p, t) over
ALL (a,b)-graphs on s vertices, so any (a,b)-graph on s vertices has its (e, p, t) in that
set. Only orders where a complete list exists are present. Every entry is computed by
`python -m flagcheck.nbhd_stats <dir> <out.json>` from B. McKay's lists "All Ramsey
(a,b)-graphs" (https://users.cecs.anu.edu.au/~bdm/data/ramsey.html, CC BY 4.0) and stored
in data/triple_points.json with the number of graphs read per order.

Soundness: a missing point would make a constraint invalid, so the data file is checked
against the published graph counts (`EXPECTED_COUNTS`) and, when the lists are present,
recomputed in the tests. A point that is present but never attained only weakens.
"""

import gzip
import json
import os
import sys
from math import comb
from typing import Dict, Iterator, Set, Tuple

HERE = os.path.dirname(os.path.abspath(__file__))
DATA_FILE = os.path.join(HERE, "data", "triple_points.json")

# (a, b) -> largest order with a complete list among McKay's files (and (4,5): only 24).
ORDERS = {(3, 4): range(1, 9), (3, 5): range(1, 14), (3, 6): range(1, 18), (4, 4): range(1, 18),
          (4, 5): range(24, 25)}

# Published numbers of (a,b)-graphs on the largest orders (McKay's page; R(4,5): 352366 on
# 24 vertices, Angeltveit-McKay 2016). The data file must agree.
EXPECTED_COUNTS = {(3, 4, 8): 3, (3, 5, 13): 1, (3, 6, 17): 7, (4, 4, 17): 1, (4, 5, 24): 352366}


def adjacency(line: bytes):
    """(s, adjacency bitmasks) of a graph6 line (s <= 62)."""
    line = line.strip()
    s = line[0] - 63
    adj = [0] * s
    bits = 0
    nbits = 0
    data = line[1:]
    pos = 0
    for j in range(1, s):
        for i in range(j):
            if nbits == 0:
                bits = data[pos] - 63
                pos += 1
                nbits = 6
            nbits -= 1
            if (bits >> nbits) & 1:
                adj[i] |= 1 << j
                adj[j] |= 1 << i
    return s, adj


def triple_stats(s: int, adj) -> Tuple[int, int, int]:
    """(e, p, t): edges, induced 2-edge paths, triangles."""
    e2 = 0
    paths = 0
    tri3 = 0
    for c in range(s):
        nb = adj[c]
        d = bin(nb).count("1")
        e2 += d
        inner2 = 0
        x = nb
        while x:
            low = x & -x
            v = low.bit_length() - 1
            inner2 += bin(adj[v] & nb).count("1")
            x ^= low
        inner = inner2 // 2
        tri3 += inner
        paths += comb(d, 2) - inner
    return e2 // 2, paths, tri3 // 3


def all_stats(e: int, p: int, t: int, s: int) -> Tuple[int, int, int, int]:
    """(n0, n1, n2, n3): numbers of 3-sets spanning 0..3 edges."""
    n3, n2 = t, p
    n1 = e * (s - 2) - 2 * n2 - 3 * n3
    n0 = comb(s, 3) - n1 - n2 - n3
    return n0, n1, n2, n3


def _lines(path) -> Iterator[bytes]:
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rb") as fh:
        for line in fh:
            if line.strip():
                yield line


def compute(data_dir: str) -> Dict:
    out = {}
    for (a, b), orders in ORDERS.items():
        for s in orders:
            base = os.path.join(data_dir, f"r{a}{b}_{s}.g6")
            path = base if os.path.exists(base) else base + ".gz"
            pts: Set[Tuple[int, int, int]] = set()
            count = 0
            for line in _lines(path):
                n, adj = adjacency(line)
                if n != s:
                    raise ValueError(f"{path}: graph on {n} vertices")
                pts.add(triple_stats(n, adj))
                count += 1
            out.setdefault(f"{a},{b}", {})[str(s)] = {"count": count, "points": sorted(pts)}
            print(f"({a},{b}) s={s}: {count} graphs, {len(pts)} points", file=sys.stderr, flush=True)
    return out


_TABLE = None


def table() -> Dict[Tuple[int, int], Dict[int, dict]]:
    global _TABLE
    if _TABLE is None:
        with open(DATA_FILE) as fh:
            raw = json.load(fh)
        _TABLE = {tuple(map(int, k.split(","))): {int(s): {"count": v["count"],
                                                             "points": {tuple(p) for p in v["points"]}}
                                                    for s, v in d.items()}
                  for k, d in raw.items()}
    return _TABLE


def triple_points(a: int, b: int, s: int):
    """Set of (e, p, t) over all (a,b)-graphs on s vertices, or None if no complete list.
    a > b uses complements: edges <-> non-edges, so n_j <-> n_{3-j}."""
    if a > b:
        pts = triple_points(b, a, s)
        if pts is None:
            return None
        out = set()
        for e, p, t in pts:
            n0, n1, n2, n3 = all_stats(e, p, t, s)
            out.add((comb(s, 2) - e, n1, n0))      # complement: 2 edges <-> 1 edge, 3 <-> 0
        return out
    entry = table().get((a, b), {}).get(s)
    return None if entry is None else entry["points"]


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    out = compute(argv[0])
    for (a, b, s), n in EXPECTED_COUNTS.items():
        got = out[f"{a},{b}"][str(s)]["count"]
        if got != n:
            raise SystemExit(f"COUNT MISMATCH ({a},{b},{s}): {got} != {n}")
    os.makedirs(os.path.dirname(argv[1]), exist_ok=True)
    with open(argv[1], "w") as fh:
        json.dump(out, fh, indent=0)
    print("WROTE", argv[1])


if __name__ == "__main__":
    main()
