"""Minimum and maximum edge counts of (a,b)-graphs (no K_a, no independent b-set) on n
vertices, used as hypotheses by the neighbourhood edge constraints in ramsey.py.

EDGE_BOUNDS[(a, b)][n] = (min, max) number of edges over all (a,b)-graphs on n vertices,
for a <= b. Orders with no (a,b)-graph are absent. Every entry was computed by
`python -m flagcheck.edge_tables <dir>` from B. McKay's Ramsey graph lists
(https://users.cecs.anu.edu.au/~bdm/data/ramsey.html, CC BY 4.0):

  * (3,4), (3,5), (3,6), (4,4): the complete lists "All Ramsey(a,b)-graphs";
  * (4,5), n = 4..23: r45extreme.tar.gz, "the complete sets of Ramsey(4,5)-graphs with
    the smallest few edge counts and the largest few edge counts" (so the minimum and the
    maximum are exact even though the middle of the range is not listed);
  * (4,5), n = 24: r45_24.g6, all 352366 graphs (McKay-Radziszowski 1995, completed by
    Angeltveit-McKay 2016).

The minima for (3,k) are cross-checked in the tests against the independently published
Table 3 of J. Goedgebeur, S. Radziszowski, "New computational upper bounds for Ramsey
numbers R(3,k)", Electron. J. Combin. 20(1) (2013) #P30.

Soundness: a constraint "the red neighbourhood has >= min edges" is only valid if min is
a true lower bound, so an entry that is too large would be unsound; an entry that is too
small only weakens the constraint (likewise for max).
"""

import gzip
import os
import sys
import tarfile
from math import comb
from typing import Dict, Iterator, Tuple

# fmt: off
EDGE_BOUNDS: Dict[Tuple[int, int], Dict[int, Tuple[int, int]]] = {
    (3, 3): {1: (0, 0), 2: (0, 1), 3: (1, 2), 4: (2, 4), 5: (5, 5)},      # brute force (`brute_bounds`)
    (3, 4): {1: (0, 0), 2: (0, 1), 3: (0, 2), 4: (1, 4), 5: (2, 6), 6: (3, 9), 7: (6, 10), 8: (10, 12)},
    (3, 5): {1: (0, 0), 2: (0, 1), 3: (0, 2), 4: (0, 4), 5: (1, 6), 6: (2, 9), 7: (3, 12), 8: (4, 16),
             9: (7, 17), 10: (10, 20), 11: (15, 22), 12: (20, 24), 13: (26, 26)},
    (3, 6): {1: (0, 0), 2: (0, 1), 3: (0, 2), 4: (0, 4), 5: (0, 6), 6: (1, 9), 7: (2, 12), 8: (3, 16),
             9: (4, 20), 10: (5, 25), 11: (8, 26), 12: (11, 29), 13: (15, 32), 14: (20, 35),
             15: (25, 37), 16: (32, 40), 17: (40, 42)},
    (4, 4): {1: (0, 0), 2: (0, 1), 3: (0, 3), 4: (1, 5), 5: (2, 8), 6: (3, 12), 7: (5, 16), 8: (7, 21),
             9: (9, 27), 10: (14, 31), 11: (19, 36), 12: (26, 40), 13: (33, 45), 14: (41, 50),
             15: (50, 55), 16: (60, 60), 17: (68, 68)},
    (4, 5): {4: (0, 5), 5: (1, 8), 6: (2, 12), 7: (3, 16), 8: (4, 21), 9: (6, 27), 10: (8, 33),
             11: (10, 40), 12: (12, 48), 13: (17, 53), 14: (22, 60), 15: (27, 66), 16: (32, 72),
             17: (41, 79), 18: (50, 85), 19: (57, 92), 20: (68, 100), 21: (77, 107), 22: (88, 114),
             23: (101, 122), 24: (116, 132)},
}
# (4,5) orders 1..3: every graph qualifies (no K4 possible, independence < 5).
EDGE_BOUNDS[(4, 5)].update({1: (0, 0), 2: (0, 1), 3: (0, 3)})
# fmt: on


def edge_bounds(a: int, b: int, n: int) -> Tuple[int, int]:
    """(min, max) number of edges of an (a,b)-graph on n vertices. Raises KeyError when no
    such graph exists (or the order is not covered by the data)."""
    if n < 0:
        raise KeyError(n)
    if a == 1 or b == 1:
        if n == 0:
            return (0, 0)
        raise KeyError(n)
    if a == 2:              # no edge at all, n < b
        if n < b:
            return (0, 0)
        raise KeyError(n)
    if b == 2:              # complete graph, n < a
        if n < a:
            return (comb(n, 2), comb(n, 2))
        raise KeyError(n)
    if n == 0:
        return (0, 0)
    if a > b:               # complements: (a,b) on n vertices <-> (b,a)
        lo, hi = edge_bounds(b, a, n)
        return (comb(n, 2) - hi, comb(n, 2) - lo)
    return EDGE_BOUNDS[(a, b)][n]


# ---------------------------------------------------------------------------
# Recomputation from McKay's files
# ---------------------------------------------------------------------------

def parse_g6(line: bytes):
    """(n, edges as a list of pairs) for a graph6 line (n <= 62)."""
    line = line.strip()
    if line.startswith(b">>graph6<<"):
        line = line[10:]
    n = line[0] - 63
    if not 0 <= n <= 62:
        raise ValueError("only n <= 62 supported")
    bits = []
    for ch in line[1:]:
        v = ch - 63
        bits.extend((v >> (5 - i)) & 1 for i in range(6))
    edges = []
    t = 0
    for j in range(1, n):
        for i in range(j):
            if bits[t]:
                edges.append((i, j))
            t += 1
    if any(bits[t:]):
        raise ValueError("nonzero padding")
    return n, edges


def g6_edge_count(line: bytes) -> Tuple[int, int]:
    """(n, number of edges), quickly: popcount of the data bytes (padding bits are 0)."""
    line = line.strip()
    n = line[0] - 63
    return n, sum(bin(ch - 63).count("1") for ch in line[1:])


def is_ab_graph(n, edges, a, b) -> bool:
    """No K_a and no independent b-set (brute force with pruning)."""
    adj = [set() for _ in range(n)]
    for i, j in edges:
        adj[i].add(j)
        adj[j].add(i)
    non = [set(range(n)) - adj[v] - {v} for v in range(n)]

    def has_clique(nb, k, size, cand):
        if size == k:
            return True
        if size + len(cand) < k:
            return False
        for v in sorted(cand):
            if has_clique(nb, k, size + 1, {u for u in cand if u > v} & nb[v]):
                return True
        return False

    return not has_clique(adj, a, 0, set(range(n))) and not has_clique(non, b, 0, set(range(n)))


def brute_bounds(a: int, b: int, n: int):
    """(min, max) edges over ALL labelled graphs on n vertices that are (a,b)-graphs, or
    None if there is none. Exhaustive; only for n <= 6."""
    from itertools import combinations as _comb

    if n > 6:
        raise ValueError("too many graphs")
    pairs = list(_comb(range(n), 2))
    found = []
    for mask in range(1 << len(pairs)):
        edges = [p for t, p in enumerate(pairs) if mask >> t & 1]
        if is_ab_graph(n, edges, a, b):
            found.append(len(edges))
    return (min(found), max(found)) if found else None


def _lines(path) -> Iterator[bytes]:
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rb") as fh:
        for line in fh:
            if line.strip():
                yield line


def compute(data_dir: str, sample_every: int = 997) -> Dict:
    """Recompute EDGE_BOUNDS from the files in data_dir. Also checks that the vertex count
    of every graph matches its file and that a sample of graphs really are (a,b)-graphs.
    Returns {(a,b): {n: (min, max, count)}}."""
    out: Dict = {(3, 3): {n: brute_bounds(3, 3, n) + (None,) for n in range(1, 6)}}
    checked = 0

    def account(key, n_expected, line, idx, a, b):
        nonlocal checked
        n, e = g6_edge_count(line)
        if n != n_expected:
            raise ValueError(f"{key}: graph with {n} vertices in the n={n_expected} file")
        if idx % sample_every == 0:
            nn, edges = parse_g6(line)
            if len(edges) != e or not is_ab_graph(nn, edges, a, b):
                raise ValueError(f"{key} n={n}: graph {idx} is not an ({a},{b})-graph")
            checked += 1
        return e

    for (a, b), top in (((3, 4), 8), ((3, 5), 13), ((3, 6), 17), ((4, 4), 17)):
        table = {}
        for n in range(1, top + 1):
            base = os.path.join(data_dir, f"r{a}{b}_{n}.g6")
            path = base if os.path.exists(base) else base + ".gz"
            es = [account((a, b), n, line, i, a, b) for i, line in enumerate(_lines(path))]
            table[n] = (min(es), max(es), len(es))
        out[(a, b)] = table

    # (4,5): extreme sets for n = 4..23 (file names r45<n>.<edges>.g6), all graphs for 24.
    table = {}
    with tarfile.open(os.path.join(data_dir, "r45extreme.tar.gz")) as tar:
        per_n: Dict[int, list] = {}
        for member in tar.getmembers():
            if not member.isfile():
                continue
            # r45<n>.<e>.g6 (all graphs with n vertices and e edges), or r45<n>.g6 for
            # n <= 9 (all graphs with n vertices)
            name = os.path.basename(member.name)
            stem = name[len("r45"):-len(".g6")]
            n_str, _, e_str = stem.partition(".")
            n_file = int(n_str)
            lines = [ln for ln in tar.extractfile(member).read().splitlines() if ln.strip()]
            for i, line in enumerate(lines):
                e = account((4, 5), n_file, line, i, 4, 5)
                if e_str and e != int(e_str):
                    raise ValueError(f"{name}: graph with {e} edges")
                per_n.setdefault(n_file, []).append(e)
        for n, es in per_n.items():
            table[n] = (min(es), max(es), None)
    es = [account((4, 5), 24, line, i, 4, 5) for i, line in enumerate(_lines(os.path.join(data_dir, "r45_24.g6")))]
    table[24] = (min(es), max(es), len(es))
    for n in (1, 2, 3):             # not in the archive: every graph on <= 3 vertices qualifies
        table[n] = (0, comb(n, 2), None)
    out[(4, 5)] = table
    print(f"sampled {checked} graphs: all are genuine (a,b)-graphs", file=sys.stderr)
    return out


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    computed = compute(argv[0])
    ok = True
    for key in sorted(computed):
        print(f"    {key}: {{")
        for n in sorted(computed[key]):
            lo, hi, cnt = computed[key][n]
            print(f"        {n}: ({lo}, {hi}),  # {cnt if cnt is not None else 'extreme sets'} graphs")
            if EDGE_BOUNDS and EDGE_BOUNDS.get(key, {}).get(n) != (lo, hi):
                ok = False
                print(f"MISMATCH {key} n={n}: table {EDGE_BOUNDS.get(key, {}).get(n)}", file=sys.stderr)
        print("    },")
    if EDGE_BOUNDS:
        extra = {(k, n) for k in EDGE_BOUNDS for n in EDGE_BOUNDS[k]} - {(k, n) for k in computed for n in computed[k]}
        if extra:
            ok = False
            print(f"MISMATCH: table entries not backed by data: {sorted(extra)}", file=sys.stderr)
    print("TABLE MATCHES DATA" if ok and EDGE_BOUNDS else "TABLE NOT CHECKED" if not EDGE_BOUNDS else "TABLE MISMATCH")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
