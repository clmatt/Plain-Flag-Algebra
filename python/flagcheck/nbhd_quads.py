"""Four-vertex statistics of neighbourhoods, used by the "quad" constraints of ramsey_run.

For a graph H on s vertices (edges = red) let c(H) in Z^11 count its 4-subsets by the
induced graph (CLASSES below: the 11 graphs on 4 vertices, named by their sorted degree
sequence, which determines a 4-vertex graph). A *quad cut* is a direction a in Q^11 and,
for each red degree d of the window, an offset b_d with

    a . c(H) >= b_d   for every (a',b')-graph H on s(d) vertices,

where H is the red neighbourhood ((k-1,l)-graph on s = d vertices) or the blue one
((k,l-1)-graph on s = N-1-d) of a root of red degree d. b_d is

  * min over McKay's complete list of all (a',b')-graphs on s vertices ("list"), or
  * C(s,4) min_{H admissible} a_H ("trivial", always valid: c >= 0, sum c = C(s,4), and
    c_H = 0 for a class H containing a red K_a' or an independent b'-set).

A cut may instead be *localised* at one degree d ("degree": d, one offset): then the
constraint is L_d(R) (a . Q - 24 b_d / N^4) >= 0 with L_d the window's Lagrange basis
polynomial, which is 1 at roots of degree d and 0 at the others (neighbourhoods.quad_cut_expr).

The directions were found by separating the n = 7 SDP solution from the convex hull of the
listed vectors (scratch tooling, recorded in docs/LOG.md); any direction gives a valid cut, so
only the offsets need checking. They are stored in data/quad_cuts.json and recomputed from
the lists by `python -m flagcheck.nbhd_quads check <dir>` and the tests (when the lists are
present). Lists: B. McKay, "All Ramsey (a,b)-graphs",
https://users.cecs.anu.edu.au/~bdm/data/ramsey.html (CC BY 4.0); the graph counts read
must equal nbhd_stats.EXPECTED_COUNTS / the published ones.
"""

import gzip
import itertools
import json
import os
import sys
from fractions import Fraction
from math import comb, gcd
from typing import Dict, List, Set, Tuple

HERE = os.path.dirname(os.path.abspath(__file__))
CUTS_FILE = os.path.join(HERE, "data", "quad_cuts.json")

CLASSES = [(0, 0, 0, 0), (1, 1, 0, 0), (1, 1, 1, 1), (2, 1, 1, 0), (2, 2, 2, 0), (2, 2, 1, 1),
           (3, 1, 1, 1), (2, 2, 2, 2), (3, 2, 2, 1), (3, 3, 2, 2), (3, 3, 3, 3)]
# A representative edge set of each class on vertices 0..3.
REPS = [[], [(0, 1)], [(0, 1), (2, 3)], [(0, 1), (1, 2)], [(0, 1), (1, 2), (0, 2)], [(0, 1), (1, 2), (2, 3)],
        [(0, 1), (0, 2), (0, 3)], [(0, 1), (1, 2), (2, 3), (0, 3)], [(0, 1), (0, 2), (1, 2), (0, 3)],
        [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3)], [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]]
PAIRS = list(itertools.combinations(range(4), 2))
# Published numbers of all (a,b)-graphs on s vertices (McKay's page).
LIST_COUNTS = {(3, 6, 15): 64732, (3, 6, 16): 2576, (3, 6, 17): 7, (4, 5, 24): 352366}


def class_of(edges) -> int:
    deg = [0] * 4
    for i, j in edges:
        deg[i] += 1
        deg[j] += 1
    return CLASSES.index(tuple(sorted(deg, reverse=True)))


def admissible_classes(a: int, b: int) -> List[int]:
    """Classes that can occur in an (a,b)-graph: no K_a, no independent b-set (4 vertices)."""
    out = []
    for h, E in enumerate(REPS):
        Es = set(E)
        if a <= 4 and any(all((x, y) in Es for x, y in itertools.combinations(S, 2))
                          for S in itertools.combinations(range(4), a)):
            continue
        if b <= 4 and any(all((x, y) not in Es for x, y in itertools.combinations(S, 2))
                          for S in itertools.combinations(range(4), b)):
            continue
        out.append(h)
    return out


def _g6_matrix(line: bytes):
    import numpy as np

    line = line.strip()
    s = line[0] - 63
    bits = []
    for ch in line[1:]:
        v = ch - 63
        bits += [(v >> k) & 1 for k in range(5, -1, -1)]
    A = np.zeros((s, s), dtype=np.int8)
    pos = 0
    for j in range(1, s):
        for i in range(j):
            if bits[pos]:
                A[i, j] = A[j, i] = 1
            pos += 1
    return s, A


def quad_vectors(path: str, s: int) -> Tuple[int, Set[Tuple[int, ...]]]:
    """(number of graphs, distinct count vectors c) of the graph6 file (optionally .gz)."""
    import numpy as np

    lut = np.zeros(64, dtype=np.int64)
    for mask in range(64):
        lut[mask] = class_of([PAIRS[b] for b in range(6) if mask >> b & 1])
    quads = np.array(list(itertools.combinations(range(s), 4)), dtype=np.int64)
    seen, count = set(), 0
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rb") as fh:
        for line in fh:
            if not line.strip():
                continue
            n, A = _g6_matrix(line)
            if n != s:
                raise ValueError(f"{path}: a graph on {n} vertices, expected {s}")
            mask = np.zeros(len(quads), dtype=np.int64)
            for bit, (i, j) in enumerate(PAIRS):
                mask |= A[quads[:, i], quads[:, j]].astype(np.int64) << bit
            seen.add(tuple(int(x) for x in np.bincount(lut[mask], minlength=11)))
            count += 1
    return count, seen


def list_path(data_dir: str, a: int, b: int, s: int):
    base = os.path.join(data_dir, f"r{a}{b}_{s}.g6")
    for p in (base, base + ".gz"):
        if os.path.exists(p):
            return p
    return None


def min_dot(a: List[Fraction], vectors) -> Fraction:
    """min_c a.c exactly (integer arithmetic over the common denominator)."""
    D = 1
    for x in a:
        D = D * x.denominator // gcd(D, x.denominator)
    ai = [int(x * D) for x in a]
    return Fraction(min(sum(p * q for p, q in zip(ai, c)) for c in vectors), D)


def trivial_offset(a: List[Fraction], fam: Tuple[int, int], s: int) -> Fraction:
    return comb(s, 4) * min(a[h] for h in admissible_classes(*fam))


def family(k: int, l: int, N: int, colour: str, d: int) -> Tuple[Tuple[int, int], int]:
    """((a',b'), order) of the red / blue neighbourhood of a root of red degree d."""
    return ((k - 1, l), d) if colour == "red" else ((k, l - 1), N - 1 - d)


def load_cuts() -> List[dict]:
    """The stored cuts with Fraction directions and offsets {d: (Fraction, source)}."""
    with open(CUTS_FILE) as fh:
        raw = json.load(fh)
    out = []
    for c in raw["cuts"]:
        cut = {"k": c["k"], "l": c["l"], "N": c["N"], "colour": c["colour"],
               "a": [Fraction(x) for x in c["a"]],
               "offsets": {int(d): (Fraction(v["b"]), v["source"]) for d, v in c["offsets"].items()}}
        if "degree" in c:
            cut["degree"] = int(c["degree"])
        out.append(cut)
    return out


def check_cuts(cuts, data_dir: str, log=lambda s: None) -> None:
    """Recompute every offset: 'trivial' ones from the formula, 'list' ones from the lists
    in data_dir (raising if a list is absent, has the wrong count, or the offset differs)."""
    cache: Dict[Tuple[int, int, int], Set] = {}
    for c in cuts:
        for d, (b, source) in c["offsets"].items():
            fam, s = family(c["k"], c["l"], c["N"], c["colour"], d)
            if source == "trivial":
                want = trivial_offset(c["a"], fam, s)
            elif source == "list":
                key = (*fam, s)
                if key not in cache:
                    path = list_path(data_dir, *key)
                    if path is None:
                        raise FileNotFoundError(f"no list r{fam[0]}{fam[1]}_{s} in {data_dir}")
                    count, vecs = quad_vectors(path, s)
                    if count != LIST_COUNTS.get(key):
                        raise ValueError(f"{path}: {count} graphs, published {LIST_COUNTS.get(key)}")
                    cache[key] = vecs
                    log(f"r{fam[0]}{fam[1]}_{s}: {count} graphs, {len(vecs)} distinct vectors")
                want = min_dot(c["a"], cache[key])
            else:
                raise ValueError(f"unknown offset source {source!r}")
            if want != b:
                raise ValueError(f"cut {c['colour']} {c['a']} d={d}: stored {b}, recomputed {want}")


def make_cuts(directions, k, l, N, window, data_dir: str) -> List[dict]:
    """Offsets for the given [(colour, a)] or [(colour, a, d)] (a cut localised at red degree
    d) from the lists in data_dir (trivial where absent)."""
    out = []
    cache = {}
    for direction in directions:
        colour, a = direction[:2]
        degree = direction[2] if len(direction) > 2 else None
        a = [Fraction(x) for x in a]
        offs = {}
        for d in ([degree] if degree is not None else window):
            fam, s = family(k, l, N, colour, d)
            key = (*fam, s)
            path = list_path(data_dir, *key) if key in LIST_COUNTS else None
            if path:
                if key not in cache:
                    count, cache[key] = quad_vectors(path, s)
                    assert count == LIST_COUNTS[key], (path, count)
                offs[d] = {"b": str(min_dot(a, cache[key])), "source": "list"}
            else:
                offs[d] = {"b": str(trivial_offset(a, fam, s)), "source": "trivial"}
        cut = {"k": k, "l": l, "N": N, "colour": colour, "a": [str(x) for x in a], "offsets": offs}
        if degree is not None:
            cut["degree"] = degree
        out.append(cut)
    return out


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    if argv[:1] == ["check"]:
        check_cuts(load_cuts(), argv[1], log=print)
        print("OK: every stored offset recomputed")
    else:
        raise SystemExit("usage: python -m flagcheck.nbhd_quads check <dir with McKay's lists>")


if __name__ == "__main__":
    main()
