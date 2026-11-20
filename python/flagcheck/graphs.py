"""Complete edge-coloured graphs and flags, with brute-force canonical forms.

This is an independent reference implementation. It favours being obviously correct over
being fast, so it is only meant for small graphs (roughly up to 7 vertices).

A flag is a triple (k, s, cols):
  k    -- number of vertices, labelled 0..k-1
  s    -- number of root (labelled) vertices; the roots are vertices 0..s-1, in order
  cols -- tuple of edge colours for the pairs (i, j), i < j, in lexicographic order

A plain graph is a flag with s = 0. Two flags are isomorphic iff some permutation that
fixes every root maps one onto the other (the type is labelled, as in the thesis).
"""

from functools import lru_cache
from itertools import combinations, permutations, product
from typing import Iterable, NamedTuple, Sequence


class Flag(NamedTuple):
    k: int
    s: int
    cols: tuple

    def color(self, i: int, j: int) -> int:
        if i == j:
            raise ValueError("no colour on a loop")
        if i > j:
            i, j = j, i
        return self.cols[pair_index(self.k)[i][j]]

    def type(self) -> "Flag":
        """The labelled type: the subgraph induced on the roots."""
        return restrict(self, range(self.s), self.s)

    def underlying(self) -> "Flag":
        return canonical(Flag(self.k, 0, self.cols))


@lru_cache(maxsize=None)
def pairs(k: int) -> tuple:
    return tuple(combinations(range(k), 2))


@lru_cache(maxsize=None)
def pair_index(k: int) -> tuple:
    idx = [[-1] * k for _ in range(k)]
    for p, (i, j) in enumerate(pairs(k)):
        idx[i][j] = p
        idx[j][i] = p
    return tuple(tuple(row) for row in idx)


def from_edges(k: int, edges: Iterable, default: int = 0, s: int = 0) -> Flag:
    """Build a flag from (i, j, colour) triples; every other pair gets `default`."""
    idx = pair_index(k)
    cols = [default] * len(pairs(k))
    for i, j, c in edges:
        if i == j or not (0 <= i < k and 0 <= j < k):
            raise ValueError(f"bad edge {(i, j)} for k={k}")
        cols[idx[i][j]] = c
    return Flag(k, s, tuple(cols))


def from_color_function(k: int, color, s: int = 0) -> Flag:
    return Flag(k, s, tuple(color(i, j) for i, j in pairs(k)))


def restrict(F: Flag, vertices: Sequence[int], s: int) -> Flag:
    """The flag induced on `vertices`, in the given order; the first `s` become the roots."""
    vertices = list(vertices)
    m = len(vertices)
    idx = pair_index(F.k)
    cols = tuple(F.cols[idx[vertices[a]][vertices[b]]] for a, b in pairs(m))
    return Flag(m, s, cols)


@lru_cache(maxsize=None)
def _relabel_maps(k: int, s: int) -> tuple:
    """For every permutation fixing the roots, the old pair position for each new pair."""
    idx = pair_index(k)
    maps = []
    for tail in permutations(range(s, k)):
        perm = tuple(range(s)) + tail  # new vertex x is old vertex perm[x]
        maps.append(tuple(idx[perm[a]][perm[b]] for a, b in pairs(k)))
    return tuple(maps)


@lru_cache(maxsize=None)
def _canonical_cols(k: int, s: int, cols: tuple) -> tuple:
    return min(tuple(cols[p] for p in m) for m in _relabel_maps(k, s))


def canonical(F: Flag) -> Flag:
    """Canonical representative: the lexicographically smallest colour tuple over all
    relabellings that fix the roots."""
    return Flag(F.k, F.s, _canonical_cols(F.k, F.s, F.cols))


def isomorphic(F: Flag, G: Flag) -> bool:
    return F.k == G.k and F.s == G.s and canonical(F) == canonical(G)


def automorphism_count(F: Flag) -> int:
    return sum(1 for m in _relabel_maps(F.k, F.s) if tuple(F.cols[p] for p in m) == F.cols)


# ---------------------------------------------------------------------------
# Forbidden subgraphs ("zeros") and generation
# ---------------------------------------------------------------------------

def contains_induced(G: Flag, H: Flag) -> bool:
    """Does the graph G (roots ignored) contain an induced copy of the graph H?"""
    target = canonical(Flag(H.k, 0, H.cols))
    for X in combinations(range(G.k), H.k):
        if canonical(restrict(G, X, 0)) == target:
            return True
    return False


def _admissible_with_new_vertex(G: Flag, zeros_by_size: dict) -> bool:
    """Check only the vertex subsets that contain the last vertex. When G minus its last
    vertex is already admissible, this decides whether G is admissible."""
    last = G.k - 1
    for size, targets in zeros_by_size.items():
        if size > G.k:
            continue
        for rest in combinations(range(last), size - 1):
            if canonical(restrict(G, rest + (last,), 0)).cols in targets:
                return False
    return True


def _zeros_by_size(zeros: Iterable) -> dict:
    out = {}
    for Z in zeros:
        out.setdefault(Z.k, set()).add(canonical(Flag(Z.k, 0, Z.cols)).cols)
    return out


def is_admissible(G: Flag, zeros: Iterable) -> bool:
    return not any(contains_induced(G, Z) for Z in zeros)


def extend_by_vertex(F: Flag, colors: Sequence[int]) -> Flag:
    """Add vertex k, with colour colors[i] on the pair (i, k)."""
    k = F.k
    new = dict(zip(pairs(k), F.cols))
    for i in range(k):
        new[(i, k)] = colors[i]
    return Flag(k + 1, F.s, tuple(new[p] for p in pairs(k + 1)))


@lru_cache(maxsize=None)
def _generate(n: int, c: int, zeros_key: tuple) -> tuple:
    zeros_by_size = {size: set(t) for size, t in zeros_key}
    if n == 0:
        return (Flag(0, 0, ()),)
    out = set()
    for G in _generate(n - 1, c, zeros_key):
        for colors in product(range(c), repeat=n - 1):
            H = extend_by_vertex(G, colors)
            if _admissible_with_new_vertex(H, zeros_by_size):
                out.add(canonical(H))
    return tuple(sorted(out))


def _zeros_key(zeros: Iterable) -> tuple:
    return tuple(sorted((size, tuple(sorted(t))) for size, t in _zeros_by_size(zeros).items()))


def generate(n: int, c: int, zeros: Iterable = ()) -> tuple:
    """All admissible complete c-edge-coloured graphs on n vertices, up to isomorphism
    (canonical representatives, sorted)."""
    return _generate(n, c, _zeros_key(zeros))


@lru_cache(maxsize=None)
def _generate_flags(sigma: Flag, m: int, c: int, zeros_key: tuple) -> tuple:
    zeros_by_size = {size: set(t) for size, t in zeros_key}
    if m == sigma.k:
        return (Flag(sigma.k, sigma.k, sigma.cols),)
    out = set()
    for F in _generate_flags(sigma, m - 1, c, zeros_key):
        for colors in product(range(c), repeat=m - 1):
            H = extend_by_vertex(F, colors)
            if _admissible_with_new_vertex(H, zeros_by_size):
                out.add(canonical(H))
    return tuple(sorted(out))


def generate_flags(sigma: Flag, m: int, c: int, zeros: Iterable = ()) -> tuple:
    """All admissible sigma-flags on m vertices (sigma is a labelled type), up to
    flag isomorphism."""
    sigma = Flag(sigma.k, 0, sigma.cols)
    if m < sigma.k:
        raise ValueError("flag smaller than its type")
    return _generate_flags(sigma, m, c, _zeros_key(zeros))
