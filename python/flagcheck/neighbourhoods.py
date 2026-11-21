"""Pointwise ("localised") constraints from the structure of neighbourhoods.

Let H be a (k,l)-graph on N vertices (no red K_k, no blue K_l) and W its balanced blow-up.
For a root tuple theta (one point in each of distinct classes, or the same class for
green pairs) every rooted density below is an exact function of H:

  green    p(x green to root 0)                = 1/N              (classes have measure 1/N)
  codeg    p(x has colours (c0, c1) to roots)  = |S| / N          (S = vertices of H with
                                                                   those colours, outside
                                                                   the two root classes)
  tcodeg   the same with three roots
  edges    p(x, y red to root, xy red)         = 2 e(N_r(v)) / N^2
           p(x, y blue to root, xy red)        = 2 e_red(N_b(v)) / N^2

A statement g >= 0 (or g = 0) that holds at EVERY root gives, for every flag F of the same
type, [[g F]] >= 0 (or = 0) after averaging, since F >= 0. These are the linear constraints
built here; `localizing_specs` gives the PSD versions [[g v v^T]] >= 0.

Hypotheses used (all recorded in Problem.assumptions):
  * S above induces a (k - a, l - b)-graph, where a (b) is the size of the largest red
    (blue) clique among the roots joined to S in red (blue); so |S| <= R(k-a, l-b) - 1.
  * N_r(v) is a (k-1, l)-graph on d(v) vertices and N_b(v) a (k, l-1)-graph on N-1-d(v),
    so their red edge counts lie in the ranges of edge_tables.edge_bounds.
  * Every red degree lies in the window [N - R(k,l-1), R(k-1,l) - 1].
"""

from fractions import Fraction
from itertools import combinations, permutations, product
from math import comb
from typing import List, Sequence, Tuple

from .algebra import identity_flag
from .edge_tables import edge_bounds
from .graphs import Flag, from_edges, generate_flags
from .ramsey import BLUE, GREEN, NUM_COLORS, RED, ROOT
from .spec import Add, Avg, Lit, Mul, Scale

RED_PAIR = from_edges(2, [(0, 1, RED)])     # type: two roots joined in red
BLUE_PAIR = from_edges(2, [(0, 1, BLUE)])   # type: two roots joined in blue
COLOR_NAME = {GREEN: "g", RED: "r", BLUE: "b"}


def one(sigma: Flag):
    return Lit(identity_flag(sigma))


def extension(sigma: Flag, to_roots: Sequence[int]) -> Flag:
    """sigma plus one vertex with the given colours to the roots."""
    s = sigma.k
    edges = [(i, j, sigma.color(i, j)) for i in range(s) for j in range(i + 1, s)]
    edges += [(i, s, c) for i, c in enumerate(to_roots)]
    return from_edges(s + 1, edges, s=s)


def double_extension(to_root: int, between: int) -> Flag:
    """1-flag: root plus x, y, both joined to the root in `to_root`, xy in `between`."""
    return from_edges(3, [(0, 1, to_root), (0, 2, to_root), (1, 2, between)], s=1)


def expr_size(e) -> int:
    """Number of vertices of the largest flag in a product expansion of e."""
    if isinstance(e, Lit):
        return e.flag.k
    if isinstance(e, dict):
        return max(F.k for F in e)
    if isinstance(e, Add):
        return max(expr_size(a) for a in e.args)
    if isinstance(e, Scale):
        return expr_size(e.arg)
    if isinstance(e, Mul):
        sizes = [expr_size(a) for a in e.args]
        s = sigma_of(e).k
        return sum(sizes) - s * (len(sizes) - 1)
    raise TypeError(e)


# ---------------------------------------------------------------------------
# Multipliers g (each valid pointwise, see the module docstring)
# ---------------------------------------------------------------------------

def green_multipliers(N: int):
    """g = p(x green to root 0) - 1/N = 0, on the 1-vertex type and on both 2-vertex types
    (on a pair type the vertex green to root 0 has root 0's colour to root 1)."""
    out = []
    for sigma in (ROOT, RED_PAIR, BLUE_PAIR):
        to_roots = [GREEN] + [sigma.color(0, j) for j in range(1, sigma.k)]
        g = Add((Lit(extension(sigma, to_roots)), Scale(Fraction(-1, N), one(sigma))))
        out.append((f"green[{sigma.k}]", g))
    return out


def codegree_family(k: int, l: int, root_color: int, to_roots: Tuple[int, int]) -> Tuple[int, int]:
    """(k - a, l - b): S = vertices with colours to_roots to a pair of roots joined in
    root_color induces a (k - a, l - b)-graph."""
    reds = [i for i, c in enumerate(to_roots) if c == RED]
    blues = [i for i, c in enumerate(to_roots) if c == BLUE]
    a = len(reds) if len(reds) < 2 or root_color == RED else 1
    b = len(blues) if len(blues) < 2 or root_color == BLUE else 1
    return k - a, l - b


def codegree_bound(k: int, l: int, root_color: int, to_roots: Tuple[int, int], ramsey_upper) -> int:
    """Upper bound on |S| (see codegree_family)."""
    return ramsey_upper(*codegree_family(k, l, root_color, to_roots)) - 1


def codegree_multipliers(k: int, l: int, N: int, ramsey_upper):
    """g = c/N - p(x has colours (c0, c1) to the roots) >= 0 for both pair types and the
    colour patterns rr, rb, bb (br is rb with the roots swapped; averaging over all root
    labellings makes the two families equivalent)."""
    out, facts = [], []
    for sigma in (RED_PAIR, BLUE_PAIR):
        rc = sigma.color(0, 1)
        for pattern in ((RED, RED), (RED, BLUE), (BLUE, BLUE)):
            c = codegree_bound(k, l, rc, pattern, ramsey_upper)
            if c >= N - 2:
                continue                       # trivially true
            name = f"codeg[{COLOR_NAME[rc]}:{COLOR_NAME[pattern[0]]}{COLOR_NAME[pattern[1]]}<={c}]"
            g = Add((Scale(Fraction(c, N), one(sigma)), Scale(Fraction(-1), Lit(extension(sigma, pattern)))))
            out.append((name, g))
            ka, lb = codegree_family(k, l, rc, pattern)
            facts.append(f"R({ka},{lb}) <= {c + 1}")
    return out, facts


def triple_type(reds: int) -> Flag:
    """3-vertex type with `reds` red edges (first) and the rest blue."""
    cols = [RED] * reds + [BLUE] * (3 - reds)
    return from_edges(3, [(0, 1, cols[0]), (0, 2, cols[1]), (1, 2, cols[2])])


def _max_clique(sigma: Flag, verts: Sequence[int], color: int) -> int:
    return max(r for r in range(len(verts) + 1) for sub in combinations(verts, r)
               if all(sigma.color(i, j) == color for i, j in combinations(sub, 2)))


def root_family(k: int, l: int, sigma: Flag, to_roots: Sequence[int]) -> Tuple[int, int]:
    """(k - a, l - b) for any root type: a (b) is the largest red (blue) clique of sigma
    among the roots joined to S in red (blue). Agrees with codegree_family on pairs."""
    reds = [i for i, c in enumerate(to_roots) if c == RED]
    blues = [i for i, c in enumerate(to_roots) if c == BLUE]
    return k - _max_clique(sigma, reds, RED), l - _max_clique(sigma, blues, BLUE)


def triple_codegree_multipliers(k: int, l: int, N: int, ramsey_upper):
    """g = c/N - p(x has colours P to the three roots) >= 0, c = R(root_family) - 1, for
    every admissible red/blue triangle type and pattern P up to automorphisms of the type.
    Skipped when trivially true (c >= N - 3), when the extension is itself forbidden
    (c = 0), or when implied by a pair bound (c >= the bound for some two of the roots:
    S here is a subset of that pair's S)."""
    out, facts = [], []
    for reds in range(4):
        sigma = triple_type(reds)
        if _max_clique(sigma, range(3), RED) >= k or _max_clique(sigma, range(3), BLUE) >= l:
            continue
        auts = [p for p in permutations(range(3))
                if all(sigma.color(p[i], p[j]) == sigma.color(i, j) for i, j in combinations(range(3), 2))]
        seen = set()
        for pattern in product((RED, BLUE), repeat=3):
            key = min(tuple(pattern[p[i]] for i in range(3)) for p in auts)
            if key in seen:
                continue
            seen.add(key)
            ka, lb = root_family(k, l, sigma, pattern)
            c = ramsey_upper(ka, lb) - 1
            pair_c = min(ramsey_upper(*root_family(k, l, from_edges(2, [(0, 1, sigma.color(i, j))]),
                                                   (pattern[i], pattern[j]))) - 1
                         for i, j in combinations(range(3), 2))
            if c == 0 or c >= N - 3 or c >= pair_c:
                continue
            tag = "".join(COLOR_NAME[sigma.color(i, j)] for i, j in ((0, 1), (0, 2), (1, 2)))
            name = f"tcodeg[{tag}:{''.join(COLOR_NAME[x] for x in pattern)}<={c}]"
            g = Add((Scale(Fraction(c, N), one(sigma)), Scale(Fraction(-1), Lit(extension(sigma, pattern)))))
            out.append((name, g))
            facts.append(f"R({ka},{lb}) <= {c + 1}")
    return out, facts


def _lagrange(values, window, N):
    """The polynomial in R (red degree / N) that equals values[d] at R = d/N for d in window."""
    from .ramsey import one_flag, red_flag

    def red_minus(c):
        return Add((Lit(red_flag()), Scale(-c, one_flag())))

    terms = []
    for d in window:
        others = [e for e in window if e != d]
        if not others:
            terms.append(Scale(values[d], one_flag()))
            continue
        denom = Fraction(1)
        for e in others:
            denom *= Fraction(d - e, N)
        terms.append(Scale(values[d] / denom, Mul(tuple(red_minus(Fraction(e, N)) for e in others))))
    return Add(tuple(terms))


def _hull_lines(points, lower: bool):
    """Lines a*x + b with a*x + b <= y (lower) or >= y (upper) at every point, through
    consecutive vertices of the lower (upper) convex hull."""
    pts = sorted(points)
    hull = []
    for p in pts:
        while len(hull) >= 2:
            (x1, y1), (x2, y2) = hull[-2], hull[-1]
            cross = (x2 - x1) * (p[1] - y1) - (y2 - y1) * (p[0] - x1)
            if (cross <= 0) if lower else (cross >= 0):
                hull.pop()
            else:
                break
        hull.append(p)
    lines = []
    for (x1, y1), (x2, y2) in zip(hull, hull[1:]):
        a = Fraction(y2 - y1, x2 - x1)
        lines.append((a, y1 - a * x1))
    if len(hull) == 1:
        lines.append((Fraction(0), Fraction(hull[0][1])))
    for a, b in lines:              # re-check validity exactly
        for x, y in points:
            assert (a * x + b <= y) if lower else (a * x + b >= y)
    return lines


def edge_multipliers(k: int, l: int, N: int, window: Sequence[int], n: int):
    """Red edge counts inside the red and the blue neighbourhood of every vertex.

    For each of the two neighbourhoods and each side (min / max) this gives
      * the exact interpolant g = T - L(R) >= 0 (or L(R) - T >= 0) if it fits in n vertices
        (L agrees with the table at every degree in the window, so g >= 0 at every root);
      * the hull lines g = T - (a R + b) >= 0 (degree 1, always fits)."""
    from .ramsey import one_flag, red_flag

    out, facts = [], []
    window = list(window)
    for label, T, (a, b), size in (
        ("red-nbhd", double_extension(RED, RED), (k - 1, l), lambda d: d),
        ("blue-nbhd", double_extension(BLUE, RED), (k, l - 1), lambda d: N - 1 - d),
    ):
        bounds = {}
        for d in window:
            try:
                bounds[d] = edge_bounds(a, b, size(d))
            except KeyError:
                raise ValueError(f"no edge data for ({a},{b})-graphs on {size(d)} vertices") from None
        facts += [f"e({a},{b},{size(d)}) in [{lo},{hi}]" for d, (lo, hi) in bounds.items()]
        for side, idx in (("min", 0), ("max", 1)):
            sign = 1 if side == "min" else -1          # g = sign * (T - bound)
            vals = {d: Fraction(2 * bounds[d][idx], N * N) for d in window}
            # exact interpolant
            L = _lagrange(vals, window, N)
            g = Add((Scale(sign, Lit(T)), Scale(-sign, L)))
            if expr_size(g) <= n:
                out.append((f"{label}-{side}-exact", g))
            # hull lines in x = R = d/N
            pts = [(Fraction(d, N), vals[d]) for d in window]
            for i, (sl, ic) in enumerate(_hull_lines(pts, lower=(side == "min"))):
                line = Add((Scale(sl, Lit(red_flag())), Scale(ic, one_flag())))
                out.append((f"{label}-{side}-line{i}", Add((Scale(sign, Lit(T)), Scale(-sign, line)))))
    return out, facts


def triple_extension(to_root: int, reds: int) -> Flag:
    """1-flag: root plus x, y, z, all joined to the root in `to_root`, with `reds` red
    edges among them (the rest blue); reds = 2 is the path x-y-z... up to isomorphism."""
    among = [(1, 2), (1, 3), (2, 3)]
    edges = [(0, x, to_root) for x in (1, 2, 3)]
    edges += [(i, j, RED if t < reds else BLUE) for t, (i, j) in enumerate(among)]
    return from_edges(4, edges, s=1)


def _plane_lines(plane, lower: bool):
    """Valid lines y >= a x + b (lower) or y <= a x + b (upper) through the hull of the points;
    a single x value gives the horizontal line through the extreme y."""
    best = {}                     # per x only the extreme y matters (and avoids vertical edges)
    for x, y in plane:
        best[x] = y if x not in best else (min(best[x], y) if lower else max(best[x], y))
    if len(best) == 1:
        (y,) = best.values()
        return [(Fraction(0), y)]
    lines = _hull_lines(list(best.items()), lower)
    for a, b in lines:            # valid for every point, not only the per-x extremes
        assert all((a * x + b <= y) if lower else (a * x + b >= y) for x, y in plane)
    return lines


def stat_multipliers(k: int, l: int, N: int, window: Sequence[int]):
    """Neighbourhood structure from complete lists (nbhd_stats.TRIPLE_POINTS).

    At a root v of degree d, N_r(v) is a (k-1, l)-graph on d vertices and N_b(v) a
    (k, l-1)-graph on N-1-d; when every such graph is listed, (e, p, t) of the
    neighbourhood is one of the listed points. In densities at the root,
        E = 2e/N^2 (double_extension),  P = 6p/N^3, T = 6t/N^3 (triple_extension),
    and each exact convex-hull line Y >= a_d X + b_d (or <=) of the listed points in the
    planes (E,P), (E,T), (P,T) holds at roots of degree d. As for the edge constraints, a_d
    and b_d are interpolated exactly in R over the window: g = +-(Y - A(R) X - B(R)) with
    A(d/N) = a_d, B(d/N) = b_d, so g >= 0 at every root. Degrees without a complete list
    get the trivially valid lines Y >= 0 and Y <= 6 C(s,3) / N^3."""
    from .nbhd_stats import table, triple_points
    from .ramsey import one_flag

    out, facts = [], []
    window = list(window)
    for label, rc, (a, b), size in (("red", RED, (k - 1, l), lambda d: d),
                                    ("blue", BLUE, (k, l - 1), lambda d: N - 1 - d)):
        if a <= 2 or b <= 2:
            continue                      # no red edges, or a red clique: nothing to say
        tabled = {d: triple_points(a, b, size(d)) for d in window if size(d) >= 3}
        tabled = {d: pts for d, pts in tabled.items() if pts}
        if not tabled:
            continue
        key = (min(a, b), max(a, b))
        facts += [f"triple stats of all {size(d)}-vertex ({key[0]},{key[1]})-graphs "
                  f"({table()[key][size(d)]['count']} graphs, McKay)" for d in tabled]
        coords = {"E": (lambda q: Fraction(2 * q[0], N * N), Lit(double_extension(rc, RED))),
                  "P": (lambda q: Fraction(6 * q[1], N ** 3), Lit(triple_extension(rc, 2))),
                  "T": (lambda q: Fraction(6 * q[2], N ** 3), Lit(triple_extension(rc, 3)))}
        planes = [("E", "P")] + ([("E", "T"), ("P", "T")] if a >= 4 else [])   # T = 0 when a = 3
        for xn, yn in planes:
            (fx, X), (fy, Y) = coords[xn], coords[yn]
            for side in ("min", "max"):
                lower = side == "min"
                sign = 1 if lower else -1               # g = sign * (Y - A X - B)
                lines = {}
                for d in window:
                    if d in tabled:
                        lines[d] = _plane_lines({(fx(q), fy(q)) for q in tabled[d]}, lower)
                    else:
                        lines[d] = [(Fraction(0), Fraction(0) if lower else Fraction(6 * comb(size(d), 3), N ** 3))]
                for i in range(max(len(v) for v in lines.values())):
                    A = _lagrange({d: v[i % len(v)][0] for d, v in lines.items()}, window, N)
                    B = _lagrange({d: v[i % len(v)][1] for d, v in lines.items()}, window, N)
                    g = Add((Scale(sign, Y), Scale(-sign, Mul((A, X))), Scale(-sign, B)))
                    out.append((f"nstats[{label}:{xn}{yn}]-{side}-{i}", g))
    return out, facts


def pair_extension(sigma: Flag, to_roots: Sequence[int], between: int) -> Flag:
    """2-flag: the pair type plus x, y, each with colours to_roots to the roots, xy in `between`."""
    rc = sigma.color(0, 1)
    edges = [(0, 1, rc), (2, 3, between)]
    for x in (2, 3):
        edges += [(0, x, to_roots[0]), (1, x, to_roots[1])]
    return from_edges(4, edges, s=2)


def pair_edge_multipliers(k: int, l: int, N: int, ramsey_upper):
    """Second level: at a pair of roots, S = vertices with a fixed colour pattern to the
    pair induces an (a', b')-graph (codegree_family), so with s = |S| the red edge count
    e(S) lies in [min, max] of edge_bounds(a', b', s). In densities,
        D = p(x has the pattern) = s/N,  T = p(x, y have the pattern, xy red) = 2 e(S)/N^2.
    s ranges over 0..R(a',b')-1 (too wide for an interpolant), so the constraints are the
    lower / upper convex-hull lines of the points (s/N, 2 e/N^2): T - (a D + b) >= 0 (or <=).
    Families with a' <= 2 or b' <= 2 are skipped (S has no red edge, or is a red clique;
    both are already forced by the forbidden cliques)."""
    out, facts = [], []
    for sigma in (RED_PAIR, BLUE_PAIR):
        rc = sigma.color(0, 1)
        for pattern in ((RED, RED), (RED, BLUE), (BLUE, BLUE)):
            a2, b2 = codegree_family(k, l, rc, pattern)
            if a2 <= 2 or b2 <= 2:
                continue
            top = ramsey_upper(a2, b2) - 1
            bounds = {s: edge_bounds(a2, b2, s) for s in range(0, top + 1)}
            facts.append(f"R({a2},{b2}) <= {top + 1}")
            facts += [f"e({a2},{b2},{s}) in [{lo},{hi}]" for s, (lo, hi) in bounds.items() if s]
            D = Lit(extension(sigma, pattern))
            T = Lit(pair_extension(sigma, pattern, RED))
            tag = f"{COLOR_NAME[rc]}:{COLOR_NAME[pattern[0]]}{COLOR_NAME[pattern[1]]}"
            for side, idx in (("min", 0), ("max", 1)):
                sign = 1 if side == "min" else -1
                pts = [(Fraction(s, N), Fraction(2 * bounds[s][idx], N * N)) for s in bounds]
                for i, (sl, ic) in enumerate(_hull_lines(pts, lower=(side == "min"))):
                    line = Add((Scale(sl, D), Scale(ic, one(sigma))))
                    out.append((f"pair-edges[{tag}]-{side}-line{i}", Add((Scale(sign, T), Scale(-sign, line)))))
    return out, facts


# ---------------------------------------------------------------------------
# Turning multipliers into constraints
# ---------------------------------------------------------------------------

def sigma_of(g) -> Flag:
    e = g
    while not isinstance(e, (Lit, dict)):
        e = e.args[0] if isinstance(e, (Add, Mul)) else e.arg
    t = (e.flag if isinstance(e, Lit) else next(iter(e))).type()
    return Flag(t.k, 0, t.cols)


def localized(g, max_size: int, zeros) -> List:
    """[[g F]] for every admissible flag F of g's type with |g| + |F| - s <= max_size."""
    sigma = sigma_of(g)
    s = sigma.k
    size = expr_size(g)
    out = []
    for m in range(s, max_size - size + s + 1):
        for F in generate_flags(sigma, m, NUM_COLORS, zeros):
            out.append(Avg(Mul((g, Lit(F)))) if m > s else Avg(g))
    return out


def localizing_specs(g, n: int):
    """(m, g) for the largest m with |g| + 2(m - s) <= n, or None."""
    s = sigma_of(g).k
    m = (n - expr_size(g)) // 2 + s
    return (m, g) if m > s else None
