"""Pointwise ("localised") constraints from the structure of neighbourhoods.

Let H be a (k,l)-graph on N vertices (no red K_k, no blue K_l) and W its balanced blow-up.
For a root tuple theta (one point in each of distinct classes, or the same class for
green pairs) every rooted density below is an exact function of H:

  green    p(x green to root 0)                = 1/N              (classes have measure 1/N)
  codeg    p(x has colours (c0, c1) to roots)  = |S| / N          (S = vertices of H with
                                                                   those colours, outside
                                                                   the two root classes)
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
from itertools import product
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
