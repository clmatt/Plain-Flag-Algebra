"""The blow-up encoding of finite graphs (Lidicky-Pfender; thesis Ch. II.2) and Ramsey
feasibility problems.

Colours: 0 = green (two points in the same blow-up class), 1 = red (edge), 2 = blue (non-edge).

Soundness argument for every constraint below. Suppose G is a red/blue colouring of K_N
with no red K_k and no blue K_l. Its balanced blow-up (each vertex becomes a class of
measure 1/N, green inside classes) is a limit object that:
  * contains none of the forbidden configurations, so every inadmissible graph has density 0;
  * has p(green K_j) = N * (1/N)^j = 1/N^(j-1) for every j >= 1;
  * gives a root in class v a red 1-flag density of exactly deg_red(v)/N.
So any statement that holds pointwise for every root (e.g. prod_{d in D} (R - d/N) = 0 when
every red degree lies in D) holds after averaging. Every function here builds such a valid
statement; the tests check them exactly on blow-ups of real Ramsey graphs.
"""

from fractions import Fraction
from itertools import combinations
from typing import Iterable, List, Sequence

from .algebra import Combo, add, average, combo, identity_flag, lift, multiply, scale
from .graphs import Flag, canonical, from_color_function, from_edges

GREEN, RED, BLUE = 0, 1, 2
NUM_COLORS = 3


def blowup_zeros() -> List[Flag]:
    """The three forbidden 3-vertex configurations that make green an equivalence relation
    with homogeneous colours between classes."""
    return [
        from_edges(3, [(0, 1, RED)], default=GREEN),            # two green edges + red
        from_edges(3, [(0, 1, BLUE)], default=GREEN),           # two green edges + blue
        from_edges(3, [(0, 1, RED), (1, 2, BLUE)], default=GREEN),  # 0~2 green, 1 sees them differently
    ]


def monochromatic_clique(k: int, color: int) -> Flag:
    return from_edges(k, [], default=color)


def ramsey_zeros(k: int, l: int) -> List[Flag]:
    """No red K_k, no blue K_l, plus the blow-up structure."""
    return blowup_zeros() + [monochromatic_clique(k, RED), monochromatic_clique(l, BLUE)]


# ---------------------------------------------------------------------------
# Finite graphs
# ---------------------------------------------------------------------------

def circulant(N: int, distances: Iterable[int]) -> Flag:
    """Red/blue colouring of K_N: red iff the cyclic distance is in `distances`."""
    D = {d % N for d in distances} | {(-d) % N for d in distances}
    return from_color_function(N, lambda i, j: RED if (j - i) % N in D else BLUE)


def paley(q: int) -> Flag:
    squares = {(x * x) % q for x in range(1, q)}
    return from_color_function(q, lambda i, j: RED if (j - i) % q in squares else BLUE)


def has_mono_clique(G: Flag, k: int, color: int) -> bool:
    """Brute force with neighbourhood pruning; fine for N <= ~40, k <= 6."""
    adj = [{j for j in range(G.k) if j != i and G.color(i, j) == color} for i in range(G.k)]

    def extend(clique_size: int, candidates: set) -> bool:
        if clique_size == k:
            return True
        if clique_size + len(candidates) < k:
            return False
        for v in sorted(candidates):
            if extend(clique_size + 1, {u for u in candidates if u > v} & adj[v]):
                return True
        return False

    return extend(0, set(range(G.k)))


def is_ramsey_graph(G: Flag, k: int, l: int) -> bool:
    return not has_mono_clique(G, k, RED) and not has_mono_clique(G, l, BLUE)


def red_degrees(G: Flag) -> List[int]:
    return [sum(1 for j in range(G.k) if j != i and G.color(i, j) == RED) for i in range(G.k)]


# Known critical graphs (the tests check their Ramsey property directly).
KNOWN_RAMSEY_GRAPHS = {
    (3, 3): lambda: circulant(5, [1]),          # C5
    (3, 4): lambda: circulant(8, [1, 4]),       # Wagner graph
    (3, 5): lambda: circulant(13, [1, 5]),      # cyclic (3,5)-graph
    (4, 4): lambda: paley(17),                  # Paley(17)
}


# ---------------------------------------------------------------------------
# Constraint expressions (all valid for balanced blow-ups of graphs on exactly N vertices)
# ---------------------------------------------------------------------------

ROOT = Flag(1, 0, ())  # the type of single-vertex flags


def red_flag() -> Flag:
    """1-flag: root plus one vertex joined to it by a red edge. Density = deg_red(v)/N."""
    return from_edges(2, [(0, 1, RED)], s=1)


def one_flag() -> Combo:
    """The constant 1 as a 1-flag (the identity flag on the single-vertex type)."""
    return combo((identity_flag(ROOT), 1))


def green_clique_constraint(j: int, N: int):
    """p(green K_j) = 1/N^(j-1)."""
    return combo((monochromatic_clique(j, GREEN), 1)), Fraction(1, N ** (j - 1))


def degree_polynomial(allowed: Sequence[int], N: int, zeros, extra_factors: Sequence[int] = ()) -> Combo:
    """[[ prod_{d in allowed + extra_factors} (R - d/N) ]] as a combination of graphs.
    It equals 0 whenever every red degree lies in `allowed` (extra factors may repeat
    allowed values; the product is still 0 pointwise).

    Products only keep admissible flags (those avoiding `zeros`). That is exact for every
    limit object with zero density on the zeros, which is the only kind we use."""
    expr = one_flag()
    for d in list(allowed) + list(extra_factors):
        factor = add(combo((red_flag(), 1)), one_flag(), Fraction(-d, N))
        expr = multiply(expr, factor, NUM_COLORS, tuple(zeros))
    return average(expr)


def rooted_degree_polynomial(allowed: Sequence[int], N: int, zeros) -> Combo:
    """prod_{d in allowed} (R - d/N) as a 1-flag combination, before averaging. It is 0 at
    every root whose red degree lies in `allowed`."""
    expr = one_flag()
    for d in allowed:
        factor = add(combo((red_flag(), 1)), one_flag(), Fraction(-d, N))
        expr = multiply(expr, factor, NUM_COLORS, tuple(zeros))
    return expr


def localized_degree_constraints(allowed: Sequence[int], N: int, zeros, n: int):
    """[[ P(R) * F ]] = 0 for every admissible 1-flag F that fits in n vertices, where
    P(R) = prod_{d in allowed} (R - d/N).

    Valid because P(R) vanishes at every root of a real graph whose degrees lie in
    `allowed`, so P(R)*F vanishes pointwise too. F = the identity flag gives the plain
    averaged constraint `degree_polynomial`."""
    from .graphs import generate_flags

    zeros = tuple(zeros)
    P = rooted_degree_polynomial(allowed, N, zeros)
    out = []
    for m in range(1, n - len(allowed) + 1):
        for F in generate_flags(ROOT, m, NUM_COLORS, zeros):
            expr = average(multiply(P, combo((F, 1)), NUM_COLORS, zeros))
            if expr:
                out.append((expr, Fraction(0)))
    return out


def lift_to(expr: Combo, n: int, zeros) -> Combo:
    return lift(expr, n, NUM_COLORS, zeros)
