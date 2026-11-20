"""Flag algebra operations written directly from the definitions in thesis Ch. I.

A linear combination of flags is a dict {canonical Flag: Fraction}. Every flag in one
combination must have the same labelled type, i.e. the same s and the same colours on
the roots.
"""

from fractions import Fraction
from itertools import combinations, permutations
from math import comb, perm
from typing import Dict, Iterable

from .graphs import Flag, canonical, generate, generate_flags, restrict

Combo = Dict[Flag, Fraction]


def combo(*terms) -> Combo:
    """combo((F, 1), (G, Fraction(1, 2)), ...) builds a canonicalised combination."""
    out: Combo = {}
    for F, coeff in terms:
        add_term(out, canonical(F), Fraction(coeff))
    return out


def add_term(out: Combo, F: Flag, coeff: Fraction) -> None:
    if coeff == 0:
        return
    v = out.get(F, Fraction(0)) + coeff
    if v == 0:
        out.pop(F, None)
    else:
        out[F] = v


def add(a: Combo, b: Combo, scale_b: Fraction = Fraction(1)) -> Combo:
    out = dict(a)
    for F, v in b.items():
        add_term(out, F, v * scale_b)
    return out


def scale(a: Combo, k) -> Combo:
    k = Fraction(k)
    return {F: v * k for F, v in a.items()} if k != 0 else {}


def same_type(F: Flag, G: Flag) -> bool:
    return F.s == G.s and F.type().cols == G.type().cols


def identity_flag(sigma: Flag) -> Flag:
    """The sigma-flag whose underlying graph is sigma itself. It is identically 1."""
    return Flag(sigma.k, sigma.k, sigma.cols)


# ---------------------------------------------------------------------------
# Densities between flags (exact, finite)
# ---------------------------------------------------------------------------

def flag_density(F: Flag, G: Flag) -> Fraction:
    """p(F, G): the probability that a random (|F|-s)-subset X of the non-roots of G gives
    (G[roots + X], roots) isomorphic to F. Both must have the same labelled type."""
    if not same_type(F, G):
        raise ValueError("flag_density needs flags of the same labelled type")
    s = F.s
    if F.k > G.k:
        raise ValueError("F larger than G")
    target = canonical(F)
    roots = tuple(range(s))
    hits = 0
    total = 0
    for X in combinations(range(s, G.k), F.k - s):
        total += 1
        if canonical(restrict(G, roots + X, s)) == target:
            hits += 1
    return Fraction(hits, total)


# ---------------------------------------------------------------------------
# Chain rule, product and averaging
# ---------------------------------------------------------------------------

def lift(a: Combo, m: int, c: int, zeros: Iterable = ()) -> Combo:
    """Chain rule (Thm 2): rewrite every flag as a combination of flags on m vertices.
    Flags with an inadmissible underlying graph have density 0 and are dropped."""
    zeros = tuple(zeros)
    out: Combo = {}
    for F, coeff in a.items():
        if F.k > m:
            raise ValueError("cannot lift to fewer vertices")
        for H in generate_flags(F.type(), m, c, zeros):
            d = flag_density(F, H)
            if d:
                add_term(out, H, coeff * d)
    return out


def product_coefficient(F1: Flag, F2: Flag, H: Flag) -> Fraction:
    """c(F1, F2, H) from Thm 3: the probability that a random split of the non-roots of H
    into parts of sizes |F1|-s and |F2|-s induces F1 on the first part and F2 on the second."""
    s = F1.s
    a, b = F1.k - s, F2.k - s
    assert H.k == s + a + b
    roots = tuple(range(s))
    t1, t2 = canonical(F1), canonical(F2)
    nonroots = range(s, H.k)
    hits = 0
    total = 0
    for X in combinations(nonroots, a):
        Z = tuple(v for v in nonroots if v not in X)
        total += 1
        if canonical(restrict(H, roots + X, s)) == t1 and canonical(restrict(H, roots + Z, s)) == t2:
            hits += 1
    return Fraction(hits, total)


def multiply(a: Combo, b: Combo, c: int, zeros: Iterable = ()) -> Combo:
    zeros = tuple(zeros)
    out: Combo = {}
    for F1, v1 in a.items():
        for F2, v2 in b.items():
            if not same_type(F1, F2):
                raise ValueError("can only multiply flags of the same labelled type")
            m = F1.k + F2.k - F1.s
            for H in generate_flags(F1.type(), m, c, zeros):
                coeff = product_coefficient(F1, F2, H)
                if coeff:
                    add_term(out, H, v1 * v2 * coeff)
    return out


def average_coefficient(F: Flag) -> Fraction:
    """q_sigma(F) from Thm 4: the probability that a uniformly random injective
    theta: [s] -> V(F) makes (F, theta) a flag isomorphic to F."""
    s = F.s
    target = canonical(F)
    hits = 0
    for theta in permutations(range(F.k), s):
        rest = tuple(v for v in range(F.k) if v not in theta)
        if canonical(restrict(F, theta + rest, s)) == target:
            hits += 1
    return Fraction(hits, perm(F.k, s))


def average(a: Combo) -> Combo:
    """The averaging (unlabelling) operator [[.]]: sigma-flags -> graphs."""
    out: Combo = {}
    for F, v in a.items():
        add_term(out, F.underlying(), v * average_coefficient(F))
    return out


# ---------------------------------------------------------------------------
# The SDP data Y_{G, sigma} (Algorithm 2, as probabilities)
# ---------------------------------------------------------------------------

def Y_matrix(G: Flag, sigma: Flag, flags: list) -> list:
    """Y[c][d] = probability that a uniformly random injective theta: [s] -> V(G),
    followed by a uniformly random split of the remaining n - s vertices into X and Z of
    size m = (n - s)/2 each, has (G[theta + X], theta) ~= flags[c] and
    (G[theta + Z], theta) ~= flags[d]. Pairs where G[theta] is not exactly the labelled
    type sigma contribute 0.

    Then [[ v^T A v ]] = sum_G <Y_G, A> G for symmetric A, where v is the vector `flags`.
    """
    n, s = G.k, sigma.k
    if (n - s) % 2:
        raise ValueError("n - |sigma| must be even")
    m = (n - s) // 2
    index = {canonical(F): i for i, F in enumerate(flags)}
    if len(index) != len(flags):
        raise ValueError("duplicate flags")
    size = len(flags)
    counts = [[0] * size for _ in range(size)]
    total = perm(n, s) * comb(n - s, m)
    for theta in permutations(range(n), s):
        if restrict(G, theta, 0).cols != sigma.cols:
            continue
        rest = [v for v in range(n) if v not in theta]
        for X in combinations(rest, m):
            Z = tuple(v for v in rest if v not in X)
            i = index[canonical(restrict(G, theta + X, s))]
            j = index[canonical(restrict(G, theta + Z, s))]
            counts[i][j] += 1
    return [[Fraction(x, total) for x in row] for row in counts]


def sdp_blocks(n: int, c: int, zeros: Iterable = (), include_empty_type: bool = False) -> list:
    """The standard maximal choice: for each type size s < n with s = n mod 2 and each
    labelled type sigma on s vertices (one representative per isomorphism class), the
    vector of all admissible sigma-flags on (n + s)/2 vertices.

    The empty type (s = 0, only for even n) is valid and can help. The original C++
    skips it (`if(sizeOfFlag > 0)`), so it is off by default to match.
    Returns a list of (sigma, flags)."""
    zeros = tuple(zeros)
    out = []
    for s in range(n % 2, n, 2):
        if s == 0 and not include_empty_type:
            continue
        m = (n + s) // 2
        for sigma in generate(s, c, zeros):
            flags = list(generate_flags(sigma, m, c, zeros))
            if flags:
                out.append((sigma, flags))
    return out
