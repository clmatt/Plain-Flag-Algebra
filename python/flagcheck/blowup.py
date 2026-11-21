"""Eliminating the clone ("green") variables of the blow-up encoding (reference version).

Let W be the balanced blow-up of a red/blue graph H on N vertices: N classes of measure
1/N, green inside a class. Sample n i.i.d. points of W. Their pattern of coincident
classes is a set partition with r blocks with probability (N)_r / N^n, and the classes
of the blocks are then r distinct uniform vertices of H. Hence, for a 3-coloured graph G
on n vertices with r green classes and green-free quotient Q_r(G),

    p_W(G) = n!/aut(G) * (N)_r/N^n * aut(Q_r)/r! * p_inj(Q_r),

where p_inj(Q) is the probability that |Q| distinct uniform vertices of H induce Q.
Writing z_Q = p_inj(Q) for the green-free graphs Q on n vertices,
p_inj(Q_r) = sum_Q d(Q_r, Q) z_Q (a uniform r-subset of a uniform n-subset is uniform), so

    y_G = sum_Q M[G][Q] z_Q,   M[G][Q] = lambda(G) d(Q_r(G), Q),

exactly, for every H with N >= n. z is a probability vector on the green-free graphs, so
an SDP over z (constraint vectors c replaced by c M) is a valid relaxation with fewer
variables. Every column of M sums to 1 (tests).
"""

from fractions import Fraction
from math import factorial
from typing import Dict, Tuple

from .algebra import lift
from .graphs import Flag, automorphism_count, canonical, generate, restrict
from .ramsey import GREEN


def falling(N: int, r: int) -> int:
    out = 1
    for i in range(r):
        out *= N - i
    return out


def quotient(G: Flag) -> Tuple[Tuple[Tuple[int, ...], ...], Flag]:
    """(green classes, canonical green-free quotient). Checks the blow-up structure: green
    is an equivalence relation and clones have the same colour to every other vertex."""
    n = G.k
    classes, seen = [], set()
    for v in range(n):
        if v in seen:
            continue
        cls = tuple(u for u in range(n) if u == v or G.color(u, v) == GREEN)
        seen.update(cls)
        classes.append(cls)
    for cls in classes:
        for u in cls:
            for w in cls:
                if u != w and G.color(u, w) != GREEN:
                    raise ValueError("green is not transitive")
    for a in classes:
        for b in classes:
            if a is not b:
                cols = {G.color(u, w) for u in a for w in b}
                if len(cols) != 1 or GREEN in cols:
                    raise ValueError("clones disagree")
    Q = canonical(restrict(G, [cls[0] for cls in classes], 0))
    return tuple(classes), Q


def coefficient(G: Flag, N: int) -> Tuple[Fraction, Flag]:
    """(lambda(G), Q_r(G)) with p_W(G) = lambda(G) p_inj(Q_r(G))."""
    classes, Q = quotient(G)
    n, r = G.k, len(classes)
    lam = Fraction(factorial(n), automorphism_count(G)) * Fraction(falling(N, r), N ** n) \
        * Fraction(automorphism_count(Q), factorial(r))
    return lam, Q


def is_green_free(G: Flag) -> bool:
    return GREEN not in G.cols


def blowup_map(n: int, zeros, N: int) -> Dict[Flag, Dict[Flag, Fraction]]:
    """M[G] = {Q: M[G][Q]} for every admissible G on n vertices; Q ranges over the
    admissible green-free graphs on n vertices."""
    zeros = tuple(zeros)
    M = {}
    for G in generate(n, 3, zeros):
        lam, Qr = coefficient(G, N)
        lifted = lift({Qr: Fraction(1)}, n, 3, zeros)
        M[G] = {Q: lam * d for Q, d in lifted.items() if is_green_free(Q)}
    return M


def reduce_vector(c: Dict[Flag, Fraction], M) -> Dict[Flag, Fraction]:
    """c . y = (c M) . z."""
    out: Dict[Flag, Fraction] = {}
    for G, v in c.items():
        for Q, m in M[G].items():
            out[Q] = out.get(Q, Fraction(0)) + v * m
    return {Q: v for Q, v in out.items() if v != 0}


def reduce_sdp(sdp, N: int, zeros):
    """The reduced FlagSDP over the green-free graphs (reference for flagsdp's `blowup`):
    Y_Q = sum_G M[G][Q] Y_G, and every objective / constraint vector c becomes c M."""
    from .sdp import FlagSDP

    M = blowup_map(sdp.n, zeros, N)
    index = {G: i for i, G in enumerate(sdp.graphs)}
    if set(index) != set(M):
        raise ValueError("graph lists differ")
    reduced = [G for G in sdp.graphs if is_green_free(G)]
    ri = {Q: i for i, Q in enumerate(reduced)}

    def vec(v):
        out = [Fraction(0)] * len(reduced)
        for G, g in index.items():
            if v[g]:
                for Q, m in M[G].items():
                    out[ri[Q]] += v[g] * m
        return out

    Y = [[{} for _ in sdp.blocks] for _ in reduced]
    for G, g in index.items():
        for Q, m in M[G].items():
            q = ri[Q]
            for b, entries in enumerate(sdp.Y[g]):
                for key, val in entries.items():
                    Y[q][b][key] = Y[q][b].get(key, Fraction(0)) + m * val
    Y = [[{k: v for k, v in d.items() if v != 0} for d in row] for row in Y]
    return FlagSDP(n=sdp.n, colors=sdp.colors, graphs=reduced, blocks=sdp.blocks, Y=Y, f=vec(sdp.f),
                   eqs=[(vec(v), r) for v, r in sdp.eqs], ineqs=[(vec(v), r) for v, r in sdp.ineqs],
                   block_kinds=sdp.block_kinds)


def column(Q: Flag, graphs, N: int) -> Dict[int, Fraction]:
    """{index g: M[graphs[g]][Q]} for one green-free Q, computed directly: lambda(G) times
    the number of r-subsets of Q inducing Q_r(G), over C(n, r). Independent of lift()."""
    from itertools import combinations
    from math import comb

    n = Q.k
    sub: Dict[Flag, int] = {}
    for r in range(1, n + 1):
        for X in combinations(range(n), r):
            F = canonical(restrict(Q, list(X), 0))
            sub[F] = sub.get(F, 0) + 1
    out = {}
    for g, G in enumerate(graphs):
        lam, Qr = coefficient(G, N)
        c = sub.get(Qr, 0)
        if c:
            out[g] = lam * Fraction(c, comb(n, Qr.k))
    return out
