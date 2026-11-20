"""Exact certificates for plain flag algebra bounds (see sdp.py for the convention).

A certificate is
  * for each SDP block, rational vectors q_1..q_r and rational weights D_1..D_r >= 0,
    so that A_sigma = sum_i D_i q_i q_i^T is PSD by construction;
  * rational multipliers y (equalities, free) and u >= 0 (inequalities).
The bound it proves is computed exactly as d = min_G r_G. Nothing numerical is trusted.

All arithmetic is exact. Each block is scaled to integers by the LCM of its
denominators, so the Gram matrices are big-integer matrices; Fractions are only formed
once per graph.
"""

from dataclasses import dataclass
from fractions import Fraction
from math import lcm
from typing import List, Tuple

from .sdp import FlagSDP, NumericSolution


@dataclass
class Certificate:
    blocks: List[Tuple[List[List[Fraction]], List[Fraction]]]  # per block: (rows q_i, weights D_i)
    y: List[Fraction]
    u: List[Fraction]


@dataclass
class CheckResult:
    bound: Fraction             # the proven d: f . s >= bound for every feasible s
    slacks: List[Fraction]      # r_G for every graph
    tight: List[int]            # graphs attaining the minimum


def _integer_gram(rows, weights, size):
    """Return (A_int, scale) with A_sigma = A_int / scale exactly (A_int: object ndarray)."""
    import numpy as np

    L = 1
    for q in rows:
        for x in q:
            L = lcm(L, Fraction(x).denominator)
    M = 1
    for w in weights:
        M = lcm(M, Fraction(w).denominator)
    A = np.zeros((size, size), dtype=object)
    A[:, :] = 0
    for q, w in zip(rows, weights):
        if len(q) != size:
            raise ValueError("certificate vector has the wrong length")
        w = Fraction(w)
        if w == 0:
            continue
        wi = w.numerator * (M // w.denominator)
        Q = np.array([int(Fraction(x) * L) for x in q], dtype=object)
        A += wi * np.outer(Q, Q)
    return A, M * L * L


def check_certificate(sdp: FlagSDP, cert: Certificate) -> CheckResult:
    if len(cert.blocks) != len(sdp.blocks):
        raise ValueError("wrong number of blocks")
    if len(cert.y) != len(sdp.eqs) or len(cert.u) != len(sdp.ineqs):
        raise ValueError("wrong number of multipliers")
    if any(Fraction(u) < 0 for u in cert.u):
        raise ValueError("inequality multipliers must be non-negative")
    grams = []
    for (sigma, flags), (rows, weights) in zip(sdp.blocks, cert.blocks):
        if len(rows) != len(weights):
            raise ValueError("rows/weights mismatch")
        if any(Fraction(w) < 0 for w in weights):
            raise ValueError("Gram weights must be non-negative")
        grams.append(_integer_gram(rows, weights, len(flags)))

    slacks = []
    for g in range(len(sdp.graphs)):
        r = Fraction(sdp.f[g])
        for b, entries in enumerate(sdp.Y[g]):
            if not entries:
                continue
            A, scale = grams[b]
            # Y entries share a denominator per block in practice; group by denominator.
            by_den = {}
            for (i, j), v in entries.items():
                v = Fraction(v)
                term = v.numerator * A[i, j] * (1 if i == j else 2)
                by_den[v.denominator] = by_den.get(v.denominator, 0) + term
            for den, num in by_den.items():
                r -= Fraction(num, den * scale)
        for (a, alpha), y in zip(sdp.eqs, cert.y):
            if y:
                r -= Fraction(y) * (a[g] - alpha)
        for (bvec, beta), u in zip(sdp.ineqs, cert.u):
            if u:
                r -= Fraction(u) * (beta - bvec[g])
        slacks.append(r)
    bound = min(slacks)
    return CheckResult(bound, slacks, [g for g, r in enumerate(slacks) if r == bound])


def round_solution(sol: NumericSolution, bits: int = 40, eig_floor: float = 0.0,
                   max_den: int = None) -> Certificate:
    """Turn a numerical CSDP solution into a rational certificate. The result is always a
    *valid* certificate; rounding only affects how strong the proven bound is.

    Entries are rounded to the dyadic grid 2^-bits, so every block has one common
    denominator and the exact check stays fast. (max_den is accepted for backwards
    compatibility and translated to the nearest number of bits.)"""
    import numpy as np

    if max_den is not None:
        bits = max(1, int(max_den).bit_length())
    D = 1 << bits

    def rat(x: float) -> Fraction:
        return Fraction(round(x * D), D)

    blocks = []
    for A in sol.blocks:
        A = (A + A.T) / 2
        vals, vecs = np.linalg.eigh(A)
        rows, weights = [], []
        for lam, v in zip(vals, vecs.T):
            if lam <= eig_floor:
                continue
            w = rat(float(lam))
            if w <= 0:
                continue
            rows.append([rat(float(x)) for x in v])
            weights.append(w)
        blocks.append((rows, weights))
    L = sol.layout
    y = [rat(sol.lp[L.y_plus(j)] - sol.lp[L.y_minus(j)]) for j in range(L.num_eqs)]
    u = [max(Fraction(0), rat(sol.lp[L.u(i)])) for i in range(L.num_ineqs)]
    return Certificate(blocks, y, u)
