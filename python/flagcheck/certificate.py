"""Exact certificates for plain flag algebra bounds (see sdp.py for the convention).

A certificate is
  * for each SDP block, rational vectors q_1..q_r and rational weights D_1..D_r >= 0,
    so that A_sigma = sum_i D_i q_i q_i^T is PSD by construction;
  * rational multipliers y (equalities, free) and u >= 0 (inequalities).
The bound it proves is computed exactly as d = min_G r_G. Nothing numerical is trusted.
"""

from dataclasses import dataclass
from fractions import Fraction
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


def _gram(rows: List[List[Fraction]], weights: List[Fraction], size: int) -> List[List[Fraction]]:
    A = [[Fraction(0)] * size for _ in range(size)]
    for q, w in zip(rows, weights):
        if len(q) != size:
            raise ValueError("certificate vector has the wrong length")
        if w == 0:
            continue
        nz = [(i, x) for i, x in enumerate(q) if x]
        for i, xi in nz:
            wi = w * xi
            for j, xj in nz:
                A[i][j] += wi * xj
    return A


def check_certificate(sdp: FlagSDP, cert: Certificate) -> CheckResult:
    if len(cert.blocks) != len(sdp.blocks):
        raise ValueError("wrong number of blocks")
    if len(cert.y) != len(sdp.eqs) or len(cert.u) != len(sdp.ineqs):
        raise ValueError("wrong number of multipliers")
    if any(u < 0 for u in cert.u):
        raise ValueError("inequality multipliers must be non-negative")
    mats = []
    for (sigma, flags), (rows, weights) in zip(sdp.blocks, cert.blocks):
        if len(rows) != len(weights):
            raise ValueError("rows/weights mismatch")
        if any(w < 0 for w in weights):
            raise ValueError("Gram weights must be non-negative")
        mats.append(_gram(rows, weights, len(flags)))

    slacks = []
    for g in range(len(sdp.graphs)):
        r = sdp.f[g]
        for b, entries in enumerate(sdp.Y[g]):
            A = mats[b]
            for (i, j), v in entries.items():
                r -= v * A[i][j] * (1 if i == j else 2)
        for (a, alpha), y in zip(sdp.eqs, cert.y):
            r -= y * (a[g] - alpha)
        for (bvec, beta), u in zip(sdp.ineqs, cert.u):
            r -= u * (beta - bvec[g])
        slacks.append(r)
    bound = min(slacks)
    return CheckResult(bound, slacks, [g for g, r in enumerate(slacks) if r == bound])


def _rat(x: float, max_den: int) -> Fraction:
    return Fraction(x).limit_denominator(max_den)


def round_solution(sol: NumericSolution, max_den: int = 10 ** 9, eig_floor: float = 0.0) -> Certificate:
    """Turn a numerical CSDP solution into a rational certificate. The result is always a
    *valid* certificate; rounding only affects how strong the proven bound is."""
    import numpy as np

    blocks = []
    for A in sol.blocks:
        A = (A + A.T) / 2
        vals, vecs = np.linalg.eigh(A)
        rows, weights = [], []
        for lam, v in zip(vals, vecs.T):
            if lam <= eig_floor:
                continue
            rows.append([_rat(float(x), max_den) for x in v])
            weights.append(_rat(float(lam), max_den))
        blocks.append((rows, weights))
    L = sol.layout
    y = [_rat(sol.lp[L.y_plus(j)] - sol.lp[L.y_minus(j)], max_den) for j in range(L.num_eqs)]
    u = [max(Fraction(0), _rat(sol.lp[L.u(i)], max_den)) for i in range(L.num_ineqs)]
    return Certificate(blocks, y, u)
