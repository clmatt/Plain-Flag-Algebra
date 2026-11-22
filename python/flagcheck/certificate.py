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


def is_psd_exact(M) -> bool:
    """Exact PSD test of a symmetric rational matrix by LDL^T with diagonal pivoting.
    A zero pivot is only allowed when its whole remaining row is zero."""
    A = [[Fraction(x) for x in row] for row in M]
    n = len(A)
    active = list(range(n))
    while active:
        p = max(active, key=lambda i: A[i][i])
        d = A[p][p]
        if d < 0:
            return False
        if d == 0:
            if any(A[p][j] != 0 for j in active):
                return False
            active.remove(p)
            continue
        active.remove(p)
        for i in active:
            if A[i][p] == 0:
                continue
            f = A[i][p] / d
            for j in active:
                A[i][j] -= f * A[p][j]
    return True


def check_flagsdp_file(path: str, cert: Certificate) -> CheckResult:
    """Exactly the computation of check_certificate, streamed from a `flagsdp` output file
    so that problems too large to load as Fractions can be checked. Plain-block Y lines
    are integer counts (sorted by graph); localising-block lines are exact rationals."""
    with open(path) as fh:
        def nxt():
            return fh.readline().split()

        if nxt() != ["FLAGSDP", "1"]:
            raise ValueError("not a FLAGSDP v1 file")
        nxt()  # n
        nxt()  # colors
        m = int(nxt()[1])
        for _ in range(m):
            fh.readline()
        nb = int(nxt()[1])
        sizes, dens, kinds = [], [], []
        for b in range(nb):
            parts = nxt()
            if parts[0] != "block" or int(parts[1]) != b:
                raise ValueError("bad block header")
            sizes.append(int(parts[3]))
            dens.append(int(parts[4]))
            kinds.append(parts[5])
            for _ in range(sizes[-1]):
                fh.readline()
        f = [Fraction(0)] * m
        parts = nxt()
        for _ in range(int(parts[1])):
            g, v = nxt()
            f[int(g)] = Fraction(v)
        # The constraint rows are folded into the slacks as they are read instead of being
        # stored: at n = 8 they hold ~10^8 entries, far too many to keep as Fractions.
        #   r_g = f_g - ysum_g - sum_j y_j (a_j[g] - alpha_j) - sum_i u_i (beta_i - b_i[g])
        ys = [Fraction(y) for y in cert.y]
        us = [Fraction(u) for u in cert.u]
        if any(u < 0 for u in us):
            raise ValueError("inequality multipliers must be non-negative")
        corr = [Fraction(0)] * m      # -sum_j y_j a_j[g] + sum_i u_i b_i[g]
        const = Fraction(0)           # sum_j y_j alpha_j - sum_i u_i beta_i
        neq = nle = 0
        parts = nxt()
        while parts[0] == "constraint":
            if parts[1] == "eq":
                if neq >= len(ys):
                    raise ValueError("certificate does not match the problem shape")
                c, neq = -ys[neq], neq + 1
            else:
                if nle >= len(us):
                    raise ValueError("certificate does not match the problem shape")
                c, nle = us[nle], nle + 1
            const -= c * Fraction(parts[2])
            for _ in range(int(parts[3])):
                g, v = nxt()
                if c:
                    corr[int(g)] += c * Fraction(v)
            parts = nxt()
        if parts[0] != "Y":
            raise ValueError("expected Y section")
        if len(cert.blocks) != nb or len(ys) != neq or len(us) != nle:
            raise ValueError("certificate does not match the problem shape")
        grams = []
        for (rows, weights), d in zip(cert.blocks, sizes):
            if any(Fraction(w) < 0 for w in weights):
                raise ValueError("Gram weights must be non-negative")
            grams.append(_integer_gram(rows, weights, d))

        ysum = [Fraction(0)] * m
        cur_g, acc = -1, {}

        def flush():
            if cur_g >= 0:
                for b, num in acc.items():
                    ysum[cur_g] += Fraction(num, dens[b] * grams[b][1])

        seen, seen_loc = set(), set()
        for _ in range(int(parts[1])):
            g, b, i, j, v = fh.readline().split()
            g, b, i, j = int(g), int(b), int(i), int(j)
            if not (0 <= i <= j < sizes[b]) or not (0 <= g < m):
                raise ValueError(f"bad Y entry {g} {b} {i} {j}")
            mult = 1 if i == j else 2
            if kinds[b] == "plain":
                if g != cur_g:
                    if g < cur_g:
                        raise ValueError("plain Y entries are not sorted by graph")
                    flush()
                    cur_g, acc = g, {}
                    seen = set()
                if (b, i, j) in seen:
                    raise ValueError(f"duplicate Y entry {g} {b} {i} {j}")
                seen.add((b, i, j))
                acc[b] = acc.get(b, 0) + int(v) * grams[b][0][i, j] * mult
            else:
                if (g, b, i, j) in seen_loc:
                    raise ValueError(f"duplicate Y entry {g} {b} {i} {j}")
                seen_loc.add((g, b, i, j))
                ysum[g] += Fraction(v) * Fraction(int(grams[b][0][i, j]) * mult, dens[b] * grams[b][1])
        flush()
        if fh.readline().strip() != "end":
            raise ValueError("missing end marker")

    slacks = [f[g] - ysum[g] + corr[g] + const for g in range(m)]
    bound = min(slacks)
    return CheckResult(bound, slacks, [g for g, r in enumerate(slacks) if r == bound])
