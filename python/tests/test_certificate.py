"""The exact certificate checker is the final arbiter of every claim, so it is checked
against a deliberately naive implementation, and against hand-made certificates."""

import random
from fractions import Fraction

import pytest

from flagcheck.certificate import Certificate, check_certificate
from flagcheck.problems import P3, goodman, inducibility, ramsey_feasibility
from flagcheck.sdp import build_reference_sdp


def naive_slacks(sdp, cert):
    out = []
    mats = []
    for (sigma, flags), (rows, weights) in zip(sdp.blocks, cert.blocks):
        d = len(flags)
        A = [[sum((w * q[i] * q[j] for q, w in zip(rows, weights)), Fraction(0)) for j in range(d)] for i in range(d)]
        mats.append(A)
    for g in range(len(sdp.graphs)):
        r = sdp.f[g]
        for b, entries in enumerate(sdp.Y[g]):
            for (i, j), v in entries.items():
                r -= v * mats[b][i][j] * (1 if i == j else 2)
        for (a, alpha), y in zip(sdp.eqs, cert.y):
            r -= y * (a[g] - alpha)
        for (bv, beta), u in zip(sdp.ineqs, cert.u):
            r -= u * (beta - bv[g])
        out.append(r)
    return out


def random_cert(rng, sdp):
    def q():
        return Fraction(rng.randint(-50, 50), rng.choice([1, 3, 7, 1024, 10 ** 9 + 7]))

    blocks = []
    for _, flags in sdp.blocks:
        k = rng.randint(0, len(flags))
        rows = [[q() for _ in flags] for _ in range(k)]
        weights = [abs(q()) for _ in range(k)]
        blocks.append((rows, weights))
    y = [q() for _ in sdp.eqs]
    u = [abs(q()) for _ in sdp.ineqs]
    return Certificate(blocks, y, u)


@pytest.mark.parametrize("seed", range(6))
def test_fast_checker_matches_naive(seed):
    rng = random.Random(seed)
    p = [goodman(3), inducibility(P3, 4), ramsey_feasibility(3, 4, 9, 5)][seed % 3]
    sdp = build_reference_sdp(p.n, p.colors, p.zeros, p.objective, p.eqs, p.ineqs)
    cert = random_cert(rng, sdp)
    res = check_certificate(sdp, cert)
    assert res.slacks == naive_slacks(sdp, cert)
    assert res.bound == min(res.slacks)


def test_rejects_negative_weights_and_multipliers():
    p = ramsey_feasibility(3, 4, 9, 5)
    sdp = build_reference_sdp(p.n, p.colors, p.zeros, p.objective, p.eqs, p.ineqs)
    cert = random_cert(random.Random(1), sdp)
    bad = Certificate([(rows, [-w for w in ws] if ws else ws) for rows, ws in cert.blocks], cert.y, cert.u)
    if any(ws for _, ws in cert.blocks):
        with pytest.raises(ValueError):
            check_certificate(sdp, bad)


def test_hand_made_p3_certificate():
    """Thesis Ch. I: A = (3/4) [[1, -1], [-1, 1]] proves max P3 <= 3/4 on n = 3, i.e.
    -P3 >= -3/4 (block flags ordered as the reference generates them)."""
    p = inducibility(P3, 3)
    sdp = build_reference_sdp(p.n, p.colors, p.zeros, p.objective)
    (sigma, flags), = sdp.blocks
    assert len(flags) == 2
    cert = Certificate([([[Fraction(1), Fraction(-1)]], [Fraction(3, 4)])], [], [])
    assert check_certificate(sdp, cert).bound == Fraction(-3, 4)
