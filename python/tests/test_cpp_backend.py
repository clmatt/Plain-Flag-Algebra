"""The refactored C++ core must reproduce the independent Python reference EXACTLY:
the same graphs, complete flag lists, identical Y entries, identical constraint vectors.
Build first:  cmake -S cpp -B ~/build/flagalg -G Ninja && cmake --build ~/build/flagalg
"""

import shutil
from fractions import Fraction

import pytest

from flagcheck import cppbackend
from flagcheck.algebra import Y_matrix, combo
from flagcheck.graphs import Flag, canonical, generate, generate_flags
from flagcheck.problems import (
    C5,
    P4,
    Problem,
    goodman,
    inducibility,
    mantel,
    ramsey_feasibility,
    ramsey_objective,
)
from flagcheck.ramsey import ROOT, monochromatic_clique, one_flag, ramsey_zeros, red_flag, GREEN
from flagcheck.ramsey_run import build_problem
from flagcheck.sdp import build_from_problem, solve_csdp
from flagcheck.spec import Add, Avg, Lit, Mul, Scale

pytestmark = pytest.mark.skipif(not cppbackend.available(), reason="flagsdp not built")


def degree_problem(N, n, window, localized_flag=None):
    """Ramsey (3,4) feasibility with a degree polynomial written as an expression tree,
    so C++'s multiply/average is exercised (and cross-checked by the reference)."""
    zeros = tuple(ramsey_zeros(3, 4))
    factors = [Add((Lit(red_flag()), Scale(Fraction(-d, N), one_flag()))) for d in range(window[0], window[1] + 1)]
    if localized_flag is not None:
        factors.append(Lit(localized_flag))
    eqs = [(combo((monochromatic_clique(2, GREEN), 1)), Fraction(1, N)),
           (Avg(Mul(tuple(factors))), Fraction(0))]
    return Problem(f"deg{N}", n, 3, zeros, {}, eqs=eqs, normalize=True)


def localizing_problem(N, n):
    """(3,4) feasibility with localising blocks for R - 2/N >= 0 (m = 2) and for the
    two-sided window (R - 2/N)(4/N - R) >= 0 (m = 1, i.e. a 1x1 block)."""
    p = ramsey_feasibility(3, 4, N, n)
    lower = Add((Lit(red_flag()), Scale(Fraction(-2, N), one_flag())))
    upper = Add((Scale(Fraction(4, N), one_flag()), Scale(-1, Lit(red_flag()))))
    p.localizing = [(2, lower), (1, Mul((lower, upper)))]
    p.name = f"loc{N}"
    return p


PROBLEMS = [
    goodman(3),
    goodman(4),
    mantel(4),
    inducibility(P4, 5),
    inducibility(C5, 5),
    Problem("empty-type-4", 4, 2, (), combo((P4, -1)), include_empty_type=True),
    Problem("empty-type-6", 6, 2, (), combo((P4, -1)), include_empty_type=True),
    Problem("three-colours-4", 4, 3, (), combo((monochromatic_clique(3, 1), 1))),
    ramsey_objective(3, 4, 5),
    ramsey_feasibility(4, 4, 18, 5),
    ramsey_feasibility(3, 5, 14, 6),
    degree_problem(8, 5, (2, 4)),
    degree_problem(8, 5, (3, 3), localized_flag=Flag(3, 1, (1, 2, 0))),
    localizing_problem(8, 5),
    # neighbourhood constraints: 2-rooted multipliers, linear and localising
    build_problem(3, 4, 8, 5, "none", extras=("green", "codeg", "edges", "lpsd")),
    build_problem(3, 4, 5, 5, "none", extras=("edges", "lpsd")),        # cubic interpolant
]


@pytest.mark.parametrize("problem", PROBLEMS, ids=lambda p: p.name)
def test_cpp_matches_reference_exactly(problem):
    ref = build_from_problem(problem)
    cpp = cppbackend.build_cpp_sdp(problem)
    assert cpp.block_kinds == (ref.block_kinds or ["plain"] * len(ref.blocks))
    plain = [b for b, k in enumerate(cpp.block_kinds) if k == "plain"]
    loc_cpp = [b for b, k in enumerate(cpp.block_kinds) if k == "localizing"]
    loc_ref = [b for b, k in enumerate(ref.block_kinds) if k == "localizing"]

    # 1. Same graphs; build the index map cpp -> ref through Python's canonical form.
    ref_index = {G: i for i, G in enumerate(ref.graphs)}
    to_ref = [ref_index[canonical(G)] for G in cpp.graphs]
    assert sorted(to_ref) == list(range(len(ref.graphs)))

    # 2. Objective and constraint vectors identical.
    for g, r in enumerate(to_ref):
        assert cpp.f[g] == ref.f[r]
    assert len(cpp.eqs) == len(ref.eqs) and len(cpp.ineqs) == len(ref.ineqs)
    for (cv, crhs), (rv, rrhs) in zip(cpp.eqs + cpp.ineqs, ref.eqs + ref.ineqs):
        assert crhs == rrhs
        for g, r in enumerate(to_ref):
            assert cv[g] == rv[r], (problem.name, g)

    # 3. Same type classes; each block is the complete, duplicate-free flag list.
    cpp_types = sorted(canonical(Flag(cpp.blocks[b][0].k, 0, cpp.blocks[b][0].cols)) for b in plain)
    ref_types = sorted(canonical(Flag(s.k, 0, s.cols)) for (s, _), k in zip(ref.blocks, ref.block_kinds) if k == "plain")
    assert cpp_types == ref_types
    for sigma, flags in cpp.blocks:
        canon = [canonical(F) for F in flags]
        assert len(set(canon)) == len(canon)
        assert all(F.type().cols == sigma.cols and F.s == sigma.k for F in flags)
        m = flags[0].k
        assert set(canon) == set(generate_flags(sigma, m, problem.colors, problem.zeros))

    # 4. Every Y entry, recomputed independently for C++'s labelled types and flag lists.
    for g, G in enumerate(cpp.graphs):
        for b in plain:
            sigma, flags = cpp.blocks[b]
            full = Y_matrix(G, sigma, flags)
            expected = {(i, j): full[i][j] for i in range(len(flags)) for j in range(i, len(flags)) if full[i][j]}
            assert cpp.Y[g][b] == expected, (problem.name, g, b)

    # 5. Localising blocks: same labelled type, and equal entries after mapping flags.
    assert len(loc_cpp) == len(loc_ref)
    for cb, rb in zip(loc_cpp, loc_ref):
        assert cpp.blocks[cb][0].cols == ref.blocks[rb][0].cols
        rindex = {F: i for i, F in enumerate(ref.blocks[rb][1])}
        perm = [rindex[canonical(F)] for F in cpp.blocks[cb][1]]
        assert sorted(perm) == list(range(len(perm)))
        for g, r in enumerate(to_ref):
            mapped = {}
            for (i, j), v in cpp.Y[g][cb].items():
                a, c = perm[i], perm[j]
                mapped[(min(a, c), max(a, c))] = v
            assert mapped == ref.Y[r][rb], (problem.name, g, cb)


@pytest.mark.skipif(shutil.which("csdp") is None, reason="needs CSDP")
@pytest.mark.parametrize("problem", [goodman(3), inducibility(P4, 5), ramsey_feasibility(4, 4, 18, 6)],
                         ids=lambda p: p.name)
def test_cpp_and_reference_solve_to_same_value(problem):
    ref = build_from_problem(problem)
    cpp = cppbackend.build_cpp_sdp(problem)
    a = solve_csdp(ref, normalize=problem.normalize).primal_objective
    b = solve_csdp(cpp, normalize=problem.normalize).primal_objective
    assert abs(a - b) < 1e-7 + 1e-6 * abs(a)
