"""The clone-elimination map y = M z must hold exactly on real blow-ups."""

import random
from fractions import Fraction
from itertools import combinations
from math import comb

import pytest

from flagcheck.blowup import blowup_map, coefficient, is_green_free, quotient, reduce_vector
from flagcheck.graphs import canonical, from_color_function, generate, restrict
from flagcheck.limits import StepGraphon
from flagcheck.ramsey import BLUE, GREEN, RED, blowup_zeros, paley, ramsey_zeros


def injective_densities(H, n):
    """z_Q = fraction of n-subsets of H inducing Q (by brute force)."""
    z = {}
    total = comb(H.k, n)
    for X in combinations(range(H.k), n):
        Q = canonical(restrict(H, list(X), 0))
        z[Q] = z.get(Q, Fraction(0)) + Fraction(1, total)
    return z


def random_graph(rng, N):
    edges = {(i, j): rng.choice((RED, BLUE)) for i, j in combinations(range(N), 2)}
    return from_color_function(N, lambda i, j: edges[(min(i, j), max(i, j))])


CASES = [(5, 3, 0), (6, 4, 1), (7, 4, 2), (6, 5, 3), (8, 5, 4)]


@pytest.mark.parametrize("N,n,seed", CASES)
def test_map_is_exact_on_random_blowups(N, n, seed):
    rng = random.Random(seed)
    zeros = blowup_zeros()
    M = blowup_map(n, zeros, N)
    for _ in range(3):
        H = random_graph(rng, N)
        y = StepGraphon.blowup(H, clone_color=GREEN).graph_distribution(n)
        z = injective_densities(H, n)
        assert set(y) <= set(M)                       # every sampled graph is admissible
        for G, row in M.items():
            assert y.get(G, Fraction(0)) == sum(m * z.get(Q, 0) for Q, m in row.items()), G


def test_map_is_exact_on_paley17_with_ramsey_zeros():
    H = paley(17)
    zeros = ramsey_zeros(4, 4)
    M = blowup_map(4, zeros, 17)
    y = StepGraphon.blowup(H, clone_color=GREEN).graph_distribution(4)
    z = injective_densities(H, 4)
    for G, row in M.items():
        assert y.get(G, Fraction(0)) == sum(m * z.get(Q, 0) for Q, m in row.items())


@pytest.mark.parametrize("n,N", [(4, 9), (5, 13), (5, 40)])
def test_columns_sum_to_one(n, N):
    M = blowup_map(n, ramsey_zeros(4, 5), N)
    cols = {}
    for row in M.values():
        for Q, m in row.items():
            cols[Q] = cols.get(Q, 0) + m
    green_free = [G for G in M if is_green_free(G)]
    assert set(cols) == set(green_free)
    assert all(v == 1 for v in cols.values())
    for G in green_free:                              # green-free graphs map to themselves
        assert M[G] == {G: Fraction(N ** 0) * coefficient(G, N)[0]} and coefficient(G, N)[0] == \
            Fraction(__import__("math").prod(range(N - n + 1, N + 1)), N ** n)


def test_quotient_rejects_non_blowups():
    bad = from_color_function(3, lambda i, j: GREEN if (i, j) == (0, 1) else (RED if (i, j) == (0, 2) else BLUE))
    with pytest.raises(ValueError):
        quotient(bad)


def test_reduce_vector_preserves_values():
    rng = random.Random(9)
    N, n = 7, 4
    M = blowup_map(n, blowup_zeros(), N)
    H = random_graph(rng, N)
    y = StepGraphon.blowup(H, clone_color=GREEN).graph_distribution(n)
    z = injective_densities(H, n)
    c = {G: Fraction(rng.randint(-5, 5), rng.randint(1, 5)) for G in M}
    cz = reduce_vector(c, M)
    assert sum(v * y.get(G, 0) for G, v in c.items()) == sum(v * z.get(Q, 0) for Q, v in cz.items())


# ---------------------------------------------------------------------------
# flagsdp's reduced output against the reference reduction
# ---------------------------------------------------------------------------

import dataclasses
import shutil

from flagcheck import cppbackend
from flagcheck.blowup import reduce_sdp
from flagcheck.problems import ramsey_feasibility
from flagcheck.ramsey_run import build_problem, run
from flagcheck.sdp import FlagSDP

ALL = ("green", "codeg", "tcodeg", "edges", "pedges", "nstats", "lpsd")
REDUCED_PROBLEMS = [
    dataclasses.replace(ramsey_feasibility(3, 4, 9, 5), blowup=9, name="R34-9-b"),
    build_problem(3, 4, 9, 5, "all", blowup=True),
    build_problem(3, 4, 8, 5, "none", extras=ALL, blowup=True),
    build_problem(4, 4, 17, 5, "local", extras=("codeg", "edges", "pedges", "lpsd"), blowup=True),
]


@pytest.mark.skipif(not cppbackend.available(), reason="flagsdp not built")
@pytest.mark.parametrize("problem", REDUCED_PROBLEMS, ids=lambda p: p.name)
def test_cpp_reduction_equals_reference(problem):
    """The unreduced C++ problem is checked exactly against the reference elsewhere; here
    its reduction by the Python map must equal flagsdp's reduced output exactly."""
    full = cppbackend.build_cpp_sdp(dataclasses.replace(problem, blowup=0))
    full = FlagSDP(**{**full.__dict__, "graphs": [canonical(G) for G in full.graphs]})
    ref = reduce_sdp(full, problem.blowup, problem.zeros)
    cpp = cppbackend.build_cpp_sdp(problem)
    assert [canonical(G) for G in cpp.graphs] == ref.graphs
    assert [s for s, _ in cpp.blocks] == [s for s, _ in ref.blocks]
    assert [f for _, f in cpp.blocks] == [f for _, f in ref.blocks]
    assert cpp.block_kinds == ref.block_kinds
    assert cpp.f == ref.f
    assert cpp.eqs == ref.eqs and cpp.ineqs == ref.ineqs
    assert cpp.Y == ref.Y
    assert len(cpp.graphs) < len(full.graphs)


@pytest.mark.skipif(not cppbackend.available() or shutil.which("csdp") is None, reason="needs flagsdp + CSDP")
@pytest.mark.parametrize("extras", [(), ALL])
def test_reduced_driver_end_to_end(extras):
    _, _, res9, _ = run(3, 4, 9, 5, "all", extras=extras, blowup=True)
    _, _, res8, _ = run(3, 4, 8, 5, "all", extras=extras, blowup=True)
    assert res9.bound > 0
    assert res8.bound <= 0


@pytest.mark.skipif(not cppbackend.available() or shutil.which("csdp") is None, reason="needs flagsdp + CSDP")
@pytest.mark.parametrize("kln", [(3, 5, 13, 6), (4, 4, 17, 6)])
def test_reduced_never_refutes_a_real_graph(kln):
    k, l, N, n = kln
    _, sol, res, _ = run(k, l, N, n, "all", extras=ALL, blowup=True)
    assert res.bound <= 0, sol.primal_objective
