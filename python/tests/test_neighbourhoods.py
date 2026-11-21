"""Soundness and correctness of the neighbourhood constraints (green / codeg / edges), the
edge tables they rely on, and the direct evaluator used by verify."""

import random
from fractions import Fraction
from math import comb

import pytest

from flagcheck import neighbourhoods as nb
from flagcheck.algebra import lift
from flagcheck.direct import GraphEvaluator
from flagcheck.edge_tables import EDGE_BOUNDS, brute_bounds, edge_bounds, is_ab_graph
from flagcheck.graphs import Flag, canonical, generate, generate_flags
from flagcheck.limits import StepGraphon
from flagcheck.ramsey import (
    BLUE,
    KNOWN_RAMSEY_GRAPHS,
    RED,
    circulant,
    is_ramsey_graph,
    paley,
    ramsey_zeros,
    red_degrees,
)
from flagcheck.ramsey_run import build_problem, degree_window, ramsey_upper
from flagcheck.spec import Add, Avg, Lit, Mul, Scale, evaluate

# ---------------------------------------------------------------------------
# Edge tables
# ---------------------------------------------------------------------------

# Goedgebeur & Radziszowski, EJC 20(1) (2013) #P30, Table 3 ("Exact values of e(3,k,n)"),
# an independent source for the minima computed from McKay's lists.
GR_TABLE3 = {
    4: {4: 1, 5: 2, 6: 3, 7: 6, 8: 10},
    5: {5: 1, 6: 2, 7: 3, 8: 4, 9: 7, 10: 10, 11: 15, 12: 20, 13: 26},
    6: {6: 1, 7: 2, 8: 3, 9: 4, 10: 5, 11: 8, 12: 11, 13: 15, 14: 20, 15: 25, 16: 32, 17: 40},
}


def test_minima_match_goedgebeur_radziszowski():
    for k, row in GR_TABLE3.items():
        for n, e in row.items():
            assert EDGE_BOUNDS[(3, k)][n][0] == e, (k, n)
        # below the table's range, e(3,k,n) = 0 (an independent set of size n < k)
        for n in range(1, k):
            assert EDGE_BOUNDS[(3, k)][n][0] == 0


def test_table_orders_end_at_the_ramsey_number():
    for (a, b), row in EDGE_BOUNDS.items():
        assert max(row) == ramsey_upper(a, b) - 1, (a, b)
        assert sorted(row) == list(range(1, max(row) + 1)), (a, b)
        for n, (lo, hi) in row.items():
            assert 0 <= lo <= hi <= comb(n, 2)


@pytest.mark.parametrize("ab", [(3, 3), (3, 4), (3, 5), (3, 6), (4, 4), (4, 5)])
def test_small_orders_by_brute_force(ab):
    a, b = ab
    for n in range(1, 7):
        bf = brute_bounds(a, b, n)
        if bf is None:
            with pytest.raises(KeyError):
                edge_bounds(a, b, n)
        else:
            assert edge_bounds(a, b, n) == bf, (ab, n)
            assert edge_bounds(b, a, n) == brute_bounds(b, a, n), (ab, n)


def test_trivial_families():
    assert edge_bounds(2, 5, 4) == (0, 0)
    assert edge_bounds(5, 2, 4) == (6, 6)
    with pytest.raises(KeyError):
        edge_bounds(2, 5, 5)
    assert edge_bounds(1, 7, 0) == (0, 0)


def _neighbourhood_edges(G, v, color):
    S = [u for u in range(G.k) if u != v and G.color(u, v) == color]
    return len(S), sum(1 for i in S for j in S if i < j and G.color(i, j) == RED)


@pytest.mark.parametrize("kl", sorted(KNOWN_RAMSEY_GRAPHS))
def test_real_neighbourhoods_lie_in_the_tables(kl):
    k, l = kl
    G = KNOWN_RAMSEY_GRAPHS[kl]()
    for v in range(G.k):
        for color, (a, b) in ((RED, (k - 1, l)), (BLUE, (k, l - 1))):
            size, e = _neighbourhood_edges(G, v, color)
            lo, hi = edge_bounds(a, b, size)
            assert lo <= e <= hi


def test_paley17_neighbourhood_attains_the_maximum():
    # Paley(17): red neighbourhoods are 3-regular on 8 vertices = the max e for (3,4).
    G = paley(17)
    assert _neighbourhood_edges(G, 0, RED) == (8, 12) and edge_bounds(3, 4, 8) == (10, 12)


def test_g6_parser():
    from flagcheck.edge_tables import g6_edge_count, parse_g6

    # C5 in graph6 is "Dhc"; K4 is "C~".
    assert parse_g6(b"C~") == (4, [(0, 1), (0, 2), (1, 2), (0, 3), (1, 3), (2, 3)])
    n, edges = parse_g6(b"Dhc")
    assert n == 5 and len(edges) == 5 and g6_edge_count(b"Dhc") == (5, 5)
    assert is_ab_graph(5, edges, 3, 3)


# ---------------------------------------------------------------------------
# Direct evaluator vs the reference algebra
# ---------------------------------------------------------------------------

def test_direct_evaluator_matches_reference():
    from flagcheck.ramsey import one_flag, red_flag

    zeros = tuple(ramsey_zeros(3, 4))
    rng = random.Random(1)
    graphs = list(generate(5, 3, zeros))
    g2 = nb.codegree_multipliers(3, 4, 8, ramsey_upper)[0][0][1]
    exprs = [
        Avg(Mul((Add((Lit(red_flag()), Scale(Fraction(-1, 3), one_flag()))), Lit(red_flag())))),
        Avg(Mul((g2, Lit(rng.choice(generate_flags(nb.sigma_of(g2), 3, 3, zeros)))))),
        Avg(nb.green_multipliers(8)[1][1]),
        Add((Lit(Flag(3, 0, (RED, RED, BLUE))), Scale(2, Avg(Lit(red_flag()))))),
        Avg(Mul((g2, Lit(generate_flags(nb.sigma_of(g2), 3, 3, zeros)[0]),
                 Lit(generate_flags(nb.sigma_of(g2), 3, 3, zeros)[-1])))),
    ]
    for e in exprs:
        ref = lift(evaluate(e, 3, zeros), 5, 3, zeros)
        for G in rng.sample(graphs, 15):
            assert GraphEvaluator(G).value(e) == ref.get(canonical(G), Fraction(0)), e


# ---------------------------------------------------------------------------
# Multipliers are valid pointwise on blow-ups of real Ramsey graphs
# ---------------------------------------------------------------------------

def _pointwise_values(W, g, s):
    gc = evaluate(g, 3, ())
    return [W.evaluate(gc, parts) for parts, _ in W.root_placements(s)]


CASES = [((3, 3), circulant(5, [1])), ((3, 4), circulant(8, [1, 4])), ((3, 5), circulant(13, [1, 5])),
         ((4, 4), paley(17)), ((3, 4), circulant(5, [1]))]   # last: C5 is a (3,4)-graph, N = 5 < R - 1


@pytest.mark.parametrize("kl,G", CASES, ids=lambda x: str(x) if isinstance(x, tuple) else str(x.k))
def test_multipliers_hold_pointwise(kl, G):
    k, l = kl
    N = G.k
    assert is_ramsey_graph(G, k, l)
    W = StepGraphon.blowup(G)
    lo, hi = degree_window(k, l, N)
    assert lo <= min(red_degrees(G)) and max(red_degrees(G)) <= hi
    for _, g in nb.green_multipliers(N):
        assert all(v == 0 for v in _pointwise_values(W, g, nb.sigma_of(g).k))
    codeg, _ = nb.codegree_multipliers(k, l, N, ramsey_upper)
    edges, _ = nb.edge_multipliers(k, l, N, range(lo, hi + 1), 4)
    for name, g in codeg + edges:
        vals = _pointwise_values(W, g, nb.sigma_of(g).k)
        assert min(vals) >= 0, name


def test_multipliers_are_tight_on_paley17():
    """Paley(17) attains the codegree bound for red pairs (3 = R(2,4) - 1) and the maximum
    edge count in red neighbourhoods; raising either bound by one would be violated."""
    G = paley(17)
    W = StepGraphon.blowup(G)
    lo, hi = degree_window(4, 4, 17)
    codeg, _ = nb.codegree_multipliers(4, 4, 17, ramsey_upper)
    edges, _ = nb.edge_multipliers(4, 4, 17, range(lo, hi + 1), 4)
    tight = {name for name, g in codeg + edges if min(_pointwise_values(W, g, nb.sigma_of(g).k)) == 0}
    assert "codeg[r:rr<=3]" in tight
    assert "red-nbhd-max-exact" in tight
    # a bound that is too strong by 1/N fails
    from flagcheck.algebra import identity_flag

    bad = Add((Scale(Fraction(2, 17), Lit(identity_flag(nb.RED_PAIR))),
               Scale(-1, Lit(nb.extension(nb.RED_PAIR, (RED, RED))))))
    assert min(_pointwise_values(W, bad, 2)) < 0


def test_lagrange_interpolant_is_exact_on_the_window():
    """C5 as a (3,4)-graph with N = 5: window [0, 3] gives a cubic interpolant; it must be
    0 where the bound is attained and detect a wrong table value."""
    G = circulant(5, [1])
    W = StepGraphon.blowup(G)
    lo, hi = degree_window(3, 4, 5)
    assert (lo, hi) == (0, 3)
    edges, _ = nb.edge_multipliers(3, 4, 5, range(lo, hi + 1), 5)
    names = [nm for nm, _ in edges]
    assert "blue-nbhd-min-exact" in names and "blue-nbhd-max-exact" in names
    vals = {nm: _pointwise_values(W, g, 1) for nm, g in edges}
    assert all(min(v) >= 0 for v in vals.values())
    # blue neighbourhood of a C5 vertex: 2 vertices joined in red; table max e(3,3,2) = 1
    assert min(vals["blue-nbhd-max-exact"]) == 0


def test_codegree_bounds():
    # R(4,6): pairs joined in red; S red to both is a (2,6)-graph, etc.
    cb = nb.codegree_bound
    assert cb(4, 6, RED, (RED, RED), ramsey_upper) == 5
    assert cb(4, 6, RED, (RED, BLUE), ramsey_upper) == 13
    assert cb(4, 6, RED, (BLUE, BLUE), ramsey_upper) == 24    # blue pair of roots not a blue K2
    assert cb(4, 6, BLUE, (BLUE, BLUE), ramsey_upper) == 17
    assert cb(4, 6, BLUE, (RED, RED), ramsey_upper) == 17
    assert cb(4, 6, BLUE, (RED, BLUE), ramsey_upper) == 13


# ---------------------------------------------------------------------------
# Whole problems
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("kl,n", [((3, 4), 5), ((3, 5), 4), ((4, 4), 4)])
def test_every_extra_constraint_holds_on_a_critical_graph(kl, n):
    k, l = kl
    G = KNOWN_RAMSEY_GRAPHS[kl]()
    W = StepGraphon.blowup(G)
    p = build_problem(k, l, G.k, n, "none", extras=("green", "codeg", "edges", "lpsd"))
    assert p.ineqs and p.assumptions
    for expr, rhs in p.eqs:
        assert W.evaluate(evaluate(expr, 3, p.zeros)) == rhs
    for expr, rhs in p.ineqs:
        assert W.evaluate(evaluate(expr, 3, p.zeros)) <= rhs
    for m, g in p.localizing:
        s = nb.sigma_of(g).k
        assert min(_pointwise_values(W, g, s)) >= 0
