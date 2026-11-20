"""Soundness of the Ramsey encoding: blow-ups of real Ramsey graphs must satisfy every
constraint we would ever add, exactly."""

from fractions import Fraction

import pytest

from flagcheck.algebra import lift
from flagcheck.graphs import generate
from flagcheck.limits import StepGraphon
from flagcheck.ramsey import (
    BLUE,
    GREEN,
    KNOWN_RAMSEY_GRAPHS,
    RED,
    blowup_zeros,
    circulant,
    degree_polynomial,
    green_clique_constraint,
    is_ramsey_graph,
    monochromatic_clique,
    ramsey_zeros,
    red_degrees,
)

RAMSEY_VALUES = {(3, 3): 6, (3, 4): 9, (3, 5): 14, (4, 4): 18}


@pytest.mark.parametrize("kl", sorted(KNOWN_RAMSEY_GRAPHS))
def test_known_graphs_are_critical_ramsey_graphs(kl):
    k, l = kl
    G = KNOWN_RAMSEY_GRAPHS[kl]()
    assert G.k == RAMSEY_VALUES[kl] - 1
    assert is_ramsey_graph(G, k, l)


def test_clique_checker_finds_cliques():
    # C6 has no triangle, but its complement does (two of them).
    C6 = circulant(6, [1])
    assert not is_ramsey_graph(C6, 3, 3)
    assert is_ramsey_graph(circulant(5, [1]), 3, 3)
    # Paley(17) has red triangles, so it is not a (3, 4)-graph; R(3,4) = 9 so no
    # 9-vertex circulant is a (3,4)-graph either.
    assert not is_ramsey_graph(KNOWN_RAMSEY_GRAPHS[(4, 4)](), 3, 4)
    assert not any(is_ramsey_graph(circulant(9, S), 3, 4) for S in ([1], [1, 4], [2, 3], [1, 2]))


def test_blowup_zeros_shape():
    # Exactly the three 3-vertex patterns with >= 1 green edge that break the class structure.
    zs = generate(3, 3, blowup_zeros())
    for G in zs:
        greens = sum(1 for c in G.cols if c == GREEN)
        if greens == 2:
            pytest.fail(f"two green edges must force the third: {G}")
        if greens == 1:
            others = [c for c in G.cols if c != GREEN]
            assert others[0] == others[1], G
    # 10 colourings of K3 in total; 3 are forbidden.
    assert len(zs) == 10 - 3


@pytest.mark.parametrize("kl,n", [((3, 3), 5), ((3, 4), 4), ((3, 5), 4), ((4, 4), 3)])
def test_blowup_satisfies_zeros_and_green_densities(kl, n):
    k, l = kl
    G = KNOWN_RAMSEY_GRAPHS[kl]()
    N = G.k
    W = StepGraphon.blowup(G)
    admissible = set(generate(n, 3, ramsey_zeros(k, l)))
    dist = W.graph_distribution(n)
    assert set(dist) <= admissible
    for j in range(2, n + 1):
        expr, rhs = green_clique_constraint(j, N)
        assert W.evaluate(expr) == rhs


@pytest.mark.parametrize("kl", [(3, 3), (3, 4), (3, 5), (4, 4)])
def test_degree_polynomial_vanishes_on_real_graphs(kl):
    k, l = kl
    G = KNOWN_RAMSEY_GRAPHS[kl]()
    N = G.k
    zeros = ramsey_zeros(k, l)
    W = StepGraphon.blowup(G)
    degrees = sorted(set(red_degrees(G)))
    # The theoretical window from smaller Ramsey numbers always contains the real degrees;
    # use the real degree set plus one neighbour to exercise a multi-factor product.
    allowed = degrees + [degrees[-1] + 1]
    expr = degree_polynomial(allowed, N, zeros)
    assert W.evaluate(expr) == 0
    # The same polynomial lifted by the chain rule must still evaluate to 0.
    assert W.evaluate(lift(expr, 4, 3, zeros)) == 0


def test_degree_polynomial_detects_wrong_degree_sets():
    """Regression for the mainRamsey.cpp bug (deg24 built from deg23Vec): a polynomial
    whose roots miss a real degree must NOT vanish on a real graph."""
    G = KNOWN_RAMSEY_GRAPHS[(3, 4)]()        # Wagner graph: 3-regular on 8 vertices
    zeros = ramsey_zeros(3, 4)
    W = StepGraphon.blowup(G)
    assert set(red_degrees(G)) == {3}
    good = degree_polynomial([2, 3, 4], 8, zeros)
    buggy = degree_polynomial([2, 4, 4], 8, zeros)   # "deg3" accidentally replaced
    assert W.evaluate(good) == 0
    assert W.evaluate(buggy) != 0
    # Value of the buggy one: (3/8-2/8)(3/8-4/8)^2 = (1/8)(1/64)
    assert W.evaluate(buggy) == Fraction(1, 8) * Fraction(1, 64)


def test_unbalanced_blowups_have_green_density_at_least_one_over_N():
    """Cauchy-Schwarz: sum w_i^2 >= 1/N with equality iff balanced (thesis Ch. II.2)."""
    G = KNOWN_RAMSEY_GRAPHS[(3, 3)]()
    expr, _ = green_clique_constraint(2, 5)
    balanced = StepGraphon.blowup(G).evaluate(expr)
    skewed = StepGraphon.blowup(G, [Fraction(1, 10)] * 4 + [Fraction(6, 10)]).evaluate(expr)
    assert balanced == Fraction(1, 5)
    assert skewed > Fraction(1, 5)
