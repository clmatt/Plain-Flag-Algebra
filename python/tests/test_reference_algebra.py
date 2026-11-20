"""Self-tests of the Python reference: known counts, and exact flag-algebra identities on
step graphons. A step graphon is a genuine limit object, so every identity must hold
exactly, with no tolerance."""

import random
from fractions import Fraction

import pytest

from flagcheck.algebra import (
    Y_matrix,
    average,
    average_coefficient,
    combo,
    flag_density,
    lift,
    multiply,
    sdp_blocks,
)
from flagcheck.graphs import (
    Flag,
    automorphism_count,
    canonical,
    from_edges,
    generate,
    generate_flags,
    restrict,
)
from flagcheck.limits import StepGraphon

# 2 colours: 0 = non-edge, 1 = edge.
EDGE = from_edges(2, [(0, 1, 1)])
K3 = from_edges(3, [(0, 1, 1), (0, 2, 1), (1, 2, 1)])
EMPTY3 = from_edges(3, [])


def random_step_graphon(rng, parts, colors, diag_colors=None):
    weights = [rng.randint(1, 5) for _ in range(parts)]
    total = sum(weights)
    weights = [Fraction(w, total) for w in weights]
    C = [[0] * parts for _ in range(parts)]
    for a in range(parts):
        for b in range(a, parts):
            pool = diag_colors if (a == b and diag_colors is not None) else range(colors)
            C[a][b] = C[b][a] = rng.choice(list(pool))
    return StepGraphon(weights, C)


# ---------------------------------------------------------------------------
# Counts
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("n,count", [(1, 1), (2, 2), (3, 4), (4, 11), (5, 34), (6, 156)])
def test_number_of_graphs(n, count):
    # OEIS A000088
    assert len(generate(n, 2)) == count


@pytest.mark.parametrize("n,count", [(1, 1), (2, 2), (3, 3), (4, 7), (5, 14), (6, 38)])
def test_number_of_triangle_free_graphs(n, count):
    # OEIS A006785
    assert len(generate(n, 2, [K3])) == count


@pytest.mark.parametrize("n,count", [(2, 3), (3, 10), (4, 66)])
def test_number_of_3_edge_colourings_of_complete_graphs(n, count):
    # OEIS A063841 (complete graphs with edges coloured in 3 colours, up to isomorphism)
    assert len(generate(n, 3)) == count


def test_flags_have_labelled_type():
    # K2-type flags on 3 vertices (type = labelled edge): the third vertex can see each
    # root independently, so there are 4 flags (no symmetry because the roots are labelled).
    assert len(generate_flags(EDGE, 3, 2)) == 4
    # Single-vertex type, flags on 3 vertices: graphs on 3 vertices with one root.
    # Rooted P3s: centre rooted or end rooted, so 2 rooted P3s, plus K3, empty, rooted
    # K2 + K1 (root in the edge / root isolated) -> 6.
    assert len(generate_flags(Flag(1, 0, ()), 3, 2)) == 6


def test_canonical_and_automorphisms():
    P3a = from_edges(3, [(0, 1, 1), (1, 2, 1)])
    P3b = from_edges(3, [(0, 2, 1), (2, 1, 1)])
    assert canonical(P3a) == canonical(P3b)
    assert automorphism_count(P3a) == 2
    assert automorphism_count(K3) == 6
    # As 1-flags, centre-rooted and end-rooted P3 differ.
    centre = from_edges(3, [(0, 1, 1), (0, 2, 1)], s=1)
    end = from_edges(3, [(0, 1, 1), (1, 2, 1)], s=1)
    assert canonical(centre) != canonical(end)


# ---------------------------------------------------------------------------
# Worked examples from the thesis
# ---------------------------------------------------------------------------

def test_thesis_average_examples():
    # Root at the centre of P3 averages to (1/3) P3; rooted K3 averages to K3.
    centre = from_edges(3, [(0, 1, 1), (0, 2, 1)], s=1)
    assert average_coefficient(centre) == Fraction(1, 3)
    assert average(combo((centre, 1))) == {canonical(Flag(3, 0, centre.cols)): Fraction(1, 3)}
    rooted_k3 = Flag(3, 1, K3.cols)
    assert average_coefficient(rooted_k3) == 1


def test_thesis_product_examples():
    root = Flag(1, 0, ())
    e = from_edges(2, [(0, 1, 1)], s=1)      # root adjacent to the other vertex
    ne = from_edges(2, [], s=1)              # root not adjacent
    # e * e = (both adjacent, edge between) + (both adjacent, no edge)
    sq = multiply(combo((e, 1)), combo((e, 1)), 2)
    expected = combo(
        (from_edges(3, [(0, 1, 1), (0, 2, 1), (1, 2, 1)], s=1), 1),
        (from_edges(3, [(0, 1, 1), (0, 2, 1)], s=1), 1),
    )
    assert sq == expected
    # e * ne = 1/2 (rooted P3 end-ish with edge) + 1/2 (...) : two flags, coefficient 1/2 each
    mixed = multiply(combo((e, 1)), combo((ne, 1)), 2)
    assert sorted(mixed.values()) == [Fraction(1, 2), Fraction(1, 2)]
    assert len(generate_flags(root, 3, 2)) == 6


def test_thesis_p3_sdp_example():
    """Ch. I example: v = (root-adjacent, root-non-adjacent), and
    [[v^T A v]] = A11 K3 + (A11/3 + 2A12/3) P3 + (2A12/3 + A22/3) (K2+K1) + A22 E3."""
    root = Flag(1, 0, ())
    flags = list(generate_flags(root, 2, 2))
    e = canonical(from_edges(2, [(0, 1, 1)], s=1))
    ne = canonical(from_edges(2, [], s=1))
    order = [flags.index(e), flags.index(ne)]
    P3 = canonical(from_edges(3, [(0, 1, 1), (1, 2, 1)]))
    K2K1 = canonical(from_edges(3, [(0, 1, 1)]))
    expected = {
        canonical(K3): [[1, 0], [0, 0]],
        P3: [[Fraction(1, 3), Fraction(1, 3)], [Fraction(1, 3), 0]],
        K2K1: [[0, Fraction(1, 3)], [Fraction(1, 3), Fraction(1, 3)]],
        canonical(EMPTY3): [[0, 0], [0, 1]],
    }
    for G in generate(3, 2):
        Y = Y_matrix(G, root, flags)
        Yo = [[Y[order[i]][order[j]] for j in range(2)] for i in range(2)]
        assert Yo == expected[G], G


# ---------------------------------------------------------------------------
# Exact identities on step graphons
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("seed", range(4))
def test_chain_rule_on_step_graphons(seed):
    rng = random.Random(seed)
    W = random_step_graphon(rng, 3, 2)
    for H in generate(3, 2):
        lifted = lift(combo((H, 1)), 5, 2)
        assert W.evaluate(lifted) == W.density(H)


def test_chain_rule_on_finite_graphs():
    rng = random.Random(7)
    for _ in range(5):
        G = from_edges(7, [(i, j, 1) for i in range(7) for j in range(i + 1, 7) if rng.random() < 0.5])
        for H in generate(3, 2):
            lhs = flag_density(H, G)
            rhs = sum(flag_density(H, F) * flag_density(F, G) for F in generate(5, 2))
            assert lhs == rhs


@pytest.mark.parametrize("seed", range(4))
def test_product_identity_on_step_graphons(seed):
    rng = random.Random(100 + seed)
    W = random_step_graphon(rng, 3, 2)
    sigma = generate(1, 2)[0]
    flags = generate_flags(sigma, 3, 2)
    for F1 in flags:
        for F2 in flags:
            prod = multiply(combo((F1, 1)), combo((F2, 1)), 2)
            for roots, _ in W.root_placements(1):
                assert W.evaluate(prod, roots) == W.flag_density(F1, roots) * W.flag_density(F2, roots)


@pytest.mark.parametrize("seed", range(4))
def test_average_identity_on_step_graphons(seed):
    rng = random.Random(200 + seed)
    W = random_step_graphon(rng, 3, 3)
    for sigma in generate(2, 3):
        for F in generate_flags(sigma, 4, 3):
            expected = sum(p * W.flag_density(F, roots) for roots, p in W.root_placements(2))
            assert W.evaluate(average(combo((F, 1)))) == expected


@pytest.mark.parametrize("n,colors,seed", [(4, 2, 0), (5, 2, 1), (4, 3, 2), (5, 3, 3)])
def test_Y_identity_on_step_graphons(n, colors, seed):
    """The identity the whole method rests on: sum_G p(G) Y_{G,sigma} = E_theta[v v^T]."""
    rng = random.Random(300 + seed)
    W = random_step_graphon(rng, 3, colors)
    graphs = generate(n, colors)
    dist = W.graph_distribution(n)
    assert sum(dist.values()) == 1
    assert set(dist) <= set(graphs)
    for sigma, flags in sdp_blocks(n, colors, include_empty_type=True):
        size = len(flags)
        M = [[Fraction(0)] * size for _ in range(size)]
        for G in graphs:
            pG = dist.get(G, Fraction(0))
            if pG:
                Y = Y_matrix(G, sigma, flags)
                for i in range(size):
                    for j in range(size):
                        M[i][j] += pG * Y[i][j]
        E = [[Fraction(0)] * size for _ in range(size)]
        for roots, p in W.root_placements(sigma.k):
            if restrict(W._sample_flag(roots, 0), range(sigma.k), 0).cols != sigma.cols:
                continue
            v = [W.flag_density(F, roots) for F in flags]
            for i in range(size):
                for j in range(size):
                    E[i][j] += p * v[i] * v[j]
        assert M == E, sigma
