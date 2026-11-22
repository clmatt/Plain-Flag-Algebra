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
    """Values of g at every root placement whose roots induce g's type (elsewhere every
    flag of that type, hence g, is trivially 0)."""
    from flagcheck.algebra import identity_flag

    sigma = nb.sigma_of(g)
    assert sigma.k == s
    ident = {identity_flag(sigma): Fraction(1)}
    gc = evaluate(g, 3, ())
    vals = [W.evaluate(gc, parts) for parts, _ in W.root_placements(s) if W.evaluate(ident, parts) == 1]
    assert vals, "no root placement of this type"
    return vals


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
    pedges, _ = nb.pair_edge_multipliers(k, l, N, ramsey_upper)
    for name, g in codeg + edges + pedges:
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


def test_pair_edge_multipliers_shape_and_tightness():
    """R(4,6): the five non-trivial pair families, with the right Ramsey classes. On
    Paley(17), the flags D and T at concrete pairs equal |S|/N and 2 e(S)/N^2 computed
    directly from the graph, and every pair multiplier is >= 0 at every pair."""
    mults, facts = nb.pair_edge_multipliers(4, 6, 40, ramsey_upper)
    tags = {name.split("]")[0] + "]" for name, _ in mults}
    assert tags == {"pair-edges[r:rb]", "pair-edges[r:bb]", "pair-edges[b:rr]", "pair-edges[b:rb]",
                    "pair-edges[b:bb]"}
    for fact in ("R(3,5) <= 14", "R(4,5) <= 25", "R(3,6) <= 18", "R(4,4) <= 18", "e(3,5,13) in [26,26]"):
        assert fact in facts
    G = paley(17)
    N = G.k
    W = StepGraphon.blowup(G)
    for u, v, sigma in ((0, 1, nb.RED_PAIR), (0, 3, nb.BLUE_PAIR)):      # 1 is a square mod 17, 3 is not
        assert G.color(u, v) == sigma.color(0, 1)
        for pattern in ((RED, RED), (RED, BLUE), (BLUE, RED), (BLUE, BLUE)):
            S = [x for x in range(N) if x not in (u, v) and (G.color(x, u), G.color(x, v)) == pattern]
            e = sum(1 for i in S for j in S if i < j and G.color(i, j) == RED)
            assert W.flag_density(nb.extension(sigma, pattern), (u, v)) == Fraction(len(S), N)
            assert W.flag_density(nb.pair_extension(sigma, pattern, RED), (u, v)) == Fraction(2 * e, N * N)
    pm, _ = nb.pair_edge_multipliers(4, 4, 17, ramsey_upper)
    assert pm
    for name, g in pm:
        assert min(_pointwise_values(W, g, 2)) >= 0, name


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
    p = build_problem(k, l, G.k, n, "none", extras=("green", "codeg", "tcodeg", "edges", "pedges", "nstats", "lpsd"))
    assert p.ineqs and p.assumptions
    for expr, rhs in p.eqs:
        assert W.evaluate(evaluate(expr, 3, p.zeros)) == rhs
    for expr, rhs in p.ineqs:
        assert W.evaluate(evaluate(expr, 3, p.zeros)) <= rhs
    for m, g in p.localizing:
        s = nb.sigma_of(g).k
        assert min(_pointwise_values(W, g, s)) >= 0


# ---------------------------------------------------------------------------
# Three roots
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("kl,G", CASES, ids=lambda x: str(x) if isinstance(x, tuple) else str(x.k))
def test_triple_codegree_holds_pointwise(kl, G):
    k, l = kl
    W = StepGraphon.blowup(G)
    mults, _ = nb.triple_codegree_multipliers(k, l, G.k, ramsey_upper)
    checked = 0
    for name, g in mults:
        try:
            vals = _pointwise_values(W, g, 3)
        except AssertionError:          # G has no root triple of this type (C5: no blue triangle)
            continue
        assert min(vals) >= 0, name
        checked += 1
    assert checked or not mults         # (3,3), N = 5: every triple bound is trivial


def test_triple_codegree_semantics_and_tightness_on_paley17():
    """At concrete root triangles of Paley(17) the extension density is |S|/N computed from
    the graph, |S| never exceeds the bound, and some bound is attained."""
    from itertools import combinations as combs

    G = paley(17)
    N = G.k
    W = StepGraphon.blowup(G)
    mults, facts = nb.triple_codegree_multipliers(4, 4, N, ramsey_upper)
    assert "R(2,3) <= 3" in facts
    attained = False
    for name, g in mults:
        sigma = nb.sigma_of(g)
        tag, rest = name[len("tcodeg["):].split(":")
        pattern = tuple(RED if ch == "r" else BLUE for ch in rest[:3])
        c = int(rest.split("<=")[1].rstrip("]"))
        for roots in combs(range(N), 3):
            for perm in ((0, 1, 2), (0, 2, 1), (1, 0, 2), (1, 2, 0), (2, 0, 1), (2, 1, 0)):
                r = tuple(roots[p] for p in perm)
                if any(G.color(r[i], r[j]) != sigma.color(i, j) for i, j in combs(range(3), 2)):
                    continue
                S = [x for x in range(N) if x not in r and all(G.color(x, r[i]) == pattern[i] for i in range(3))]
                assert W.flag_density(nb.extension(sigma, pattern), r) == Fraction(len(S), N)
                assert len(S) <= c, name
                attained |= len(S) == c
    assert attained


def test_triple_codegree_families():
    m, facts = nb.triple_codegree_multipliers(4, 6, 40, ramsey_upper)
    names = {nm for nm, _ in m}
    # common blue neighbourhood of a blue triangle is a (4,3)-graph: at most 8 vertices
    assert "tcodeg[bbb:bbb<=8]" in names
    # x red to a red edge and blue to a third root: (2,5)-graph, beats the pair bound 5
    assert "tcodeg[rrr:rrb<=4]" in names
    # implied by pairs, hence absent: x blue to all of a red triangle ((4,5), 24 = pair bound)
    assert not any(nm.startswith("tcodeg[rrr:bbb") for nm in names)
    assert "R(3,4) <= 9" in facts and "R(2,5) <= 5" in facts


# ---------------------------------------------------------------------------
# Neighbourhood structure (nstats)
# ---------------------------------------------------------------------------

import os
from itertools import combinations as _combs

from flagcheck import nbhd_stats
from flagcheck.graphs import from_edges

DATA = os.path.expanduser("~/data/ramsey")


def _value_at(W, e, root):
    """Exact value of an expression at one root of a blow-up: rooted densities multiply
    pointwise, so this is independent of the flag-product machinery."""
    if isinstance(e, Lit):
        return W.flag_density(e.flag, root)
    if isinstance(e, dict):
        return sum(v * W.flag_density(F, root) for F, v in e.items())
    if isinstance(e, Add):
        return sum(_value_at(W, a, root) for a in e.args)
    if isinstance(e, Scale):
        return e.q * _value_at(W, e.arg, root)
    if isinstance(e, Mul):
        out = Fraction(1)
        for a in e.args:
            out *= _value_at(W, a, root)
        return out
    raise TypeError(e)


def _stats(G, verts):
    """(e, p, t) of the red graph induced on verts, by brute force."""
    cnt = [0, 0, 0, 0]
    for T in _combs(verts, 3):
        cnt[sum(G.color(x, y) == RED for x, y in _combs(T, 2))] += 1
    e = sum(G.color(x, y) == RED for x, y in _combs(verts, 2))
    return e, cnt[2], cnt[3]


def _graph_from_g6(line):
    from flagcheck.edge_tables import parse_g6

    n, edges = parse_g6(line)
    E = set(edges)
    return from_edges(n, [(i, j, RED if (i, j) in E else BLUE) for i, j in _combs(range(n), 2)])


def test_triple_table_matches_published_counts_and_edge_tables():
    T = nbhd_stats.table()
    for (a, b, s), n in nbhd_stats.EXPECTED_COUNTS.items():
        assert T[(a, b)][s]["count"] == n
    for (a, b), d in T.items():
        for s, v in d.items():
            es = [p[0] for p in v["points"]]
            assert (min(es), max(es)) == EDGE_BOUNDS[(a, b)][s], (a, b, s)


@pytest.mark.parametrize("ab", [(3, 4), (3, 5), (3, 6), (4, 4)])
def test_triple_table_complete_on_small_orders_by_brute_force(ab):
    """Independent of McKay's files: every labelled (a,b)-graph on s <= 6 vertices."""
    a, b = ab
    for s in range(3, 7):
        pairs = list(_combs(range(s), 2))
        seen = set()
        for mask in range(1 << len(pairs)):
            edges = [p for i, p in enumerate(pairs) if mask >> i & 1]
            if is_ab_graph(s, edges, a, b):
                G = from_edges(s, [(i, j, RED if (i, j) in set(edges) else BLUE) for i, j in pairs])
                seen.add(_stats(G, range(s)))
        assert seen == nbhd_stats.triple_points(a, b, s), (a, b, s)


def test_complement_points():
    """(4,3)-graphs are complements of (3,4)-graphs: checked on the Paley(17) blue
    neighbourhoods, which are (4,3)-graphs on 8 vertices."""
    G = paley(17)
    for v in range(17):
        blue = [x for x in range(17) if x != v and G.color(v, x) == BLUE]
        red = [x for x in range(17) if x != v and G.color(v, x) == RED]
        assert _stats(G, blue) in nbhd_stats.triple_points(4, 3, 8)
        assert _stats(G, red) in nbhd_stats.triple_points(3, 4, 8)


def test_stat_flags_mean_what_they_say():
    G = paley(17)
    N = G.k
    W = StepGraphon.blowup(G)
    for v in (0, 5):
        for rc in (RED, BLUE):
            S = [x for x in range(N) if x != v and G.color(v, x) == rc]
            e, p, t = _stats(G, S)
            assert W.flag_density(nb.double_extension(rc, RED), (v,)) == Fraction(2 * e, N * N)
            assert W.flag_density(nb.triple_extension(rc, 2), (v,)) == Fraction(6 * p, N ** 3)
            assert W.flag_density(nb.triple_extension(rc, 3), (v,)) == Fraction(6 * t, N ** 3)


def _real_graphs():
    yield (4, 4), paley(17)
    yield (3, 5), circulant(13, [1, 5])
    if os.path.exists(os.path.join(DATA, "r45_24.g6")):
        with open(os.path.join(DATA, "r45_24.g6"), "rb") as fh:
            for i, line in enumerate(fh):
                if i in (0, 176183, 352365):
                    yield (4, 5), _graph_from_g6(line)
    for s in (16, 15):
        path = os.path.join(DATA, f"r44_{s}.g6")
        if os.path.exists(path):
            with open(path, "rb") as fh:
                yield (4, 4), _graph_from_g6(fh.readline())


@pytest.mark.parametrize("kl,G", list(_real_graphs()), ids=lambda x: str(x) if isinstance(x, tuple) else str(x.k))
def test_stat_multipliers_hold_at_every_vertex(kl, G):
    k, l = kl
    N = G.k
    assert is_ramsey_graph(G, k, l)
    lo, hi = degree_window(k, l, N)
    mults, facts = nb.stat_multipliers(k, l, N, range(lo, hi + 1))
    W = StepGraphon.blowup(G)
    tight = 0
    for name, g in mults:
        vals = [_value_at(W, g, (v,)) for v in range(N)]
        assert min(vals) >= 0, name
        tight += min(vals) == 0
    assert tight, "some hull line should be attained by a real graph"
    assert facts
    if kl == (4, 5):     # neighbourhoods of a (4,5)-graph are (3,5)- and (4,4)-graphs
        assert any("(3,5)-graphs" in f for f in facts) and any("(4,4)-graphs" in f for f in facts)


def test_stat_multiplier_detects_a_wrong_point():
    """A table holding only a point with a different path count than the actual
    neighbourhood must make some multiplier negative at that vertex: the constraints
    really use the data."""
    G = paley(17)
    N = G.k
    W = StepGraphon.blowup(G)
    real = nbhd_stats.triple_points(3, 4, 8)
    S = [x for x in range(1, N) if G.color(0, x) == RED]
    mine = _stats(G, S)
    others = {next(q for q in real if q[1] != mine[1])}
    saved = nbhd_stats.table()[(3, 4)][8]["points"]
    try:
        nbhd_stats.table()[(3, 4)][8]["points"] = others
        mults, _ = nb.stat_multipliers(4, 4, N, [8])
        assert min(_value_at(W, g, (0,)) for _, g in mults) < 0
    finally:
        nbhd_stats.table()[(3, 4)][8]["points"] = saved


# ---------------------------------------------------------------------------
# mrtri: triangles against edges in (4,5)-neighbourhoods
# ---------------------------------------------------------------------------

def _extreme_45(orders):
    """{s: [graph6 lines]} of r45extreme.tar.gz for the given orders (all graphs with the
    file's (s, e))."""
    import tarfile

    out = {}
    with tarfile.open(os.path.join(DATA, "r45extreme.tar.gz")) as tar:
        for m in tar.getmembers():
            if m.isfile():
                for line in tar.extractfile(m).read().splitlines():
                    if line.strip() and line.strip()[0] - 63 in orders:
                        out.setdefault(line.strip()[0] - 63, []).append(line)
    return out


@pytest.mark.skipif(not os.path.exists(os.path.join(DATA, "r45extreme.tar.gz")), reason="McKay's lists absent")
def test_triangle_tables_against_every_listed_graph():
    """EXTREME_TRIANGLES is exactly what the archive gives, and no listed (4,5)-graph on 22-24
    vertices violates MR5 Table IV (a check of the transcription)."""
    seen = {}
    for s, lines in _extreme_45((22, 23)).items():
        for line in lines:
            e, _, t = nbhd_stats.triple_stats(*nbhd_stats.adjacency(line))
            a, b = seen.setdefault(s, {}).get(e, (t, t))
            seen[s][e] = (min(a, t), max(b, t))
    assert seen == nbhd_stats.EXTREME_TRIANGLES
    for s, table in nbhd_stats.MR5_TRIANGLES.items():
        pts = nbhd_stats.triple_points(4, 5, s) or set()
        pts |= {(e, None, t) for e, rng in seen.get(s, {}).items() for t in rng}
        assert pts
        for e, _, t in pts:
            assert table[e][0] <= t <= table[e][1], (s, e, t)
    for s in (22, 23, 24):
        bounds = nbhd_stats.triangle_bounds(s)
        assert set(bounds) <= set(range(edge_bounds(4, 5, s)[0], edge_bounds(4, 5, s)[1] + 1))
        for e, rng in seen.get(s, {}).items():
            assert bounds[e] == rng


def _rooted(H, N, red_extra):
    """A graph on N vertices whose root 0 is blue to a copy of H (red edges of H kept) and
    red to the other N-1-|H| vertices; red_extra decides colours among the rest."""
    s = H.k
    edges = []
    for x in range(1, N):
        edges.append((0, x, BLUE if x <= s else RED))
    for x, y in _combs(range(1, N), 2):
        if y <= s:
            edges.append((x, y, H.color(x - 1, y - 1)))
        else:
            edges.append((x, y, RED if red_extra(x, y) else BLUE))
    return from_edges(N, edges)


def _sample_45(s):
    if s == 24:
        with open(os.path.join(DATA, "r45_24.g6"), "rb") as fh:
            lines = [ln for i, ln in enumerate(fh) if i in (0, 1000, 176183, 352365)]
    else:
        lines = _extreme_45((s,))[s]
        lines = lines[:2] + lines[-2:]
    return [_graph_from_g6(ln) for ln in lines]


@pytest.mark.skipif(not os.path.exists(os.path.join(DATA, "r45_24.g6")), reason="McKay's lists absent")
@pytest.mark.parametrize("s", [22, 23, 24])
def test_triangle_multipliers_hold_at_real_neighbourhoods(s):
    """R(4,6), N = 40: a root of red degree 39 - s whose blue neighbourhood is a real
    (4,5)-graph on s vertices satisfies every mrtri multiplier."""
    N = 40
    lo, hi = degree_window(4, 6, N)
    mults, facts = nb.triangle_multipliers(4, 6, N, range(lo, hi + 1))
    assert facts and all("blue" in name for name, _ in mults)
    for H in _sample_45(s):
        W = StepGraphon.blowup(_rooted(H, N, lambda x, y: (x + y) % 3 == 0))
        vals = [_value_at(W, g, (0,)) for _, g in mults]
        assert min(vals) >= 0


@pytest.mark.skipif(not os.path.exists(os.path.join(DATA, "r45_24.g6")), reason="McKay's lists absent")
def test_triangle_multiplier_detects_too_few_triangles():
    """Raising the lower bound above a real graph's triangle count must make a multiplier
    negative at its root: the constraints really use the table."""
    N = 40
    lo, hi = degree_window(4, 6, N)
    H = _sample_45(23)[0]
    e, _, t = _stats(H, range(23))
    saved = dict(nbhd_stats.EXTREME_TRIANGLES[23])
    try:
        nbhd_stats.EXTREME_TRIANGLES[23][e] = (t + 1, t + 1)
        mults, _ = nb.triangle_multipliers(4, 6, N, range(lo, hi + 1))
        W = StepGraphon.blowup(_rooted(H, N, lambda x, y: (x + y) % 3 == 0))
        assert min(_value_at(W, g, (0,)) for _, g in mults) < 0
    finally:
        nbhd_stats.EXTREME_TRIANGLES[23] = saved


# ---------------------------------------------------------------------------
# eint: integrality of neighbourhood edge counts
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("kl,G", CASES[1:4] + [((4, 4), paley(17))], ids=lambda x: str(x) if isinstance(x, tuple) else str(x.k))
def test_edge_integrality_counts_and_cuts(kl, G):
    """X = (N^2/2) T at each root is exactly the number of red edges inside the root's red
    (blue) neighbourhood, and every cut g_a >= 0 holds at every root."""
    from itertools import combinations as _c

    from flagcheck.ramsey import one_flag

    k, l = kl
    N = G.k
    W = StepGraphon.blowup(G)
    lo, hi = degree_window(k, l, N)
    mults, facts = nb.edge_integrality_multipliers(k, l, N, range(lo, hi + 1))
    assert mults and facts
    half = Fraction(N * N, 2)
    for v in range(N):
        for color, T in ((RED, nb.double_extension(RED, RED)), (BLUE, nb.double_extension(BLUE, RED))):
            nbhd = [u for u in range(N) if u != v and G.color(u, v) == color]
            count = sum(G.color(x, y) == RED for x, y in _c(nbhd, 2))
            assert half * W.flag_density(T, (v,)) == count
        for name, g in mults:
            assert _value_at(W, g, (v,)) >= 0, (name, v)
    # a cut at a half-integer offset is not valid: it fails at some vertex
    T = nb.double_extension(RED, RED)
    X = Scale(half, Lit(T))
    xs = {half * W.flag_density(T, (v,)) for v in range(N)}
    c = min(xs)
    bad = Mul((Add((X, Scale(-(c - Fraction(1, 2)), one_flag()))),
               Add((X, Scale(-(c + Fraction(1, 2)), one_flag())))))
    assert min(_value_at(W, bad, (v,)) for v in range(N)) < 0
