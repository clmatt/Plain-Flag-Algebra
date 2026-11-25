"""The quad constraints (nbhd_quads, neighbourhoods.quad_multipliers): the flag encoding,
the offsets, and validity at roots whose neighbourhood is a real listed graph."""
import gzip
import os
import random
from fractions import Fraction
from itertools import combinations

import pytest

from flagcheck import nbhd_quads as nq
from flagcheck import neighbourhoods as nb
from flagcheck.graphs import from_edges
from flagcheck.limits import StepGraphon
from flagcheck.ramsey import BLUE, RED
from flagcheck.ramsey_run import build_problem, degree_window
from flagcheck.spec import Add, Lit, Scale

DATA = os.path.expanduser("~/data/ramsey")
HAVE_LISTS = all(nq.list_path(DATA, *t) for t in nq.LIST_COUNTS)


def _counts(adj, verts):
    """c in Z^11 of the graph induced on verts (adj: set of frozenset edges), by brute force."""
    c = [0] * 11
    for S in combinations(verts, 4):
        E = [(i, j) for i, j in combinations(range(4), 2) if frozenset((S[i], S[j])) in adj]
        c[nq.class_of(E)] += 1
    return c


def _g6_edges(line):
    from flagcheck.edge_tables import parse_g6

    n, edges = parse_g6(line)
    return n, {frozenset(e) for e in edges}


def test_classes_and_admissibility():
    assert [nq.class_of(E) for E in nq.REPS] == list(range(11))
    assert len(set(nq.CLASSES)) == 11
    assert nq.admissible_classes(3, 6) == [0, 1, 2, 3, 5, 6, 7]          # triangle-free
    assert nq.admissible_classes(4, 5) == list(range(10))                 # no K4
    assert nq.admissible_classes(4, 4) == list(range(1, 10))              # no K4, no empty 4-set


@pytest.mark.parametrize("seed", range(3))
@pytest.mark.parametrize("rc", [RED, BLUE])
def test_quad_flag_density_is_24_c_over_M4(rc, seed):
    """At a root of the balanced blow-up of an M-vertex graph, Q[rc,h] has density
    24 c_h / M^4, c = counts of the root's rc-neighbourhood (green clone pairs excluded)."""
    rng = random.Random(seed)
    M = 9
    cols = {}
    for x, y in combinations(range(M), 2):
        cols[(x, y)] = rng.choice([RED, BLUE])
    G = from_edges(M, [(x, y, c) for (x, y), c in cols.items()])
    W = StepGraphon.blowup(G)
    nbhd = [x for x in range(1, M) if cols[(0, x)] == rc]
    adj = {frozenset(e) for e, c in cols.items() if c == RED}
    c = _counts(adj, nbhd)
    for h in range(11):
        assert W.flag_density(nb.quad_flag(rc, h), (0,)) == Fraction(24 * c[h], M ** 4)


def _cut_value(g, quad_dens, W):
    """Value of g = sum a_h Q[rc,h] - B(R) at root 0, with the Q densities supplied."""
    from tests.test_neighbourhoods import _value_at

    total = Fraction(0)
    assert isinstance(g, Add)
    for t in g.args:
        if isinstance(t, Scale) and isinstance(t.arg, Lit) and t.arg.flag.k == 5:
            total += t.q * quad_dens[t.arg.flag]
        else:
            total += _value_at(W, t, (0,))
    return total


def _root_graph(N, red_nbhd_adj, s_red, blue_nbhd_adj, s_blue):
    """Root 0 red to 1..s_red (inducing red_nbhd_adj), blue to the rest (inducing
    blue_nbhd_adj on s_red+1..N-1, relabelled from 0); red between the two parts."""
    assert 1 + s_red + s_blue == N
    edges = []
    for x in range(1, N):
        edges.append((0, x, RED if x <= s_red else BLUE))
    for x, y in combinations(range(1, N), 2):
        if y <= s_red:
            col = RED if frozenset((x - 1, y - 1)) in red_nbhd_adj else BLUE
        elif x > s_red:
            col = RED if frozenset((x - 1 - s_red, y - 1 - s_red)) in blue_nbhd_adj else BLUE
        else:
            col = RED if (x + y) % 2 else BLUE
        edges.append((x, y, col))
    return from_edges(N, edges)


@pytest.mark.skipif(not HAVE_LISTS, reason="McKay's lists absent")
def test_quad_cuts_hold_at_real_neighbourhoods():
    """R(4,6), N = 40: every stored cut is >= 0 at a root whose red (blue) neighbourhood is a
    listed (3,6)- ((4,5)-) graph of the right order, and equals 24 (a.c - b_d) / N^4 there."""
    N = 40
    lo, hi = degree_window(4, 6, N)
    mults, facts = nb.quad_multipliers(4, 6, N, range(lo, hi + 1))
    cuts = [c for c in nq.load_cuts() if (c["k"], c["l"], c["N"]) == (4, 6, N)]
    assert len(mults) == len(cuts) > 0 and facts
    from flagcheck.spec import Mul
    from tests.test_neighbourhoods import _value_at
    with open(os.path.join(DATA, "r45_24.g6"), "rb") as fh:
        blue24 = [ln for i, ln in enumerate(fh) if i in (0, 1000, 176183, 352365)]
    for d in range(lo, hi + 1):
        path = nq.list_path(DATA, 3, 6, d)
        with (gzip.open if path.endswith(".gz") else open)(path, "rb") as fh:
            reds = [ln for i, ln in enumerate(fh) if i < 3]
        for red_line in reds:
            s_red, red_adj = _g6_edges(red_line)
            blue_lines = blue24 if d == 15 else [None]
            for blue_line in blue_lines:
                s_blue = N - 1 - d
                blue_adj = _g6_edges(blue_line)[1] if blue_line else set()
                G = _root_graph(N, red_adj, s_red, blue_adj, s_blue)
                W = StepGraphon.blowup(G)
                c = {RED: _counts(red_adj, range(s_red)), BLUE: _counts(blue_adj, range(s_blue))}
                for (name, g), cut in zip(mults, cuts):
                    rc = RED if cut["colour"] == "red" else BLUE
                    dens = {nb.quad_flag(rc, h): Fraction(24 * c[rc][h], N ** 4) for h in range(11)}
                    if "degree" in cut:                 # L_e(R) (a.Q - 24 b / N^4)
                        assert isinstance(g, Mul)
                        L, inner = g.args
                        val = _value_at(W, L, (0,)) * _cut_value(inner, dens, W)
                        e = cut["degree"]
                        want = (Fraction(24, N ** 4) * (sum(a * x for a, x in zip(cut["a"], c[rc]))
                                                        - cut["offsets"][e][0]) if d == e else 0)
                    else:
                        val = _cut_value(g, dens, W)
                        want = Fraction(24, N ** 4) * (sum(a * x for a, x in zip(cut["a"], c[rc]))
                                                       - cut["offsets"][d][0])
                    assert val == want, (name, d)
                    if rc == RED or blue_line:          # a real listed neighbourhood
                        assert val >= 0, (name, d)


A_TEST = [Fraction(x) for x in (1, -1, 2, 0, 0, Fraction(1, 3), -1, 1, 0, 0, 0)]


def _loaded(raw):
    """A make_cuts entry in load_cuts form."""
    out = {**raw, "a": [Fraction(x) for x in raw["a"]],
           "offsets": {int(d): (Fraction(v["b"]), v["source"]) for d, v in raw["offsets"].items()}}
    return out


@pytest.mark.skipif(not os.path.exists(os.path.join(DATA, "r36_17.g6")), reason="lists absent")
def test_check_cuts_catches_a_wrong_offset():
    cut = _loaded(nq.make_cuts([("red", A_TEST, 17)], 4, 6, 40, range(15, 18), DATA)[0])
    b, src = cut["offsets"][17]
    assert src == "list"
    nq.check_cuts([cut], DATA)
    with pytest.raises(ValueError):
        nq.check_cuts([dict(cut, offsets={17: (b + 1, src)})], DATA)


def test_quad_multipliers_reject_a_wrong_trivial_offset(monkeypatch):
    cut = _loaded(nq.make_cuts([("blue", A_TEST)], 4, 6, 40, range(15, 18), "/nonexistent")[0])
    assert all(s == "trivial" for _, s in cut["offsets"].values())
    monkeypatch.setattr(nq, "load_cuts", lambda: [cut])
    nb.quad_multipliers(4, 6, 40, range(15, 18))
    cut["offsets"][16] = (cut["offsets"][16][0] + 1, "trivial")
    with pytest.raises(ValueError):
        nb.quad_multipliers(4, 6, 40, range(15, 18))


def test_quad_extra_builds():
    p = build_problem(4, 6, 40, 7, "local", False, ("quad",))
    assert any("McKay's list" in a for a in p.assumptions)
    assert len(p.ineqs) > 0


@pytest.mark.slow
@pytest.mark.skipif(not HAVE_LISTS, reason="McKay's lists absent")
def test_stored_offsets_recomputed_from_lists():
    nq.check_cuts(nq.load_cuts(), DATA)


@pytest.mark.parametrize("rc", [RED, BLUE])
def test_quad_flags_are_the_eleven_classes(rc):
    """Regression: a representative edge written as (j, i) with j > i was once silently
    dropped (C4 became P4), so check the flag's leaf graph itself, not just the class list."""
    from flagcheck.graphs import canonical

    flags = [canonical(nb.quad_flag(rc, h)) for h in range(11)]
    assert len(set(flags)) == 11
    for h in range(11):
        F = nb.quad_flag(rc, h)
        E = [(i, j) for i, j in combinations(range(4), 2) if F.color(i + 1, j + 1) == RED]
        assert nq.class_of(E) == h
        assert all(F.color(0, x) == rc for x in range(1, 5))


@pytest.mark.skipif(not HAVE_LISTS, reason="McKay's lists absent")
def test_localised_quad_cut_is_zero_off_its_degree_and_valid_on_it():
    """L_d(R)(a.Q - 24 b/N^4) at real roots of every window degree: 0 unless the degree is
    d, and 24 (a.c - b)/N^4 >= 0 there (red neighbourhood a listed (3,6)-graph)."""
    from flagcheck.spec import Mul
    from tests.test_neighbourhoods import _value_at

    N = 40
    window = list(range(15, 18))
    a = [Fraction(x) for x in (1, -1, 2, 0, 0, Fraction(1, 3), -1, 1, 0, 0, 0)]
    for d in window:
        cut = nq.make_cuts([("red", a, d)], 4, 6, N, window, DATA)[0]
        cut = {**cut, "a": a, "offsets": {d: (Fraction(cut["offsets"][d]["b"]), "list")}}
        g, facts = nb.quad_cut_expr(4, 6, N, window, cut)
        assert facts and isinstance(g, Mul) and nb.expr_size(g) == 7
        L, inner = g.args
        b = cut["offsets"][d][0]
        for e in window:
            path = nq.list_path(DATA, 3, 6, e)
            with (gzip.open if path.endswith(".gz") else open)(path, "rb") as fh:
                line = next(ln for ln in fh if ln.strip())
            s, adj = _g6_edges(line)
            W = StepGraphon.blowup(_root_graph(N, adj, s, set(), N - 1 - s))
            c = _counts(adj, range(s))
            dens = {nb.quad_flag(RED, h): Fraction(24 * c[h], N ** 4) for h in range(11)}
            val = _value_at(W, L, (0,)) * _cut_value(inner, dens, W)
            want = Fraction(24, N ** 4) * (sum(x * y for x, y in zip(a, c)) - b) if e == d else 0
            assert val == want
            assert val >= 0
