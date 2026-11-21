"""Soundness and end-to-end tests of the Ramsey driver (flagcheck.ramsey_run)."""

import shutil
from fractions import Fraction

import pytest

from flagcheck import cppbackend
from flagcheck.limits import StepGraphon
from flagcheck.ramsey import KNOWN_RAMSEY_GRAPHS, red_degrees
from flagcheck.ramsey_run import build_problem, degree_window, ramsey_upper, run
from flagcheck.spec import evaluate


@pytest.mark.parametrize("kl,n", [((3, 3), 4), ((3, 4), 5), ((3, 5), 4)])
def test_every_constraint_holds_on_a_real_critical_graph(kl, n):
    """For N = R(k,l) - 1 a real graph exists, so its balanced blow-up must satisfy every
    equality exactly, and every localising multiplier must be >= 0 at every root."""
    k, l = kl
    G = KNOWN_RAMSEY_GRAPHS[kl]()
    N = G.k
    lo, hi = degree_window(k, l, N)
    assert lo <= min(red_degrees(G)) and max(red_degrees(G)) <= hi
    W = StepGraphon.blowup(G)
    p = build_problem(k, l, N, n, "all")
    assert p.eqs
    for expr, rhs in p.eqs:
        assert W.evaluate(evaluate(expr, 3, p.zeros)) == rhs
    for m, g in p.localizing:
        gc = evaluate(g, 3, p.zeros)
        assert all(W.evaluate(gc, (v,)) >= 0 for v in range(N))


def test_degree_window_values():
    assert ramsey_upper(2, 5) == 5 and ramsey_upper(5, 4) == 25
    assert degree_window(4, 5, 24) == (24 - 18, 13)
    assert degree_window(4, 6, 39) == (39 - 25, 17)   # thesis: red degree in [14, 17]


@pytest.mark.skipif(not cppbackend.available() or shutil.which("csdp") is None, reason="needs flagsdp + CSDP")
@pytest.mark.parametrize("degrees", ["none", "avg", "local", "psd", "all"])
def test_driver_end_to_end_small(degrees):
    # R(3,4) <= 9 must be proved with every option, and R(3,4) <= 8 never.
    _, _, res9, _ = run(3, 4, 9, 5, degrees)
    _, _, res8, _ = run(3, 4, 8, 5, degrees)
    assert res9.bound > 0
    assert res8.bound <= 0


ALL_EXTRAS = ("green", "codeg", "tcodeg", "edges", "pedges", "nstats", "lpsd")


@pytest.mark.skipif(not cppbackend.available() or shutil.which("csdp") is None, reason="needs flagsdp + CSDP")
@pytest.mark.parametrize("extras", [("green",), ("codeg",), ("edges",), ALL_EXTRAS])
def test_driver_with_extras_small(extras):
    _, _, res9, _ = run(3, 4, 9, 5, "none", extras=extras)
    _, _, res8, _ = run(3, 4, 8, 5, "none", extras=extras)
    assert res9.bound > 0
    assert res8.bound <= 0


@pytest.mark.skipif(not cppbackend.available() or shutil.which("csdp") is None or shutil.which("sdpa") is None,
                    reason="needs flagsdp + CSDP + SDPA")
@pytest.mark.parametrize("args", [(3, 4, 9, 5, "all", ()), (3, 4, 8, 5, "all", ()),
                                  (3, 4, 9, 5, "none", ALL_EXTRAS), (4, 4, 17, 6, "none", ())])
def test_sdpa_agrees_with_csdp(args):
    """Same input file, same verdict from the exact checker, and objectives that agree to
    solver accuracy (SDPA's standard form is CSDP's with primal and dual swapped)."""
    k, l, N, n, deg, extras = args
    _, sc, rc, _ = run(k, l, N, n, deg, extras=extras, solver="csdp")
    _, ss, rs, _ = run(k, l, N, n, deg, extras=extras, solver="sdpa")
    assert (rc.bound > 0) == (rs.bound > 0)
    assert abs(sc.primal_objective - ss.primal_objective) < 1e-6
    assert len(ss.densities) == len(sc.densities)


@pytest.mark.skipif(not cppbackend.available() or shutil.which("csdp") is None, reason="needs flagsdp + CSDP")
@pytest.mark.parametrize("kln", [(3, 5, 13, 5), (4, 4, 17, 5), (4, 4, 17, 6)])
def test_extras_never_refute_a_real_graph(kln):
    """C13 and Paley(17) exist (and are tight for several of the bounds), so no option may
    ever produce a positive certificate here."""
    k, l, N, n = kln
    _, sol, res, _ = run(k, l, N, n, "all", extras=ALL_EXTRAS)
    assert res.bound <= 0, sol.primal_objective
