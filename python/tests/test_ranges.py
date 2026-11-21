"""The certified range of a statistic must contain its value on every real blow-up, and
adding valid constraints must not widen it."""

import shutil
from fractions import Fraction

import pytest

from flagcheck import cppbackend
from flagcheck.graphs import from_color_function
from flagcheck.limits import StepGraphon
from flagcheck.ramsey import BLUE, GREEN, RED
from flagcheck.ranges import STATS, stat_range

pytestmark = pytest.mark.skipif(not cppbackend.available() or shutil.which("csdp") is None,
                                reason="needs flagsdp + CSDP")


def wagner():
    """C8 plus long diagonals: triangle-free, independence number 3, so a (3,4)-graph on 8."""
    red = {(i, (i + d) % 8) for i in range(8) for d in (1, 4)}
    return from_color_function(8, lambda i, j: RED if (i, j) in red or (j, i) in red else BLUE)


def value(H, stat):
    dist = StepGraphon.blowup(H, clone_color=GREEN).graph_distribution(max(F.k for F in STATS[stat]))
    from flagcheck.algebra import lift

    total = Fraction(0)
    for F, c in STATS[stat].items():
        for G, d in lift({F: Fraction(1)}, max(F2.k for F2 in STATS[stat]), 3, ()).items():
            total += c * d * dist.get(G, 0)
    return total


@pytest.mark.parametrize("stat", ["red2", "blue3"])
def test_range_contains_real_graph(stat):
    v = value(wagner(), stat)
    lo, hi, _ = stat_range(3, 4, 8, 5, "all", stat=stat)
    assert lo <= v <= hi, (float(lo), float(v), float(hi))
    lo2, hi2, _ = stat_range(3, 4, 8, 5, "all", extras=("green", "codeg", "edges"), stat=stat)
    assert lo2 <= v <= hi2
    assert lo2 >= lo - Fraction(1, 10 ** 7) and hi2 <= hi + Fraction(1, 10 ** 7)   # never (materially) wider


def expr_value(H, expr):
    """Exact value of a spec expression (s = 0) in the balanced blow-up of H."""
    from flagcheck.algebra import lift
    from flagcheck.spec import evaluate

    c = evaluate(expr, 3, ())
    m = max(F.k for F in c)
    dist = StepGraphon.blowup(H, clone_color=GREEN).graph_distribution(m)
    return sum(v * dist.get(G, 0) for G, v in lift(c, m, 3, ()).items())


@pytest.mark.parametrize("stat,count", [("red2", 12), ("deg3", 8), ("deg2", 0), ("odd", 4), ("blue3", None)])
def test_counts_are_whole_numbers(stat, count):
    from flagcheck.ranges import statistic

    S, scale = statistic(stat, 3, 4, 8)
    c = expr_value(wagner(), S) * scale
    assert c.denominator == 1
    if count is not None:
        assert c == count


def test_propagation_keeps_a_real_graph(tmp_path):
    from flagcheck.ranges import propagate, statistic

    status, bounds = propagate(3, 4, 8, 5, "all", (), None, False, ["red2", "deg3", "odd"], chain_dir=str(tmp_path))
    assert status == "fixpoint"
    for stat, (L, U) in bounds.items():
        S, scale = statistic(stat, 3, 4, 8)
        assert L <= expr_value(wagner(), S) * scale <= U
    assert (tmp_path / "chain.jsonl").exists() and (tmp_path / "final_phase1.json.gz").exists()


def test_propagation_detects_an_empty_region():
    from flagcheck.ranges import propagate

    status, _ = propagate(3, 4, 9, 5, "all", (), None, False, ["red2"])
    assert status == "infeasible"
