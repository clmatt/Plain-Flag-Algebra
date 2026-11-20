"""The legacy (thesis-era) C++ must agree with the independent reference on every
problem small enough for the reference. Build the harness first:
    bash legacy_harness/build.sh ~/build/legacy
"""

import os
import tempfile

import pytest

from flagcheck.certificate import check_certificate, round_solution
from flagcheck.legacy import harness_available, run_legacy
from flagcheck.problems import C5, P3, P4, goodman, inducibility, mantel, ramsey_feasibility, ramsey_objective
from flagcheck.sdp import build_reference_sdp, solve_csdp

pytestmark = pytest.mark.skipif(not harness_available(), reason="legacy harness not built")


def reference_value(problem):
    sdp = build_reference_sdp(problem.n, problem.colors, problem.zeros, problem.objective,
                              problem.eqs, problem.ineqs, problem.include_empty_type)
    return solve_csdp(sdp, normalize=problem.normalize).primal_objective


def legacy(args, normalize=False):
    with tempfile.TemporaryDirectory() as d:
        return run_legacy(args, d, normalize)


@pytest.mark.parametrize("args,problem,maximize", [
    (["goodman", 3], goodman(3), False),
    (["mantel", 3], mantel(3), True),
    (["ind_p3", 3], inducibility(P3, 3), True),
    (["ind_p4", 4], inducibility(P4, 4), True),
    (["ind_p4", 5], inducibility(P4, 5), True),
    (["ind_c5", 5], inducibility(C5, 5), True),
    (["ramsey_obj", 3, 3, 3], ramsey_objective(3, 3, 3), False),
    (["ramsey_obj", 4, 4, 5], ramsey_objective(4, 4, 5), False),
])
def test_objective_values_agree(args, problem, maximize):
    ref = reference_value(problem)
    r = legacy(args)
    # Reference proves objective >= d with objective = -G for maximisation, so the bound
    # on G is -d. The legacy maximisation bound is 1 - d_legacy.
    ref_bound = -ref if maximize else ref
    assert abs(r.bound - ref_bound) < 1e-5, (r.bound, ref_bound)


@pytest.mark.parametrize("k,l,N,n", [(3, 3, 6, 3), (3, 3, 5, 3), (3, 4, 9, 5), (3, 4, 8, 5)])
def test_feasibility_margins_agree(k, l, N, n):
    """Same phase-1 margin, up to solver accuracy (both normalise the same certificate
    cone, with the bound variable excluded from the budget)."""
    ref = reference_value(ramsey_feasibility(k, l, N, n))
    r = legacy(["ramsey_feas", k, l, N, n], normalize=True)
    margin = r.value - 1
    assert abs(margin - ref) < 1e-6 + 1e-3 * abs(ref), (margin, ref)
    assert (margin > 1e-7) == (ref > 1e-7)
