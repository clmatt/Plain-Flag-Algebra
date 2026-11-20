"""High-level tests through the full pipeline: problem spec -> SDP (Python reference) ->
CSDP -> rounding -> exact certificate check.

Every test asserts both STRENGTH (the exact proven bound is close to the known value) and
SOUNDNESS (the exact proven bound never beats a value a real construction attains)."""

import shutil
from fractions import Fraction

import pytest

from flagcheck.certificate import check_certificate, round_solution
from flagcheck.problems import (
    C5,
    P3,
    P4,
    goodman,
    inducibility,
    mantel,
    ramsey_feasibility,
    ramsey_objective,
)
from flagcheck.sdp import build_reference_sdp, solve_csdp

pytestmark = pytest.mark.skipif(shutil.which("csdp") is None, reason="needs CSDP")

TOL = 1e-6


def prove(problem):
    sdp = build_reference_sdp(problem.n, problem.colors, problem.zeros, problem.objective,
                              problem.eqs, problem.ineqs, problem.include_empty_type)
    sol = solve_csdp(sdp, normalize=problem.normalize)
    cert = round_solution(sol)
    result = check_certificate(sdp, cert)
    return sol, result


def test_goodman():
    sol, res = prove(goodman(3))
    assert abs(sol.primal_objective - 0.25) < TOL
    assert Fraction(1, 4) - Fraction(1, 10 ** 6) < res.bound <= Fraction(1, 4)  # G(n,1/2) attains 1/4


def test_mantel():
    sol, res = prove(mantel(3))
    assert abs(-sol.primal_objective - 0.5) < TOL
    assert Fraction(-1, 2) - Fraction(1, 10 ** 6) < res.bound <= Fraction(-1, 2)  # K_{n/2,n/2}


def test_p3_inducibility():
    sol, res = prove(inducibility(P3, 3))
    assert abs(-sol.primal_objective - 0.75) < TOL
    assert Fraction(-3, 4) - Fraction(1, 10 ** 6) < res.bound <= Fraction(-3, 4)  # K_{n/2,n/2}


@pytest.mark.parametrize("n,thesis_value", [(4, 0.3750000043), (5, 0.2185840013)])
def test_p4_matches_thesis_table(n, thesis_value):
    # Thesis Appendix A, Table A.1 (ordered method). Construction lower bound ~0.2014.
    sol, res = prove(inducibility(P4, n))
    assert abs(-sol.primal_objective - thesis_value) < 1e-5
    assert -res.bound >= Fraction(2014, 10000)
    assert -res.bound < thesis_value + 1e-5


def test_c5_matches_thesis_table():
    # Table A.2, n = 5 ordered: 0.0408376585. True value 1/26 (iterated blow-up of C5).
    sol, res = prove(inducibility(C5, 5))
    assert abs(-sol.primal_objective - 0.0408376585) < 1e-5
    assert -res.bound >= Fraction(1, 26)


def test_r33_objective_formulation():
    # Thesis Ch. II reports ~0.1786 here. The reference gets exactly the true optimum 1/5
    # (C5's blow-up attains it), whatever subset of the blow-up zeros is used. The
    # discrepancy with the thesis is investigated against the legacy code (docs/LOG.md).
    sol, res = prove(ramsey_objective(3, 3, 3))
    assert abs(sol.primal_objective - 0.2) < 1e-6
    assert res.bound > Fraction(1, 6)
    assert res.bound <= Fraction(1, 5)        # C5's balanced blow-up has green density 1/5


def test_r33_le_6_by_infeasibility():
    sol, res = prove(ramsey_feasibility(3, 3, 6, 3))
    assert res.bound > 0                      # exact proof that R(3,3) <= 6


def test_r33_no_false_proof_at_5():
    sol, res = prove(ramsey_feasibility(3, 3, 5, 3))
    assert res.bound <= 0                     # C5 exists, so this must not be "proved"
    assert abs(sol.primal_objective) < 1e-6


@pytest.mark.slow
def test_r34_le_9_by_infeasibility():
    sol, res = prove(ramsey_feasibility(3, 4, 9, 5))
    assert res.bound > 0


@pytest.mark.slow
def test_r34_no_false_proof_at_8():
    sol, res = prove(ramsey_feasibility(3, 4, 8, 5))
    assert res.bound <= 0
