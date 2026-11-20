"""The large-problem path (C++ writes the SDPA file; the certificate is checked by
streaming the flagsdp file) must agree with the in-memory path on small problems."""

import os
import random
import tempfile
from fractions import Fraction

import pytest

from flagcheck import cppbackend
from flagcheck.certificate import check_certificate, check_flagsdp_file
from flagcheck.problems import goodman, ramsey_feasibility
from flagcheck.ramsey import RED, monochromatic_clique
from flagcheck.ramsey_run import build_problem
from flagcheck.sdp import LPLayout, write_sdpa
from flagcheck.algebra import combo
from flagcheck.graphs import from_edges

from test_certificate import random_cert

pytestmark = pytest.mark.skipif(not cppbackend.available(), reason="flagsdp not built")


def with_inequality():
    p = ramsey_feasibility(3, 4, 9, 5)
    red_edge = from_edges(2, [(0, 1, RED)])
    p.ineqs = [(combo((red_edge, 1)), Fraction(1, 2)), (combo((red_edge, -1)), Fraction(-1, 9))]
    p.name = "with-ineq"
    return p


PROBLEMS = [goodman(4), ramsey_feasibility(4, 4, 18, 5), build_problem(3, 4, 9, 5, "all"), with_inequality()]


def parse_sdpa(path):
    with open(path) as fh:
        m = int(fh.readline())
        nb = int(fh.readline())
        sizes = fh.readline().split()
        rhs = [float(x) for x in fh.readline().split()]
        entries = {}
        for line in fh:
            a, b, i, j, v = line.split()
            key = (int(a), int(b), int(i), int(j))
            assert key not in entries, key
            entries[key] = float(v)
    return m, nb, sizes, rhs, entries


@pytest.mark.parametrize("problem", PROBLEMS, ids=lambda p: p.name)
@pytest.mark.parametrize("normalize", [False, True])
def test_cpp_sdpa_equals_python_sdpa(problem, normalize):
    with tempfile.TemporaryDirectory() as d:
        out, dat = cppbackend.run_flagsdp(problem, d, sdpa=True, normalize=normalize)
        sdp = cppbackend.load(out)
        py = os.path.join(d, "python.dat-s")
        write_sdpa(sdp, py, normalize)
        a, b = parse_sdpa(dat), parse_sdpa(py)
    assert a[:3] == b[:3]
    assert a[3] == pytest.approx(b[3], rel=1e-15, abs=0)
    assert set(a[4]) == set(b[4])
    for k, v in a[4].items():
        assert v == pytest.approx(b[4][k], rel=1e-14, abs=1e-300), k


@pytest.mark.parametrize("problem", PROBLEMS, ids=lambda p: p.name)
@pytest.mark.parametrize("seed", range(3))
def test_streaming_checker_equals_in_memory(problem, seed):
    with tempfile.TemporaryDirectory() as d:
        out, _ = cppbackend.run_flagsdp(problem, d)
        sdp = cppbackend.load(out)
        cert = random_cert(random.Random(seed), sdp)
        streamed = check_flagsdp_file(out, cert)
        assert cppbackend.read_shape(out) == (len(sdp.graphs), [len(f) for _, f in sdp.blocks],
                                             len(sdp.eqs), len(sdp.ineqs))
    in_memory = check_certificate(sdp, cert)
    assert streamed.slacks == in_memory.slacks
    assert streamed.bound == in_memory.bound
