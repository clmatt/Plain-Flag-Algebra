"""End to end: a saved result must re-verify from scratch (C++ path and pure reference
path), and any tampering with the certificate must be rejected."""

import json
import os
import shutil
import tempfile

import pytest

from flagcheck import cppbackend
from flagcheck.ramsey_run import run
from flagcheck.verify import main as verify_main

pytestmark = pytest.mark.skipif(not cppbackend.available() or shutil.which("csdp") is None,
                                reason="needs flagsdp + CSDP")


@pytest.fixture(scope="module")
def result_file():
    d = tempfile.mkdtemp()
    run(3, 4, 9, 5, "all", results_dir=d)
    yield os.path.join(d, "R34_N9_n5_all.json")
    shutil.rmtree(d, ignore_errors=True)


@pytest.fixture(scope="module")
def extras_result_file():
    d = tempfile.mkdtemp()
    run(3, 4, 9, 5, "none", results_dir=d, extras=("green", "codeg", "edges", "lpsd"))
    yield os.path.join(d, "R34_N9_n5_none_x-codeg-edges-green-lpsd.json")
    shutil.rmtree(d, ignore_errors=True)


def test_verify_with_extras(extras_result_file, capsys):
    assert verify_main([extras_result_file, "--spot-check", "8"]) == 0
    out = capsys.readouterr().out
    assert "e(3,3," in out and "R(2,4) <= 4" in out       # hypotheses are printed
    assert verify_main([extras_result_file, "--reference"]) == 0


def test_verify_rejects_dropped_extras(extras_result_file, tmp_path):
    """The certificate uses the extra constraints, so it cannot verify without them."""
    data = json.load(open(extras_result_file))
    data["params"]["extras"] = []
    bad = tmp_path / "bad.json"
    bad.write_text(json.dumps(data))
    assert verify_main([str(bad), "--spot-check", "0"]) == 1


def test_verify_cpp_path(result_file):
    assert verify_main([result_file, "--spot-check", "10"]) == 0


def test_verify_reference_path(result_file):
    assert verify_main([result_file, "--reference"]) == 0


@pytest.mark.parametrize("field", ["weight", "multiplier", "params"])
def test_tampering_is_rejected(result_file, field, tmp_path):
    data = json.load(open(result_file))
    if field == "weight":
        blk = next(b for b in data["certificate"]["blocks"] if b["weights"])
        blk["weights"][0] = "1/1"
    elif field == "multiplier":
        data["certificate"]["y"][0] = "1000/1"
    else:
        data["params"]["N"] = 8          # claim a false bound R(3,4) <= 8 with the same certificate
    bad = tmp_path / "bad.json"
    bad.write_text(json.dumps(data))
    assert verify_main([str(bad), "--spot-check", "0"]) == 1


@pytest.fixture(scope="module")
def reduced_result_file():
    d = tempfile.mkdtemp()
    run(3, 4, 9, 5, "none", results_dir=d, extras=("green", "codeg", "edges", "nstats", "lpsd"), blowup=True)
    yield os.path.join(d, "R34_N9_n5_none_x-codeg-edges-green-lpsd-nstats_b.json")
    shutil.rmtree(d, ignore_errors=True)


def test_verify_reduced(reduced_result_file):
    assert json.load(open(reduced_result_file))["params"]["blowup"] is True
    assert verify_main([reduced_result_file, "--spot-check", "12"]) == 0
    assert verify_main([reduced_result_file, "--reference"]) == 0


@pytest.mark.parametrize("change", [("blowup", False), ("N", 8)])
def test_reduced_tampering_is_rejected(reduced_result_file, change, tmp_path):
    data = json.load(open(reduced_result_file))
    data["params"][change[0]] = change[1]
    bad = tmp_path / "bad.json"
    bad.write_text(json.dumps(data))
    assert verify_main([str(bad), "--spot-check", "0"]) == 1


@pytest.fixture(scope="module")
def levels_result_file():
    d = tempfile.mkdtemp()
    run(3, 4, 9, 5, "none", results_dir=d, extras=("green", "codeg", "edges", "nstats", "lpsd"), blowup="levels")
    yield os.path.join(d, "R34_N9_n5_none_x-codeg-edges-green-lpsd-nstats_bl.json")
    shutil.rmtree(d, ignore_errors=True)


def test_verify_levels(levels_result_file):
    assert json.load(open(levels_result_file))["params"]["blowup"] == "levels"
    assert verify_main([levels_result_file, "--spot-check", "12"]) == 0
    assert verify_main([levels_result_file, "--reference"]) == 0


@pytest.mark.parametrize("change", [("blowup", True), ("blowup", False), ("N", 8)])
def test_levels_tampering_is_rejected(levels_result_file, change, tmp_path):
    data = json.load(open(levels_result_file))
    data["params"][change[0]] = change[1]
    bad = tmp_path / "bad.json"
    bad.write_text(json.dumps(data))
    assert verify_main([str(bad), "--spot-check", "0"]) == 1
