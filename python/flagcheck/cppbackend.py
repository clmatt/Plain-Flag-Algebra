"""Run the C++ `flagsdp` tool and load its exact output as a FlagSDP."""

import os
import subprocess
import tempfile
from fractions import Fraction
from typing import Optional

from .graphs import Flag
from .sdp import FlagSDP
from .spec import write_spec

FLAGSDP = os.environ.get("FLAGSDP", os.path.expanduser("~/build/flagalg/flagsdp"))


def available() -> bool:
    return os.path.exists(FLAGSDP)


def _flag(tokens) -> Flag:
    k, s = int(tokens[0]), int(tokens[1])
    cols = tuple(int(x) for x in tokens[2:])
    if len(cols) != k * (k - 1) // 2:
        raise ValueError("bad flag line")
    return Flag(k, s, cols)


def load(path: str) -> FlagSDP:
    with open(path) as fh:
        it = iter(fh.read().splitlines())

    def expect(word):
        parts = next(it).split()
        if parts[0] != word:
            raise ValueError(f"expected {word}, got {parts}")
        return parts

    if next(it).split() != ["FLAGSDP", "1"]:
        raise ValueError("not a FLAGSDP v1 file")
    n = int(expect("n")[1])
    colors = int(expect("colors")[1])
    graphs = [_flag(next(it).split()) for _ in range(int(expect("graphs")[1]))]
    blocks, denominators = [], []
    for b in range(int(expect("blocks")[1])):
        parts = expect("block")
        assert int(parts[1]) == b
        nflags, den = int(parts[3]), int(parts[4])
        sigma = _flag(parts[5:])
        flags = [_flag(next(it).split()) for _ in range(nflags)]
        blocks.append((sigma, flags))
        denominators.append(den)
    f = [Fraction(0)] * len(graphs)
    for _ in range(int(expect("objective")[1])):
        g, v = next(it).split()
        f[int(g)] = Fraction(v)
    eqs, ineqs = [], []
    line = next(it).split()
    while line[0] == "constraint":
        kind, rhs, nnz = line[1], Fraction(line[2]), int(line[3])
        vec = [Fraction(0)] * len(graphs)
        for _ in range(nnz):
            g, v = next(it).split()
            vec[int(g)] = Fraction(v)
        (eqs if kind == "eq" else ineqs).append((vec, rhs))
        line = next(it).split()
    if line[0] != "Y":
        raise ValueError(f"expected Y, got {line}")
    Y = [[{} for _ in blocks] for _ in graphs]
    for _ in range(int(line[1])):
        g, b, i, j, c = (int(x) for x in next(it).split())
        Y[g][b][(i, j)] = Fraction(c, denominators[b])
    if next(it).strip() != "end":
        raise ValueError("missing end marker")
    return FlagSDP(n=n, colors=colors, graphs=graphs, blocks=blocks, Y=Y, f=f, eqs=eqs, ineqs=ineqs)


def build_cpp_sdp(problem, workdir: Optional[str] = None, threads: Optional[int] = None) -> FlagSDP:
    own = workdir is None
    workdir = workdir or tempfile.mkdtemp(prefix="flagsdp_cpp_")
    spec = os.path.join(workdir, "problem.spec")
    out = os.path.join(workdir, "problem.flagsdp")
    write_spec(problem, spec)
    cmd = [FLAGSDP, spec, out] + (["--threads", str(threads)] if threads else [])
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"flagsdp failed:\n{proc.stderr[-3000:]}")
    sdp = load(out)
    if own:
        for name in os.listdir(workdir):
            os.remove(os.path.join(workdir, name))
        os.rmdir(workdir)
    return sdp
