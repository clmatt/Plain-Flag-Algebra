"""Run the legacy (thesis-era) C++ through the test harness and solve its CBF output with
CSDP, so its results can be compared with the independent reference.

The legacy plainFlagAlgebra writes a CBF problem:
    max  c . x      s.t.  sum_j <F_rj, X_j> + sum_k a_rk x_k + b_r  <=  0   for each row r
with PSD variables X_j, and scalar variables x_k that are free (F) or non-negative (L+).
Every row r is scaled by an integer, and the objective by `mult` ("Divide answers by").
For maximisation problems the legacy proves 1 - f >= d, so the bound on f is 1 - d.
"""

import os
import re
import shutil
import subprocess
from collections import defaultdict
from dataclasses import dataclass
from typing import Dict, List, Optional

HARNESS = os.environ.get("LEGACY_HARNESS", os.path.expanduser("~/build/legacy/harness"))


@dataclass
class CBF:
    objsense: str
    psd_dims: List[int]
    var_cones: List[tuple]          # [(cone, count)], e.g. [("F", 3), ("L+", 2)]
    con_cones: List[tuple]          # [("L-", m)]
    obj_a: Dict[int, float]
    fcoord: Dict[tuple, float]      # (row, psdvar, i, j) with i >= j -> value
    acoord: Dict[tuple, float]      # (row, var) -> value
    bcoord: Dict[int, float]
    duplicates: int

    @property
    def num_rows(self):
        return sum(c for _, c in self.con_cones)

    @property
    def num_vars(self):
        return sum(c for _, c in self.var_cones)


def parse_cbf(path: str) -> CBF:
    with open(path) as fh:
        lines = [ln.strip() for ln in fh]
    i = 0
    sections = {}
    while i < len(lines):
        key = lines[i]
        i += 1
        if not key:
            continue
        body = []
        while i < len(lines) and lines[i]:
            body.append(lines[i])
            i += 1
        sections[key] = body
    dup = 0

    def sparse(name, nidx):
        nonlocal dup
        out = defaultdict(float)
        body = sections.get(name, [])
        if not body:
            return {}
        count = int(body[0])
        if count != len(body) - 1:
            raise ValueError(f"{name}: header says {count} entries, found {len(body) - 1}")
        for ln in body[1:]:
            parts = ln.split()
            key = tuple(int(p) for p in parts[:nidx])
            if key in out:
                dup += 1
            out[key] += float(parts[nidx])
        return dict(out)

    var_body = sections["VAR"]
    nvar, ncones = map(int, var_body[0].split())
    var_cones = [(ln.split()[0], int(ln.split()[1])) for ln in var_body[1:1 + ncones]]
    assert sum(c for _, c in var_cones) == nvar, "VAR cone sizes do not add up"
    con_body = sections["CON"]
    ncon, nccones = map(int, con_body[0].split())
    con_cones = [(ln.split()[0], int(ln.split()[1])) for ln in con_body[1:1 + nccones]]
    assert sum(c for _, c in con_cones) == ncon, "CON cone sizes do not add up"
    psd = sections.get("PSDVAR", ["0"])
    psd_dims = [int(x) for x in psd[1:1 + int(psd[0])]]
    obj = {k[0]: v for k, v in sparse("OBJACOORD", 1).items()}
    bco = {k[0]: v for k, v in sparse("BCOORD", 1).items()}
    return CBF(sections["OBJSENSE"][0], psd_dims, var_cones, con_cones, obj,
               sparse("FCOORD", 4), sparse("ACOORD", 2), bco, dup)


def cbf_to_sdpa(cbf: CBF, path: str, normalize: bool = False) -> None:
    """Write the CBF problem in CSDP's primal form (max C.X s.t. A_k . X = a_k)."""
    if cbf.objsense != "MAX":
        raise ValueError("only MAX problems are produced by the legacy code")
    if any(c != "L-" for c, _ in cbf.con_cones):
        raise ValueError("only L- constraints are produced by the legacy code")
    # LP block layout: scalar vars (free -> two entries), then one slack per row, then nu.
    lp_of_var = []  # list of [(lp_index, sign)]
    idx = 0
    for cone, count in cbf.var_cones:
        for _ in range(count):
            if cone == "F":
                lp_of_var.append([(idx, 1.0), (idx + 1, -1.0)])
                idx += 2
            elif cone == "L+":
                lp_of_var.append([(idx, 1.0)])
                idx += 1
            else:
                raise ValueError(f"unsupported variable cone {cone}")
    slack0 = idx
    m = cbf.num_rows
    lp_size = slack0 + m + (1 if normalize else 0)
    lp_block = len(cbf.psd_dims) + 1
    rows = m + (1 if normalize else 0)
    rhs = [-cbf.bcoord.get(r, 0.0) for r in range(m)] + ([1.0] if normalize else [])
    with open(path, "w") as fh:
        fh.write(f"{rows}\n{len(cbf.psd_dims) + 1}\n")
        fh.write(" ".join(str(d) for d in cbf.psd_dims) + f" -{lp_size}\n")
        fh.write(" ".join(repr(x) for x in rhs) + "\n")
        for k, c in cbf.obj_a.items():
            for li, sign in lp_of_var[k]:
                fh.write(f"0 {lp_block} {li + 1} {li + 1} {sign * c!r}\n")
        for (r, j, a, b), v in cbf.fcoord.items():
            i, jj = min(a, b), max(a, b)
            fh.write(f"{r + 1} {j + 1} {i + 1} {jj + 1} {v!r}\n")
        for (r, k), v in cbf.acoord.items():
            for li, sign in lp_of_var[k]:
                fh.write(f"{r + 1} {lp_block} {li + 1} {li + 1} {sign * v!r}\n")
        for r in range(m):
            fh.write(f"{r + 1} {lp_block} {slack0 + r + 1} {slack0 + r + 1} 1.0\n")
        if normalize:
            for j, d in enumerate(cbf.psd_dims):
                for i in range(d):
                    fh.write(f"{m + 1} {j + 1} {i + 1} {i + 1} 1.0\n")
            # Scalar variable 0 is the bound d itself. It must stay out of the budget,
            # otherwise d alone can use all of it and hide every certificate.
            bound_var_lp = {li for li, _ in lp_of_var[0]}
            for li in range(slack0):
                if li not in bound_var_lp:
                    fh.write(f"{m + 1} {lp_block} {li + 1} {li + 1} 1.0\n")
            fh.write(f"{m + 1} {lp_block} {lp_size} {lp_size} 1.0\n")


@dataclass
class LegacyResult:
    mult: int
    maximize: bool
    csdp_objective: float
    value: float          # d = objective / mult  (proves f >= d, or 1 - f >= d if maximize)
    bound: float          # the bound on f: d, or 1 - d when maximising
    cbf: CBF
    csdp_log: str
    harness_log: str


def run_legacy(args: List[str], workdir: str, normalize: bool = False) -> LegacyResult:
    os.makedirs(workdir, exist_ok=True)
    proc = subprocess.run([HARNESS] + [str(a) for a in args], cwd=workdir, capture_output=True, text=True)
    hlog = proc.stdout + proc.stderr
    if proc.returncode != 0:
        raise RuntimeError(f"legacy harness failed ({proc.returncode}):\n{hlog[-3000:]}")
    mult = int(re.search(r"Divide answers by (-?\d+)", hlog).group(1))
    maximize = re.search(r"HARNESS maximize=(\d)", hlog).group(1) == "1"
    cbf = parse_cbf(os.path.join(workdir, "sdp.cbf"))
    dat = os.path.join(workdir, "legacy.dat-s")
    sol = os.path.join(workdir, "legacy.sol")
    cbf_to_sdpa(cbf, dat, normalize)
    cp = subprocess.run(["csdp", dat, sol], cwd=workdir, capture_output=True, text=True)
    clog = cp.stdout + cp.stderr
    obj = float("nan")
    for line in clog.splitlines():
        if line.startswith("Primal objective value:"):
            obj = float(line.split(":")[1])
    value = obj / mult
    return LegacyResult(mult, maximize, obj, value, 1 - value if maximize else value, cbf, clog, hlog)


def harness_available() -> bool:
    return os.path.exists(HARNESS) and shutil.which("csdp") is not None
