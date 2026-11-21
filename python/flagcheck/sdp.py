"""Building plain flag algebra SDPs, solving them with CSDP, and rounding the result to an
exact certificate.

Convention (one formulation for everything). Let s range over density vectors on the
admissible n-vertex graphs: s >= 0, sum s = 1, M_sigma(s) = sum_G s_G Y_{G,sigma} PSD,
a_j . s = alpha_j (equalities), b_i . s <= beta_i (inequalities). We prove

    f . s >= d   for every such s

by exhibiting PSD A_sigma, free y and u >= 0 such that for every graph G

    r_G := f_G - sum_sigma <Y_{G,sigma}, A_sigma> - y . (a_G - alpha) - u . (beta - b_G) >= d.

Summing r_G * s_G over a feasible s gives f.s - (>= 0) - 0 - (>= 0) >= d. To prove
infeasibility, use f = 0: any certificate with d > 0 is a contradiction.

The proven d for a rational certificate is exactly min_G r_G, computed in `certificate.py`.
"""

import os
import shutil
import subprocess
import tempfile
from dataclasses import dataclass, field
from fractions import Fraction
from typing import Dict, List, Optional, Sequence, Tuple

from .algebra import Combo, Y_matrix, lift, sdp_blocks
from .graphs import Flag, canonical, generate

SparseSym = Dict[Tuple[int, int], Fraction]  # upper triangle (i <= j) of a symmetric matrix


@dataclass
class FlagSDP:
    n: int
    colors: int
    graphs: List[Flag]
    blocks: List[Tuple[Flag, List[Flag]]]
    Y: List[List[SparseSym]]               # Y[g][b]
    f: List[Fraction]
    eqs: List[Tuple[List[Fraction], Fraction]] = field(default_factory=list)
    ineqs: List[Tuple[List[Fraction], Fraction]] = field(default_factory=list)
    block_kinds: List[str] = field(default_factory=list)  # "plain" / "localizing"; empty = all plain

    def graph_index(self) -> Dict[Flag, int]:
        return {G: i for i, G in enumerate(self.graphs)}


def combo_to_vector(expr: Combo, n: int, colors: int, zeros, graphs: List[Flag]) -> List[Fraction]:
    """Chain-rule a combination of graphs (s = 0, at most n vertices) up to n vertices and
    return its coefficient vector over `graphs`."""
    if not isinstance(expr, dict):
        from .spec import evaluate  # expression tree (spec.Expr)

        expr = evaluate(expr, colors, zeros)
    index = {G: i for i, G in enumerate(graphs)}
    vec = [Fraction(0)] * len(graphs)
    for F in expr:
        if F.s != 0:
            raise ValueError("constraints and objectives must be unlabelled (average first)")
    for G, v in lift(expr, n, colors, zeros).items():
        vec[index[G]] += v
    return vec


def build_reference_sdp(
    n: int,
    colors: int,
    zeros: Sequence[Flag] = (),
    objective: Optional[Combo] = None,
    eqs: Sequence[Tuple[Combo, Fraction]] = (),
    ineqs: Sequence[Tuple[Combo, Fraction]] = (),
    include_empty_type: bool = False,
    localizing: Sequence[Tuple[int, object]] = (),
) -> FlagSDP:
    """Everything is computed by the Python reference from definitions.

    `localizing` lists (m, g): a block [[ g v v^T ]]_sigma >= 0 over all sigma-flags v on m
    vertices, where sigma is the type of g. It is only valid if g >= 0 pointwise at every
    root of every object considered; the caller is responsible for that (see the
    soundness tests)."""
    zeros = tuple(zeros)
    graphs = list(generate(n, colors, zeros))
    blocks = sdp_blocks(n, colors, zeros, include_empty_type)
    Y = []
    for G in graphs:
        row = []
        for sigma, flags in blocks:
            full = Y_matrix(G, sigma, flags)
            row.append({(i, j): full[i][j] for i in range(len(flags)) for j in range(i, len(flags)) if full[i][j]})
        Y.append(row)
    kinds = ["plain"] * len(blocks)
    for m, g in localizing:
        sigma, flags, entries = localizing_block(g, m, n, colors, zeros, graphs)
        blocks.append((sigma, flags))
        kinds.append("localizing")
        for gi in range(len(graphs)):
            Y[gi].append(entries[gi])
    f = combo_to_vector(objective or {}, n, colors, zeros, graphs)
    return FlagSDP(
        block_kinds=kinds,
        n=n,
        colors=colors,
        graphs=graphs,
        blocks=blocks,
        Y=Y,
        f=f,
        eqs=[(combo_to_vector(e, n, colors, zeros, graphs), Fraction(rhs)) for e, rhs in eqs],
        ineqs=[(combo_to_vector(e, n, colors, zeros, graphs), Fraction(rhs)) for e, rhs in ineqs],
    )


def localizing_block(g, m: int, n: int, colors: int, zeros, graphs: List[Flag]):
    """Entries lift([[ g F_i F_j ]], n) over all admissible sigma-flags F on m vertices.
    Returns (sigma, flags, per-graph sparse upper triangles)."""
    from .algebra import average, combo, multiply
    from .graphs import generate_flags
    from .spec import evaluate

    zeros = tuple(zeros)
    gc = evaluate(g, colors, zeros)
    if not gc:
        raise ValueError("localizing multiplier is zero")
    first = next(iter(gc))
    if any(F.s != first.s or F.type().cols != first.type().cols for F in gc):
        raise ValueError("localizing multiplier mixes labelled types")
    sigma = first.type()
    s = first.s
    if m < s or max(F.k for F in gc) + 2 * (m - s) > n:
        raise ValueError("[[g F F]] does not fit in n vertices")
    flags = list(generate_flags(sigma, m, colors, zeros))
    index = {G: i for i, G in enumerate(graphs)}
    entries = [{} for _ in graphs]
    for i, Fi in enumerate(flags):
        gFi = multiply(gc, combo((Fi, 1)), colors, zeros)
        for j in range(i, len(flags)):
            e = lift(average(multiply(gFi, combo((flags[j], 1)), colors, zeros)), n, colors, zeros)
            for G, v in e.items():
                entries[index[G]][(i, j)] = v
    return Flag(sigma.k, 0, sigma.cols), flags, entries


def build_from_problem(problem) -> FlagSDP:
    return build_reference_sdp(problem.n, problem.colors, problem.zeros, problem.objective, problem.eqs,
                               problem.ineqs, problem.include_empty_type, problem.localizing)


# ---------------------------------------------------------------------------
# CSDP (SDPA sparse format)
# ---------------------------------------------------------------------------
#
# CSDP solves   max tr(C X)  s.t.  tr(A_k X) = a_k,  X PSD (block diagonal).
# X = diag(A_sigma blocks..., LP block). LP block layout:
#   lambda_G (one per graph), d+, d-, y+_j, y-_j (per equality), u_i (per inequality),
#   and nu (normalisation slack, only when normalize=True).
# Row k (graph G):  <Y_G, A> + lambda_G + d+ - d- + (y+ - y-).(a_G - alpha) + u.(beta - b_G) = f_G
# Optional last row: sum tr(A_sigma) + sum (y+ + y-) + sum u + nu = 1.
# Objective: max d+ - d-.
# The CSDP dual vector restricted to the graph rows is the density vector s.
#
# Why the normalisation row: for an infeasibility proof (f = 0) the certificate cone is
# unbounded, so CSDP would report "dual infeasible" instead of a usable certificate. With
# the row, the optimum is finite. Its SDP dual is the phase-1 problem "min w such that
# M_sigma(s) + w I is PSD, |a.s - alpha| <= w, b.s - beta <= w, s in the simplex", so
# d* > 0 exactly when the original problem is infeasible. The row only guides the solver:
# check_certificate ignores it, because any certificate is valid whatever its scale.


@dataclass
class LPLayout:
    num_graphs: int
    num_eqs: int
    num_ineqs: int
    normalized: bool = False

    @property
    def nu(self):
        if not self.normalized:
            raise AttributeError("no normalisation slack")
        return self.num_graphs + 2 + 2 * self.num_eqs + self.num_ineqs

    @property
    def d_plus(self):
        return self.num_graphs

    @property
    def d_minus(self):
        return self.num_graphs + 1

    def y_plus(self, j):
        return self.num_graphs + 2 + 2 * j

    def y_minus(self, j):
        return self.num_graphs + 3 + 2 * j

    def u(self, i):
        return self.num_graphs + 2 + 2 * self.num_eqs + i

    @property
    def size(self):
        return self.num_graphs + 2 + 2 * self.num_eqs + self.num_ineqs + (1 if self.normalized else 0)


def write_sdpa(sdp: FlagSDP, path: str, normalize: bool = False) -> LPLayout:
    m = len(sdp.graphs)
    layout = LPLayout(m, len(sdp.eqs), len(sdp.ineqs), normalize)
    lp_block = len(sdp.blocks) + 1
    rows = m + (1 if normalize else 0)
    with open(path, "w") as fh:
        fh.write(f"{rows}\n{len(sdp.blocks) + 1}\n")
        fh.write(" ".join(str(len(flags)) for _, flags in sdp.blocks) + f" -{layout.size}\n")
        rhs = [repr(float(x)) for x in sdp.f] + (["1.0"] if normalize else [])
        fh.write(" ".join(rhs) + "\n")
        fh.write(f"0 {lp_block} {layout.d_plus + 1} {layout.d_plus + 1} 1.0\n")
        fh.write(f"0 {lp_block} {layout.d_minus + 1} {layout.d_minus + 1} -1.0\n")
        for g in range(m):
            row = g + 1
            for b, entries in enumerate(sdp.Y[g]):
                for (i, j), v in entries.items():
                    fh.write(f"{row} {b + 1} {i + 1} {j + 1} {float(v)!r}\n")
            lp = [(g, 1.0), (layout.d_plus, 1.0), (layout.d_minus, -1.0)]
            for j, (a, alpha) in enumerate(sdp.eqs):
                c = float(a[g] - alpha)
                if c:
                    lp += [(layout.y_plus(j), c), (layout.y_minus(j), -c)]
            for i, (bvec, beta) in enumerate(sdp.ineqs):
                c = float(beta - bvec[g])
                if c:
                    lp.append((layout.u(i), c))
            for idx, c in lp:
                fh.write(f"{row} {lp_block} {idx + 1} {idx + 1} {c!r}\n")
        if normalize:
            row = m + 1
            for b, (_, flags) in enumerate(sdp.blocks):
                for i in range(len(flags)):
                    fh.write(f"{row} {b + 1} {i + 1} {i + 1} 1.0\n")
            lp_idx = [layout.y_plus(j) for j in range(layout.num_eqs)]
            lp_idx += [layout.y_minus(j) for j in range(layout.num_eqs)]
            lp_idx += [layout.u(i) for i in range(layout.num_ineqs)] + [layout.nu]
            for idx in lp_idx:
                fh.write(f"{row} {lp_block} {idx + 1} {idx + 1} 1.0\n")
    return layout


@dataclass
class NumericSolution:
    status: int
    primal_objective: float
    dual_objective: float
    blocks: list            # list of numpy arrays A_sigma
    lp: list                # LP block diagonal values
    densities: list         # CSDP dual vector y = density vector s
    layout: LPLayout
    log: str


def solve_csdp(sdp: FlagSDP, workdir: Optional[str] = None, params: Optional[dict] = None,
               normalize: bool = False) -> NumericSolution:
    own = workdir is None
    workdir = workdir or tempfile.mkdtemp(prefix="flagsdp_")
    prob = os.path.join(workdir, "problem.dat-s")
    layout = write_sdpa(sdp, prob, normalize)
    try:
        return run_csdp(prob, [len(flags) for _, flags in sdp.blocks], layout, workdir, params)
    finally:
        if own:
            for name in os.listdir(workdir):
                os.remove(os.path.join(workdir, name))
            os.rmdir(workdir)


SDPA_PARAMS = {"maxIteration": "200", "epsilonStar": "1.0E-8", "lambdaStar": "1.0E2", "omegaStar": "2.0",
               "lowerBound": "-1.0E5", "upperBound": "1.0E5", "betaStar": "0.1", "betaBar": "0.2",
               "gammaStar": "0.9", "epsilonDash": "1.0E-8"}


def run_sdpa(prob: str, block_sizes: List[int], layout: LPLayout, workdir: str,
             params: Optional[dict] = None, threads: Optional[int] = None) -> NumericSolution:
    """Run SDPA (double precision) on the same SDPA file as run_csdp.

    SDPA's standard form is CSDP's with primal and dual swapped: SDPA's xVec is CSDP's y
    (the densities), SDPA's yMat is CSDP's X (the certificate matrices), and objValDual
    (max F_0 . Y) is CSDP's primal objective. Solutions are printed with 17 significant
    digits (SDPA's default is 4, too coarse for exact rounding)."""
    import re

    import numpy as np

    p = dict(SDPA_PARAMS, **(params or {}))
    keys = ["maxIteration", "epsilonStar", "lambdaStar", "omegaStar", "lowerBound", "upperBound",
            "betaStar", "betaBar", "gammaStar", "epsilonDash"]
    with open(os.path.join(workdir, "param.sdpa"), "w") as fh:
        for k in keys:
            fh.write(f"{p[k]}\t{k}\n")
        for k in ("xPrint", "XPrint", "YPrint", "infPrint"):
            fh.write(f"%+.17e\t{k}\n")
    out = os.path.join(workdir, "problem.out")
    cmd = ["sdpa", prob, out] + (["-numThreads", str(threads)] if threads else [])
    proc = subprocess.run(cmd, cwd=workdir, capture_output=True, text=True)
    log = proc.stdout + proc.stderr
    if not os.path.exists(out):
        raise RuntimeError(f"SDPA produced no output (exit {proc.returncode}):\n{log}")
    text = open(out).read()
    if "param  is ./param.sdpa" not in text and "param  is param.sdpa" not in text:
        raise RuntimeError("SDPA did not read our param.sdpa (precision would be too low)")

    def scalar(name):
        m = re.search(rf"^{name}\s*=\s*(\S+)", text, re.M)
        return m.group(1) if m else None

    phase = scalar("phase.value")
    pobj = float(scalar("objValDual") or "nan")        # CSDP primal
    dobj = float(scalar("objValPrimal") or "nan")      # CSDP dual
    num = r"[-+]?\d+\.\d+e[-+]\d+"
    xvec = text.split("xVec =", 1)[1].split("xMat =", 1)[0]
    y = [float(t) for t in re.findall(num, xvec)]
    ymat = text.split("yMat =", 1)[1].split("main loop time", 1)[0]
    vals = [float(t) for t in re.findall(num, ymat)]
    expected = sum(d * d for d in block_sizes) + layout.size
    if len(vals) != expected:
        raise RuntimeError(f"SDPA yMat has {len(vals)} numbers, expected {expected}")
    blocks, pos = [], 0
    for d in block_sizes:
        M = np.array(vals[pos:pos + d * d]).reshape(d, d)
        blocks.append((M + M.T) / 2)
        pos += d * d
    lp = vals[pos:pos + layout.size]
    status = 0 if phase == "pdOPT" else 1
    return NumericSolution(status, pobj, dobj, blocks, lp, y, layout, log + f"\nSDPA phase: {phase}\n")


def run_csdp(prob: str, block_sizes: List[int], layout: LPLayout, workdir: str,
             params: Optional[dict] = None) -> NumericSolution:
    """Run CSDP on an existing SDPA file (written by write_sdpa or by `flagsdp --sdpa`)."""
    import numpy as np

    sol = os.path.join(workdir, "problem.sol")
    if params:
        with open(os.path.join(workdir, "param.csdp"), "w") as fh:
            for k, v in params.items():
                fh.write(f"{k}={v}\n")
    # Output goes to workdir/csdp.log as it is produced (line-buffered when stdbuf exists),
    # so the progress of a long solve can be watched.
    cmd = (["stdbuf", "-oL"] if shutil.which("stdbuf") else []) + ["csdp", prob, sol]
    log_path = os.path.join(workdir, "csdp.log")
    with open(log_path, "w") as fh:
        proc = subprocess.run(cmd, cwd=workdir, stdout=fh, stderr=subprocess.STDOUT)
    with open(log_path) as fh:
        log = fh.read()
    pobj = dobj = float("nan")
    for line in log.splitlines():
        if line.startswith("Primal objective value:"):
            pobj = float(line.split(":")[1])
        elif line.startswith("Dual objective value:"):
            dobj = float(line.split(":")[1])
    if not os.path.exists(sol):
        raise RuntimeError(f"CSDP produced no solution (exit {proc.returncode}):\n{log}")
    with open(sol) as fh:
        y = [float(x) for x in fh.readline().split()]
        blocks = [np.zeros((d, d)) for d in block_sizes]
        lp = [0.0] * layout.size
        for line in fh:
            parts = line.split()
            if not parts:
                continue
            mat, blk, i, j, v = int(parts[0]), int(parts[1]), int(parts[2]) - 1, int(parts[3]) - 1, float(parts[4])
            if mat != 2:
                continue  # 1 = dual slack Z, 2 = primal X
            if blk <= len(blocks):
                blocks[blk - 1][i, j] = v
                blocks[blk - 1][j, i] = v
            else:
                lp[i] = v
    return NumericSolution(proc.returncode, pobj, dobj, blocks, lp, y, layout, log)
