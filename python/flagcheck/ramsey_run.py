"""Build, solve and certify Ramsey feasibility problems "is there a (k,l)-graph on exactly
N vertices?" with optional strengthening constraints.

Usage:
    python -m flagcheck.ramsey_run K L N n [--degrees none|avg|local|psd|all] [--threads T]

A result is saved to results/R<k><l>_N<N>_n<n>_<degrees>.json with the problem
parameters, the solver output and the exact certificate. `python -m flagcheck.verify`
re-checks it from scratch.

Every constraint added here is valid for the balanced blow-up of every (k,l)-graph on
exactly N vertices (see ramsey.py for the argument). The degree window uses only
proven *upper* bounds on smaller Ramsey numbers:
    red degree   <= R(k-1, l) - 1   (the red neighbourhood is a (k-1, l)-graph)
    blue degree  <= R(k, l-1) - 1   (the blue neighbourhood is a (k, l-1)-graph)
    => red degree in [N - R(k, l-1), R(k-1, l) - 1].
"""

import argparse
import json
import os
import time
from fractions import Fraction

from .algebra import combo
from .certificate import check_certificate, round_solution
from .graphs import generate_flags
from .problems import Problem
from .ramsey import (
    GREEN,
    NUM_COLORS,
    ROOT,
    green_clique_constraint,
    monochromatic_clique,
    one_flag,
    ramsey_zeros,
    red_flag,
)
from .spec import Add, Avg, Lit, Mul, Scale

# Proven upper bounds (exact values where known), S. Radziszowski, "Small Ramsey Numbers",
# Electron. J. Combin. DS1, revision #18 (April 24, 2026):
#   exact values R(3,3..9), R(4,4), R(4,5) = 25           -- Table Ia / items (b), (c), (g)
#   R(4,6) <= 40, R(4,7) <= 58, R(5,5) <= 46, R(5,6) <= 85 -- Table Ib (Angeltveit-McKay)
# An entry BELOW the true value would be unsound (the degree window would exclude real
# graphs); an entry above it only weakens the window. Trivial cases R(1,l) = 1, R(2,l) = l.
RAMSEY_UPPER = {
    (3, 3): 6, (3, 4): 9, (3, 5): 14, (3, 6): 18, (3, 7): 23, (3, 8): 28, (3, 9): 36,
    (4, 4): 18, (4, 5): 25, (4, 6): 40, (4, 7): 58, (5, 5): 46, (5, 6): 85,
}


def ramsey_upper(k: int, l: int) -> int:
    k, l = min(k, l), max(k, l)
    if k == 1:
        return 1
    if k == 2:
        return l
    return RAMSEY_UPPER[(k, l)]


def degree_window(k: int, l: int, N: int):
    return max(0, N - ramsey_upper(k, l - 1)), min(N - 1, ramsey_upper(k - 1, l) - 1)


def red_minus(c: Fraction):
    """The 1-flag expression R - c (R = red degree / N)."""
    return Add((Lit(red_flag()), Scale(-c, one_flag())))


def build_problem(k: int, l: int, N: int, n: int, degrees: str = "none", empty_type: bool = False) -> Problem:
    zeros = tuple(ramsey_zeros(k, l))
    eqs = [green_clique_constraint(j, N) for j in range(2, n + 1)]
    localizing = []
    lo, hi = degree_window(k, l, N)
    if lo > hi:
        raise ValueError(f"empty degree window [{lo}, {hi}]: already impossible")
    factors = [red_minus(Fraction(d, N)) for d in range(lo, hi + 1)]
    width = len(factors)
    if degrees in ("avg", "local", "all") and width + 1 > n:
        if degrees != "all":
            raise ValueError(f"degree polynomial of width {width} does not fit in n = {n}")
    if degrees in ("avg", "all") and width + 1 <= n:
        eqs.append((Avg(Mul(tuple(factors))), Fraction(0)))
    if degrees in ("local", "all") and width + 1 <= n:
        # [[P(R) F]] = 0 for every admissible 1-flag F with |F| + width <= n
        P = Mul(tuple(factors))
        for m in range(2, n - width + 1):
            for F in generate_flags(ROOT, m, NUM_COLORS, zeros):
                eqs.append((Avg(Mul((P, Lit(F)))), Fraction(0)))
    if degrees in ("psd", "all"):
        lower = red_minus(Fraction(lo, N))                                   # R - lo/N >= 0
        upper = Add((Scale(Fraction(hi, N), one_flag()), Scale(Fraction(-1), Lit(red_flag()))))  # hi/N - R >= 0
        m = (n - 2) // 2 + 1        # |g| = 2, so 2 + 2(m - 1) <= n
        localizing += [(m, lower), (m, upper)]
        m2 = (n - 3) // 2 + 1       # |g| = 3 for the product
        localizing.append((m2, Mul((lower, upper))))
    name = f"R{k}{l}_N{N}_n{n}_{degrees}" + ("_e" if empty_type else "")
    return Problem(name, n, NUM_COLORS, zeros, {}, eqs=eqs, normalize=True, localizing=localizing,
                   include_empty_type=empty_type)


def _frac(x) -> str:
    x = Fraction(x)
    return f"{x.numerator}/{x.denominator}"


def _git_revision() -> str:
    import subprocess

    try:
        here = os.path.dirname(os.path.abspath(__file__))
        rev = subprocess.run(["git", "-C", here, "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
        dirty = subprocess.run(["git", "-C", here, "status", "--porcelain"], capture_output=True, text=True).stdout
        return rev + ("-dirty" if dirty.strip() else "")
    except OSError:
        return "unknown"


def save_result(path, params, sol, cert, res, timings):
    from .cppbackend import binary_fingerprint

    data = {
        "params": params,
        "git_revision": _git_revision(),
        "flagsdp_sha256": binary_fingerprint(),
        "csdp_primal_objective": sol.primal_objective,
        "csdp_dual_objective": sol.dual_objective,
        "proven_bound": _frac(res.bound),
        "proven_bound_float": float(res.bound),
        "proved_infeasible": res.bound > 0,
        "timings": timings,
        "certificate": {
            "blocks": [{"rows": [[_frac(x) for x in q] for q in rows], "weights": [_frac(w) for w in ws]}
                       for rows, ws in cert.blocks],
            "y": [_frac(v) for v in cert.y],
            "u": [_frac(v) for v in cert.u],
        },
    }
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        json.dump(data, fh)


def run(k, l, N, n, degrees="none", threads=None, results_dir=None, bits=50, empty_type=False):
    from .cppbackend import build_cpp_sdp
    from .sdp import solve_csdp

    t = {}
    t0 = time.time()
    problem = build_problem(k, l, N, n, degrees, empty_type)
    sdp = build_cpp_sdp(problem, threads=threads)
    t["build"] = time.time() - t0
    t0 = time.time()
    sol = solve_csdp(sdp, normalize=True)
    t["solve"] = time.time() - t0
    t0 = time.time()
    cert = round_solution(sol, bits=bits)
    res = check_certificate(sdp, cert)
    t["check"] = time.time() - t0
    params = {"k": k, "l": l, "N": N, "n": n, "degrees": degrees, "empty_type": empty_type,
              "window": degree_window(k, l, N)}
    if results_dir:
        save_result(os.path.join(results_dir, f"{problem.name}.json"), params, sol, cert, res, t)
    return sdp, sol, res, t


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("k", type=int)
    ap.add_argument("l", type=int)
    ap.add_argument("N", type=int)
    ap.add_argument("n", type=int)
    ap.add_argument("--degrees", default="none", choices=["none", "avg", "local", "psd", "all"])
    ap.add_argument("--threads", type=int, default=None)
    ap.add_argument("--empty-type", action="store_true", help="include the empty type (even n only)")
    ap.add_argument("--results", default=os.path.join(os.path.dirname(__file__), "..", "..", "results"))
    a = ap.parse_args()
    sdp, sol, res, t = run(a.k, a.l, a.N, a.n, a.degrees, a.threads, a.results, empty_type=a.empty_type)
    print(f"R({a.k},{a.l})<={a.N} n={a.n} degrees={a.degrees} window={degree_window(a.k, a.l, a.N)}: "
          f"graphs={len(sdp.graphs)} blocks={len(sdp.blocks)} csdp={sol.primal_objective:.3e} "
          f"exact={float(res.bound):.3e} PROVED={res.bound > 0} "
          f"(build {t['build']:.0f}s, solve {t['solve']:.0f}s, check {t['check']:.0f}s)", flush=True)


if __name__ == "__main__":
    main()
