"""How strong is a relaxation that does not prove infeasibility? Bound a statistic's range,
and use the integrality of counts to cut the region down.

    python -m flagcheck.ranges 4 6 40 7 --degrees local --extras green,codeg,edges [--stats red3,blue3]
    python -m flagcheck.ranges 4 6 40 7 --degrees local --extras green,codeg,edges \
        --stats deg15,deg16,deg17,odd --propagate --chain results/chain_R46_N40

Ranges. When R(k,l) <= N is not proved, the feasible region of the SDP is non-empty and
the phase-1 margin is ~0 whatever is added, so the margin cannot tell a useful constraint
from a useless one. Instead, for a statistic S (a linear combination of densities) we
prove lo <= S <= hi on the whole region: lo is a certified lower bound for min S
(objective S, no normalisation row) and hi = -(certified lower bound for min -S). Both are
exact (checked certificates), so [lo, hi] is a rigorous outer bound of the region's range.

Integrality. For a (k,l)-graph H on N vertices, scale * S(blow-up of H) is a whole number
for every statistic here (see statistic()). So every H satisfies
ceil(lo * scale) <= scale * S <= floor(hi * scale); if that is empty, no H exists.
Propagation adds the rounded bounds as constraints and repeats until nothing changes:
each new bound is valid for every H given the constraints already added, so the result
(and a final infeasibility, if reached) is a chain of certificates, saved in order.
"""

import argparse
import gzip
import json
import math
import os
import shutil
import tempfile
import time
from fractions import Fraction

from .certificate import check_flagsdp_file, round_solution
from .cppbackend import read_shape, run_flagsdp
from .graphs import Flag
from .problems import combo
from .ramsey import BLUE, GREEN, RED
from .ramsey_run import build_problem, degree_window
from .sdp import LPLayout, run_csdp
from .spec import Avg, Scale

STATS = {
    "red2": combo((Flag(2, 0, (RED,)), 1)),
    "red3": combo((Flag(3, 0, (RED,) * 3), 1)),
    "blue3": combo((Flag(3, 0, (BLUE,) * 3), 1)),
    "blue4": combo((Flag(4, 0, (BLUE,) * 6), 1)),
    "green2": combo((Flag(2, 0, (GREEN,)), 1)),
}


def statistic(name, k, l, N):
    """(expression, count scale): scale * value is a whole number for every (k,l)-graph H
    on N vertices, the value being measured in the balanced blow-up of H.
      red2, red3, blue3, blue4: a green-free F on j vertices has blow-up density
        j! c(F) / N^j with c(F) = number of j-subsets of H inducing F; scale N^j / j!.
      green2: 1/N always; no scale.
      degD: fraction of vertices of red degree D, as the average over the root of the
        Lagrange polynomial in (red degree / N) equal to 1 at D and 0 at the other degrees
        of the window (exact because every degree lies in the window); scale N.
      odd: fraction of vertices of odd red degree; their number is even, so scale N / 2."""
    if name in STATS:
        S = STATS[name]
        j = max(F.k for F in S)
        return S, (None if name == "green2" else Fraction(N ** j, math.factorial(j)))
    from .neighbourhoods import _lagrange

    lo, hi = degree_window(k, l, N)
    window = range(lo, hi + 1)
    if name.startswith("deg"):
        d = int(name[3:])
        if d not in window:
            raise ValueError(f"degree {d} is outside the window {lo}..{hi}")
        return Avg(_lagrange({e: Fraction(int(e == d)) for e in window}, window, N)), Fraction(N)
    if name == "odd":
        return Avg(_lagrange({e: Fraction(e % 2) for e in window}, window, N)), Fraction(N, 2)
    raise ValueError(f"unknown statistic {name}")


def lower_bound(problem, objective, threads=None, keep=None, normalize=False):
    """(certified lower bound for min objective over the problem's feasible region, CSDP
    objective, CSDP exit code). With keep=path the certificate is saved there. With
    objective {} and normalize=True this is the phase-1 feasibility problem (bound > 0
    proves the region empty)."""
    problem.objective = objective
    problem.normalize = normalize
    workdir = tempfile.mkdtemp(prefix="range_")
    try:
        out, dat = run_flagsdp(problem, workdir, threads=threads, sdpa=True, normalize=normalize)
        m, sizes, neq, nle = read_shape(out)
        sol = run_csdp(dat, sizes, LPLayout(m, neq, nle, normalize), workdir)
        cert = round_solution(sol, bits=50)
        res = check_flagsdp_file(out, cert)
        if keep:
            from .ramsey_run import _frac

            with gzip.open(keep, "wt") as fh:
                json.dump({"bound": _frac(res.bound),
                           "certificate": {"blocks": [{"rows": [[_frac(x) for x in q] for q in rows],
                                                       "weights": [_frac(w) for w in ws]}
                                                      for rows, ws in cert.blocks],
                                           "y": [_frac(v) for v in cert.y], "u": [_frac(v) for v in cert.u]}},
                          fh)
        return res.bound, sol.primal_objective, sol.status
    finally:
        shutil.rmtree(workdir, ignore_errors=True)


def stat_range(k, l, N, n, degrees="none", extras=(), loc_size=None, blowup=False, stat="red3", threads=None,
               cuts=(), keep=None):
    """(lo, hi, (csdp lo, csdp hi)) for the statistic over the region with the extra cuts
    [(expr, rhs) meaning expr <= rhs]."""
    S, _ = statistic(stat, k, l, N)
    lo, lo_f, c1 = lower_bound(_with_cuts(k, l, N, n, degrees, extras, loc_size, blowup, cuts), S, threads,
                               keep and keep + "_min.json.gz")
    mhi, hi_f, c2 = lower_bound(_with_cuts(k, l, N, n, degrees, extras, loc_size, blowup, cuts),
                                Scale(Fraction(-1), S), threads, keep and keep + "_max.json.gz")
    return lo, -mhi, (lo_f, -hi_f, c1, c2)


def _with_cuts(k, l, N, n, degrees, extras, loc_size, blowup, cuts):
    p = build_problem(k, l, N, n, degrees, False, extras, loc_size, blowup)
    p.ineqs = list(p.ineqs) + list(cuts)
    return p


def feasibility(k, l, N, n, degrees, extras, loc_size, blowup, cuts, threads=None, keep=None):
    """Certified phase-1 bound with the cuts: > 0 proves that no (k,l)-graph on N exists."""
    return lower_bound(_with_cuts(k, l, N, n, degrees, extras, loc_size, blowup, cuts), {}, threads, keep,
                       normalize=True)[0]


def propagate(k, l, N, n, degrees, extras, loc_size, blowup, stats, threads=None, chain_dir=None, rounds=5,
              log=print):
    """Repeat: bound every statistic on the current region, round the bounds to whole
    counts, add the rounded bounds as constraints. Stops when no bound improves or a
    statistic has no whole count left (then no (k,l)-graph on N vertices exists, given the
    hypotheses of build_problem). Returns (status, bounds) with bounds[stat] = (L, U) counts."""
    if chain_dir:
        os.makedirs(chain_dir, exist_ok=True)
    bounds = {}
    cuts = []
    step = 0
    for rnd in range(rounds):
        improved = False
        for stat in stats:
            S, scale = statistic(stat, k, l, N)
            keep = chain_dir and os.path.join(chain_dir, f"step{step:03d}_{stat}")
            plain = [(e, r) for e, r, _ in cuts]
            try:
                lo, hi, num = stat_range(k, l, N, n, degrees, extras, loc_size, blowup, stat, threads, plain, keep)
            except RuntimeError:                # no CSDP solution at all
                lo, hi, num = Fraction(-10 ** 6), Fraction(10 ** 6), (None, None, 1, 1)
            if num[2] != 0 or num[3] != 0:      # CSDP trouble: perhaps the region is already empty
                b = feasibility(k, l, N, n, degrees, extras, loc_size, blowup, plain, threads,
                                chain_dir and os.path.join(chain_dir, f"step{step:03d}_phase1.json.gz"))
                log(json.dumps({"step": step, "phase1_bound": float(b)}))
                if b > 0:
                    return "infeasible", bounds
            L, U = math.ceil(lo * scale), math.floor(hi * scale)
            old = bounds.get(stat)
            rec = {"round": rnd, "step": step, "stat": stat, "lo": str(lo), "hi": str(hi), "count_lo": float(lo * scale),
                   "count_hi": float(hi * scale), "L": L, "U": U, "cuts_before": len(cuts)}
            log(json.dumps(rec))
            if chain_dir:
                with open(os.path.join(chain_dir, "chain.jsonl"), "a") as fh:
                    fh.write(json.dumps(rec) + "\n")
            step += 1
            if L > U:
                log(f"INFEASIBLE: {stat} has no whole count in [{float(lo * scale)}, {float(hi * scale)}]")
                return "infeasible", bounds
            if old is None or L > old[0] or U < old[1]:
                L, U = max(L, old[0]) if old else L, min(U, old[1]) if old else U
                bounds[stat] = (L, U)
                # S >= L / scale  and  S <= U / scale, i.e. -S <= -L/scale, S <= U/scale
                cuts = [c for c in cuts if c[2] != stat]
                cuts.append((Scale(Fraction(-1), S), -Fraction(L) / scale, stat))
                cuts.append((S, Fraction(U) / scale, stat))
                improved = True
        if not improved:
            break
    b = feasibility(k, l, N, n, degrees, extras, loc_size, blowup, [(e, r) for e, r, _ in cuts], threads,
                    chain_dir and os.path.join(chain_dir, "final_phase1.json.gz"))
    log(json.dumps({"final_phase1_bound": float(b), "cuts": len(cuts)}))
    return ("infeasible" if b > 0 else "fixpoint"), bounds


def main(argv=None):
    ap = argparse.ArgumentParser()
    for a in ("k", "l", "N", "n"):
        ap.add_argument(a, type=int)
    ap.add_argument("--degrees", default="none")
    ap.add_argument("--extras", default="")
    ap.add_argument("--loc-size", type=int, default=None)
    ap.add_argument("--blowup", action="store_true")
    ap.add_argument("--stats", default="red3,blue3")
    ap.add_argument("--threads", type=int, default=None)
    ap.add_argument("--out", default=None, help="append one JSON line per statistic")
    ap.add_argument("--propagate", action="store_true", help="add rounded bounds as constraints and repeat")
    ap.add_argument("--chain", default=None, help="directory for the propagation's certificates")
    ap.add_argument("--rounds", type=int, default=5)
    a = ap.parse_args(argv)
    extras = tuple(x for x in a.extras.split(",") if x)
    stats = a.stats.split(",")
    if a.propagate:
        status, bounds = propagate(a.k, a.l, a.N, a.n, a.degrees, extras, a.loc_size, a.blowup, stats, a.threads,
                                   a.chain, a.rounds)
        print(json.dumps({"status": status, "bounds": bounds}), flush=True)
        return
    for stat in stats:
        t0 = time.time()
        lo, hi, num = stat_range(a.k, a.l, a.N, a.n, a.degrees, extras, a.loc_size, a.blowup, stat, a.threads)
        lo_f, hi_f = num[:2]
        rec = {"k": a.k, "l": a.l, "N": a.N, "n": a.n, "degrees": a.degrees, "extras": sorted(extras),
               "loc_size": a.loc_size, "blowup": a.blowup, "stat": stat, "lo": float(lo), "hi": float(hi),
               "width": float(hi - lo), "csdp": [lo_f, hi_f], "seconds": round(time.time() - t0)}
        scale = statistic(stat, a.k, a.l, a.N)[1]
        if scale:   # the whole counts left inside the certified range
            rec.update(count_lo=float(lo * scale), count_hi=float(hi * scale),
                       whole_counts=max(0, math.floor(hi * scale) - math.ceil(lo * scale) + 1))
        print(json.dumps(rec), flush=True)
        if a.out:
            with open(a.out, "a") as fh:
                fh.write(json.dumps(rec) + "\n")


if __name__ == "__main__":
    main()
