"""Re-check a saved Ramsey result from scratch.

    python -m flagcheck.verify results/R45_N29_n7_none.json [--reference] [--spot-check K]

Only the problem parameters and the certificate are read from the file. The SDP data is
rebuilt (by the C++ core, or with --reference by the independent Python reference),
the certificate is re-checked in exact arithmetic, and the proven bound is compared with
the stored one. --spot-check K recomputes, with the Python reference, the plain-block Y
matrices and the objective/constraint coefficients of K random graphs, and compares them
exactly with the C++ data. Use it at sizes where a full reference rebuild is too slow.
"""

import argparse
import json
import random
import sys
from fractions import Fraction

from .algebra import Y_matrix
from .certificate import Certificate, check_certificate
from .ramsey_run import build_problem


def load_certificate(data) -> Certificate:
    c = data["certificate"]
    blocks = [([[Fraction(x) for x in q] for q in b["rows"]], [Fraction(w) for w in b["weights"]]) for b in c["blocks"]]
    return Certificate(blocks, [Fraction(v) for v in c["y"]], [Fraction(v) for v in c["u"]])


def spot_check(sdp, k: int, seed: int = 0) -> int:
    rng = random.Random(seed)
    graphs = rng.sample(range(len(sdp.graphs)), min(k, len(sdp.graphs)))
    kinds = sdp.block_kinds or ["plain"] * len(sdp.blocks)
    for g in graphs:
        G = sdp.graphs[g]
        for b, (sigma, flags) in enumerate(sdp.blocks):
            if kinds[b] != "plain":
                continue
            full = Y_matrix(G, sigma, flags)
            expected = {(i, j): full[i][j] for i in range(len(flags)) for j in range(i, len(flags)) if full[i][j]}
            if sdp.Y[g][b] != expected:
                raise AssertionError(f"Y mismatch for graph {g}, block {b}")
    return len(graphs)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("result")
    ap.add_argument("--reference", action="store_true", help="rebuild everything with the Python reference")
    ap.add_argument("--spot-check", type=int, default=0)
    a = ap.parse_args(argv)
    with open(a.result) as fh:
        data = json.load(fh)
    p = data["params"]
    problem = build_problem(p["k"], p["l"], p["N"], p["n"], p["degrees"])
    if a.reference:
        from .sdp import build_from_problem

        sdp = build_from_problem(problem)
    else:
        from .cppbackend import build_cpp_sdp

        sdp = build_cpp_sdp(problem)
    if a.spot_check:
        print(f"spot-checked {spot_check(sdp, a.spot_check)} graphs against the reference: OK")
    res = check_certificate(sdp, load_certificate(data))
    stored = Fraction(data["proven_bound"])
    ok = res.bound == stored
    print(f"R({p['k']},{p['l']}) <= {p['N']} (n={p['n']}, degrees={p['degrees']}): proven bound {res.bound} "
          f"(~{float(res.bound):.3e}); stored bound matches: {ok}; PROVES INFEASIBILITY: {res.bound > 0}")
    return 0 if ok and res.bound > 0 else 1


if __name__ == "__main__":
    sys.exit(main())
