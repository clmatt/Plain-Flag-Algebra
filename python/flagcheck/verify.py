"""Re-check a saved Ramsey result from scratch.

    python -m flagcheck.verify results/R45_N29_n7_none.json [--reference] [--spot-check K]

Only the problem parameters and the certificate are read from the file. Then:

1. The SDP data is rebuilt by the C++ core (or, with --reference, entirely by the
   independent Python reference, which makes every check below redundant).
2. Independent checks against the Python reference (default mode):
   * the C++ graph list equals the reference's own enumeration of admissible n-vertex
     graphs (completeness matters for soundness: a missing graph could hide a false proof);
   * every plain block is the complete, duplicate-free list of sigma-flags, and the
     plain-block types are one per isomorphism class;
   * for K random graphs (--spot-check, default 25): every plain-block Y matrix, every
     objective/constraint coefficient and every localising-block entry is recomputed from
     definitions and compared exactly.
   For a clone-eliminated result (params "blowup"), the graph list must be the green-free
   admissible graphs, and each sampled entry for Q is recomputed as sum_G M[G][Q] times the
   value at G, with the column M[.][Q] counted directly (flagcheck.blowup.column).
3. The certificate is re-checked in exact arithmetic, and the proven bound must equal
   the stored one and be positive.
"""

import argparse
import json
import random
import sys
from fractions import Fraction

from .algebra import Y_matrix
from .certificate import Certificate, check_certificate
from .graphs import Flag, canonical, generate, generate_flags
from .ramsey_run import build_problem
from .spec import Avg, Lit, Mul


def load_certificate(data) -> Certificate:
    c = data["certificate"]
    blocks = [([[Fraction(x) for x in q] for q in b["rows"]], [Fraction(w) for w in b["weights"]]) for b in c["blocks"]]
    return Certificate(blocks, [Fraction(v) for v in c["y"]], [Fraction(v) for v in c["u"]])


def check_structure(sdp, problem) -> None:
    """Graph list and plain flag lists against the independent reference."""
    ref_graphs = set(generate(problem.n, problem.colors, problem.zeros))
    if problem.blowup:
        from .blowup import is_green_free

        ref_graphs = {G for G in ref_graphs if is_green_free(G)}
    cpp_graphs = [canonical(G) for G in sdp.graphs]
    if len(set(cpp_graphs)) != len(cpp_graphs):
        raise AssertionError("duplicate graphs in the C++ list")
    if set(cpp_graphs) != ref_graphs:
        raise AssertionError(f"graph list differs from the reference ({len(cpp_graphs)} vs {len(ref_graphs)})")
    kinds = sdp.block_kinds or ["plain"] * len(sdp.blocks)
    seen_types = set()
    for (sigma, flags), kind in zip(sdp.blocks, kinds):
        canon = [canonical(F) for F in flags]
        if len(set(canon)) != len(canon):
            raise AssertionError("duplicate flags in a block")
        if any(F.s != sigma.k or F.type().cols != sigma.cols for F in flags):
            raise AssertionError("block flag of the wrong type")
        if set(canon) != set(generate_flags(sigma, flags[0].k, problem.colors, problem.zeros)):
            raise AssertionError("block is not the complete list of sigma-flags")
        if kind == "plain":
            t = canonical(Flag(sigma.k, 0, sigma.cols))
            if t in seen_types:
                raise AssertionError("two plain blocks for the same type class")
            seen_types.add(t)


def spot_check(sdp, problem, sample) -> int:
    """Recompute, for the graphs in `sample`, every plain Y matrix (Python reference) and
    every objective / constraint coefficient and localising entry (directly from the
    definitions, flagcheck.direct: no symbolic products or averages)."""
    from .direct import GraphEvaluator

    kinds = sdp.block_kinds or ["plain"] * len(sdp.blocks)
    exprs = [problem.objective] + [e for e, _ in problem.eqs + problem.ineqs]
    vectors = [sdp.f] + [v for v, _ in sdp.eqs + sdp.ineqs]
    if len(sdp.eqs) != len(problem.eqs) or len(sdp.ineqs) != len(problem.ineqs):
        raise AssertionError("number of constraints differs from the problem")
    for j, ((_, r), (_, r2)) in enumerate(zip(sdp.eqs + sdp.ineqs, problem.eqs + problem.ineqs)):
        if Fraction(r) != Fraction(r2):
            raise AssertionError(f"constraint {j}: right-hand side differs")
    loc_specs = list(problem.localizing)
    loc_blocks = [b for b, kind in enumerate(kinds) if kind == "localizing"]
    if len(loc_specs) != len(loc_blocks):
        raise AssertionError("number of localising blocks differs from the problem")
    def values(G):
        """Reference values at the graph G: plain Y dicts, vector coefficients, localising dicts."""
        ev = GraphEvaluator(G)
        Ys = {}
        for b, (sigma, flags) in enumerate(sdp.blocks):
            if kinds[b] == "plain":
                full = Y_matrix(G, sigma, flags)
                Ys[b] = {(i, j): full[i][j] for i in range(len(flags)) for j in range(i, len(flags)) if full[i][j]}
        coeffs = [ev.value(e) if e else Fraction(0) for e in exprs]
        for (m, gexpr), b in zip(loc_specs, loc_blocks):
            flags = sdp.blocks[b][1]
            Ys[b] = {}
            for i in range(len(flags)):
                for j in range(i, len(flags)):
                    v = ev.value(Avg(Mul((gexpr, Lit(flags[i]), Lit(flags[j])))))
                    if v:
                        Ys[b][(i, j)] = v
        return Ys, coeffs

    if problem.blowup:
        from .blowup import column

        all_graphs = list(generate(problem.n, problem.colors, problem.zeros))
        cache = {}
    for g in sample:
        if problem.blowup:      # values at Q = sum_G M[G][Q] * values at G
            Ys, coeffs = {b: {} for b in range(len(sdp.blocks))}, [Fraction(0)] * len(exprs)
            for h, m in column(canonical(sdp.graphs[g]), all_graphs, problem.blowup).items():
                if h not in cache:
                    cache[h] = values(all_graphs[h])
                hY, hc = cache[h]
                for b, d in hY.items():
                    for key, v in d.items():
                        Ys[b][key] = Ys[b].get(key, Fraction(0)) + m * v
                coeffs = [c + m * v for c, v in zip(coeffs, hc)]
            Ys = {b: {k: v for k, v in d.items() if v} for b, d in Ys.items()}
        else:
            Ys, coeffs = values(sdp.graphs[g])
        for b in range(len(sdp.blocks)):
            if sdp.Y[g][b] != Ys.get(b, {}):
                raise AssertionError(f"Y mismatch for graph {g}, block {b} ({kinds[b]})")
        for j, (want, vec) in enumerate(zip(coeffs, vectors)):
            if want != vec[g]:
                raise AssertionError(f"objective/constraint {j - 1} coefficient mismatch for graph {g}")
    return len(sample)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("result")
    ap.add_argument("--reference", action="store_true", help="rebuild everything with the Python reference")
    ap.add_argument("--spot-check", type=int, default=25)
    a = ap.parse_args(argv)
    try:
        return _verify(a)
    except (ValueError, AssertionError) as ex:
        print(f"VERIFICATION FAILED: {ex}")
        return 1


def _verify(a):
    with open(a.result) as fh:
        data = json.load(fh)
    p = data["params"]
    problem = build_problem(p["k"], p["l"], p["N"], p["n"], p["degrees"], p.get("empty_type", False),
                            p.get("extras", ()), p.get("loc_size"), p.get("blowup", False))
    if problem.assumptions:
        print("hypotheses (published / computed Ramsey data, see ramsey_run.py and edge_tables.py):")
        for a_ in problem.assumptions:
            print(f"  {a_}")
    if a.reference:
        from .sdp import build_from_problem

        sdp = build_from_problem(problem)
        if problem.blowup:
            from .blowup import reduce_sdp

            sdp = reduce_sdp(sdp, problem.blowup, problem.zeros)
        print("rebuilt entirely with the Python reference")
        res = check_certificate(sdp, load_certificate(data))
    else:
        import shutil
        import tempfile

        from .certificate import check_flagsdp_file
        from .cppbackend import binary_fingerprint, load, read_shape, run_flagsdp

        fp = binary_fingerprint()
        if data.get("flagsdp_sha256") and data["flagsdp_sha256"] != fp:
            print("note: current flagsdp binary differs from the one that produced the result")
        workdir = tempfile.mkdtemp(prefix="verify_")
        try:
            out, _ = run_flagsdp(problem, workdir)
            m = read_shape(out)[0]
            sample = random.Random(0).sample(range(m), min(a.spot_check, m))
            sdp = load(out, keep_graphs=set(sample))  # header + Y of the sampled graphs only
            check_structure(sdp, problem)
            print(f"graph list ({len(sdp.graphs)}) and all {len(sdp.blocks)} flag lists match the reference")
            if sample:
                print(f"spot-checked {spot_check(sdp, problem, sample)} graphs (Y, constraints, localising "
                      f"blocks) against the reference: OK")
            res = check_flagsdp_file(out, load_certificate(data))
        finally:
            shutil.rmtree(workdir, ignore_errors=True)
    stored = Fraction(data["proven_bound"])
    ok = res.bound == stored
    print(f"R({p['k']},{p['l']}) <= {p['N']} (n={p['n']}, degrees={p['degrees']}): proven bound {res.bound} "
          f"(~{float(res.bound):.3e}); stored bound matches: {ok}; PROVES INFEASIBILITY: {res.bound > 0}")
    return 0 if ok and res.bound > 0 else 1


if __name__ == "__main__":
    sys.exit(main())
