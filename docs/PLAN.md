# Refactor and Ramsey-bound plan

This plan and `docs/LOG.md` are the working state of the revival of the thesis code
(C. Mattes, *Optimization Techniques in Extremal Graph Theory*, 2023, Ch. I–III).
Anyone (human or agent) resuming the work should read both files first.

## Ground rules

1. **Nothing counts as proved without an exact certificate** that passes
   `python -m flagcheck.verify`. Solver status messages are never proof.
2. **Soundness tests are mandatory.** Every constraint family is checked exactly on
   blow-ups of real Ramsey graphs (N = R(k,l) - 1), where it must hold.
3. **Two independent implementations.** C++ (`cpp/`) builds the SDPs; the Python
   reference (`python/flagcheck/`) recomputes them from definitions. They must agree
   exactly wherever the reference is fast enough.
4. **Tests gate progress.** Run `~/build/flagalg/unit_tests` and
   `cd python && ~/venvs/flag/bin/python -m pytest -q` before every commit.

## Mathematical setup (summary)

- Objects: complete graphs with `c` edge colours. Ramsey problems use 3 colours:
  0 = "clone" (green), 1 = red, 2 = blue.
- The blow-up zeros make green an equivalence relation with homogeneous colours between
  classes, so any finite graph with any positive class weights is a valid limit object.
- An admissible graph on exactly N vertices gives the balanced blow-up, with
  `p(green K_k) = 1/N^(k-1)`, and red 1-flag density deg_red(v)/N at a root in class v.
- Certificate convention: see the `python/flagcheck/sdp.py` docstring. Infeasibility =
  objective 0 and proven d > 0. The normalisation row makes the SDP bounded; its dual is
  the phase-1 "minimum uniform violation" problem.

## Phases

### Phase 0: environment [done]
- [x] WSL2 Ubuntu 26.04: g++ 15, CMake 4, nauty 2.9.3 (TLS), GMP, CSDP 6.2, SDPA
- [x] `refactor` branch (pushed); Windows-safe filenames; LF line endings

### Phase 1: tests and an independent reference [done]
- [x] Python reference (graphs, canonical forms, densities, chain rule, product,
      averaging, Y), verified by exact identities on step graphons, with mutation checks
- [x] Exact certificate checker (PSD by construction), checked against a naive version
- [x] Known Ramsey graphs and soundness tests of all constraint families
- [x] High-level tests: Goodman, Mantel, P3 = 3/4 (exact), P4 and C5 values from Appendix A,
      R(3,3) <= 6, R(3,4) <= 9 with no false proofs at 5 and 8
- [x] Legacy code harness: legacy agrees with the reference on all of the above

### Phase 2: refactor [done; review pending]
- [x] C++17 core: exact integer Y counts, GMP algebra, nauty canonical forms, OpenMP
- [x] `flagsdp` CLI + S-expression specs; Python loader
- [x] Exact equivalence with the reference (13+ problems, n <= 6; n = 7 job)
- [x] Localising PSD blocks
- [x] Legacy code moved to `legacy/` with its known issues documented
- [ ] Independent review of the soundness of the whole pipeline (reviewer agent)
- [ ] Unordered (thesis III.2) blocks, to reach larger n if needed

### Phase 3: Ramsey bounds [in progress]
- [x] Normalised phase-1 formulation; exact certificates at n = 7 (margins ~1e-7 survive rounding)
- [x] R(4,5) <= 29 on n = 7, exactly certified (objective form 1/d = 28.0856)
- [x] R(4,5) / R(4,6) with degree constraints only at n = 7: no improvement (N = 28 / N = 42, 41
      not provable; margins ~ -1e-8, i.e. the relaxation is feasible)
- [x] Neighbourhood constraints (green, codeg, edges; published edge data, user-approved),
      with soundness tests on real Ramsey graphs and a direct evaluator in verify
- [x] R(4,6) <= 42 (N = 42 infeasible) certified and verified, margin 1.6e-3
- [ ] R(4,5) calibration with extras: N = 28..25 should prove, N = 24 must not
- [ ] R(4,6): N = 41, 40 (= best known), 39 (new) with extras; then lpsd / n = 8
- [ ] R(4,5) n = 8 plain run (N = 28) result
- [ ] Every claimed bound has a certificate in `results/` that passes `flagcheck.verify`
