# Refactor and Ramsey-bound plan

This plan and `docs/LOG.md` are the working state of the revival of the thesis code
(C. Mattes, *Optimization Techniques in Extremal Graph Theory*, 2023, Ch. I–III).
Anyone (human or agent) resuming the work should read both files first.

## Ground rules

1. **Nothing counts as proved without an exact certificate.** A bound is claimed only
   when an exact rational certificate has been checked by the independent Python
   verifier (`python/flagcheck/certificate.py`). Solver status messages are never proof.
2. **Soundness tests are mandatory.** Every Ramsey setup is also run on the parameters
   just below the known value, where real Ramsey graphs exist. There it must be feasible,
   and the blow-up of a real graph must satisfy every constraint exactly.
3. **Two independent implementations.** The C++ core builds the SDPs. The Python
   reference recomputes the same objects from definitions (enumeration over step
   graphons) and cross-checks them on small instances.
4. **Tests gate progress.** A phase is complete only when the whole suite passes.

## Mathematical setup (summary)

- Objects: complete graphs with `c` edge colours. Ramsey problems use 3 colours:
  0 = "clone" (green), 1 = red, 2 = blue.
- Blow-up constraints (forbidden 3-vertex configurations) make colour 0 an
  equivalence relation with homogeneous colours between classes. Any finite graph
  with arbitrary positive part weights is then a valid limit object.
- An admissible graph on exactly N vertices gives the balanced blow-up, with
  `p(green K_k) = 1/N^(k-1)` for all k.
- Plain flag algebra on n vertices: for each labelled type σ with |σ| ≡ n (mod 2),
  `M_σ(s) = Σ_G s_G Y_{G,σ} ⪰ 0` for the density vector s of every limit object.
- Infeasibility certificate: PSD `A_σ`, multipliers `y` (equalities) and `u ≥ 0`
  (inequalities `b·s ≤ β`) such that for every admissible n-vertex G,
  `Σ_σ <Y_{G,σ}, A_σ> + y·a_G - u·b_G  <  y·α - u·β`.
  Summing against any feasible s gives a contradiction.

## Phases

### Phase 0: environment
- [ ] WSL2 Ubuntu with g++, cmake, nauty, GMP, CSDP, SDPA-GMP, Python 3 + pytest
- [x] `refactor` branch; Windows-safe filename for `R(4,6)<=39.cpp`; LF line endings

### Phase 1: tests and an independent reference (before touching the C++)
- [ ] Python reference: coloured graphs, brute-force canonical forms, induced densities
- [ ] Exact limit densities of weighted step graphons / blow-ups
- [ ] Reference `Y_{G,σ}`; identity check `Σ_G s_G Y_{G,σ} = E_θ[v_θ v_θ^T]` on step graphons
- [ ] Exact certificate verifier (PSD by construction: `A = Qᵀ D Q`, `D ≥ 0`)
- [ ] Known Ramsey graphs: C5, the (3,4)-graph on 8 vertices, the (3,5)-graph on 13, Paley(17)
- [ ] High-level tests against the original C++ (where it builds):
  - [ ] Goodman: min K3 + empty-K3 = 1/4 (n = 3)
  - [ ] Mantel: triangle-free implies edge density ≤ 1/2
  - [ ] P3 inducibility 3/4 (n = 3)
  - [ ] C5 inducibility approaches 1/26; P4 values from Appendix A
  - [ ] R(3,3) ≤ 6 (n = 3), R(3,4) ≤ 9 (n = 5), R(4,4) ≤ 18 (n = 6)
  - [ ] Never prove false bounds: R(3,3) ≤ 5, R(3,4) ≤ 8, R(4,4) ≤ 17 must NOT be proved
  - [ ] Building blocks (chain rule, multiplication, averaging) match the reference

### Phase 2: refactor (C++17 core + Python driver)
- [ ] Exact big-integer rationals (GMP) throughout; no fixed-width overflow
- [ ] Modules: coloured graph + canonical labelling (nauty), generation with forbidden
      subgraphs, flag enumeration per labelled type, `Y` construction (ordered and
      unordered), equations (add, multiply, average, chain rule), problem export
- [ ] Problem export to SDPA sparse format, scaled for floating-point solvers, plus
      the exact data the verifier needs
- [ ] Python driver: problem specs (YAML/py), solver runs (CSDP, SDPA-GMP),
      rounding to exact certificates
- [ ] The Phase 1 tests pass against the new core
- [ ] Original code moved unchanged to `legacy/`

### Phase 3: Ramsey bounds
- [ ] Robust phase-1 formulation: minimise the constraint violation δ; δ* > 0 together with
      a rounded dual gives an exact infeasibility certificate
- [ ] Reproduce R(3,3), R(3,4), R(4,4) with verified certificates
- [ ] R(4,5): the objective formulation (thesis: 29 on n = 7), then feasibility with exact-N
      constraints and degree constraints; push N down as far as possible
- [ ] R(4,6): the thesis degree-polynomial idea (with the `deg24` bug fixed)
- [ ] Every claimed bound has a certificate in `results/` and passes the verifier
