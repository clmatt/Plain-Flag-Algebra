# Progress log

Newest entries at the bottom. Record decisions, findings and anything surprising.

## 2026-11-20: initial review

Decisions (agreed with Connor): WSL2 Ubuntu; CSDP + SDPA-GMP; C++17 core + Python driver
and an independent verifier; commit to the `refactor` branch and push it to GitHub;
independent reviewer agents at phase boundaries.

Bugs found in the original code by reading it (not yet confirmed by tests):

1. `mainRamsey.cpp` (R(4,6) ≤ 40): `Equation deg24(deg23Vec, ...)` should use `deg24Vec`.
   As written, the constraint ⟦(d-22)(d-23)^2⟧ = 0 is false for real graphs, so this is an
   **unsound** over-constraint.
2. `general.h` `gcd`: `int s = b;` truncates 64-bit values, so it returns a wrong gcd once
   operands exceed 2^31.
3. `fractions.h` `Frac::reduce`: `int c = gcd(num, den);` has the same truncation.
4. `plainFlagAlgebra` ACOORD: `static_cast<int>(num * (constraintMult/den))` can overflow.
5. Scaling: each constraint row is multiplied by an integer LCM of denominators, and the
   objective by the LCM of the right-hand sides (up to 39^5). The resulting SDPs are very
   badly conditioned.

Evidence of a numerical or code problem in the thesis: the feasibility formulation
(Ch. II.2) reached only R(4,5) ≤ 31 on 8 vertices. The objective formulation, which has
a subset of its constraints, reached R(4,5) ≤ 29 on 7 vertices (Appendix A). Solved
correctly, the feasibility formulation can only be stronger.

Checked and believed correct: the ½ factor on off-diagonal entries in the CBF output; the
averaging normalisation `1/(C(n,j) C(j,s) s!)`; the sign conventions of the dual
multipliers; and the blow-up forbidden configurations.

## 2026-11-20: reference pipeline working

Environment: WSL2 Ubuntu 26.04, g++ 15.2, CMake 4.2, CSDP 6.2.0, SDPA, nauty 2.9.3. Python
venv at `~/venvs/flag` (numpy, scipy, pytest). Run the tests with
`cd python && ~/venvs/flag/bin/python -m pytest -q`.

The pipeline is: problem spec -> Python reference SDP -> CSDP -> rounding -> **exact**
certificate check (`flagcheck/certificate.py`). The `normalize=True` row makes
infeasibility problems bounded; its dual is the phase-1 "min max violation" problem.

Results (all certificates checked exactly):

| claim        | n | exact margin | status                                  |
|--------------|---|--------------|-----------------------------------------|
| R(3,3) <= 6  | 3 | 9.9e-3       | proved                                  |
| R(3,4) <= 9  | 5 | 1.6e-4       | proved                                  |
| R(3,5) <= 14 | 5 | 1.4e-5       | proved                                  |
| R(4,4) <= 18 | 6 | 1.1e-5       | proved                                  |
| R(4,4) <= 17 | 6 | -2.9e-10     | correctly not proved (Paley 17)         |
| R(4,5) <= 31 | 6 | 1.1e-6       | proved; N = 30..25 not proved           |

Thesis numbers reproduced: P4 (n = 4: 0.375, n = 5: 0.218584), C5 (n = 5: 0.0408377), and the
R(4,5) objective formulation at n = 6 (1/d = 30.25, so <= 31, as in Table A.6).

Discrepancy: the thesis gives min green K2 ≈ 0.1786 for R(3,3) on n = 3. The reference gets
exactly 0.2 (the true optimum, attained by C5's blow-up) with every subset of the blow-up
zeros. To investigate against the legacy C++.

Observation: the thesis-style averaged degree polynomial [[prod (R - d/N)]] = 0 is a
single linear equation and is weak, because the polynomial changes sign inside the window.
Pointwise validity gives [[P(R) F]] = 0 for every rooted F (`localized_degree_constraints`),
and the inequality form gives PSD localising blocks [[(R - lo/N) v v^T]] ⪰ 0 (not yet
implemented). At n = 6 neither the averaged nor the localised form improves R(4,5)
beyond 31. n = 6 is probably too small; we need n = 7-8, hence the C++ core.

## 2026-11-20: legacy code at n = 7; new C++ core

Legacy code (via harness + normalised phase-1 + CSDP), R(4,5), n = 7:
- objective form: 1/d = 28.086, so R(4,5) <= 29 (matches Table A.6);
- feasibility form: N = 31, 30, 29 numerically infeasible (margins 7.4e-7, 3.6e-7, 1.5e-7);
  N = 28: CSDP "stuck at edge of primal feasibility".

So the thesis-era code, posed as a normalised phase-1 problem, already reaches 29 on n = 7,
where the thesis (Mosek, unnormalised) reached only 31 on n = 8. The mathematics of the
legacy SDP generation is fine for n <= 7 (it agrees with the reference); the failure was
in how infeasibility was posed to the solver and read back. The margins are ~1e-7, so
exact certificates need care.

New C++ core (`cpp/`, CMake + Ninja, GMP, nauty with TLS, OpenMP):
- `Flag` (k <= 10) with canonical forms: brute force over non-root permutations when
  there are <= 4 non-roots (every flag in Y), nauty with root cells otherwise. Unit tests
  check 3000 random cases where brute force and nauty must induce the same classes.
- `Theory`: generation of admissible graphs and flags, plus density / product / average /
  lift in exact GMP rationals.
- `computeY`: Algorithm 2 with exact integer counts (no fractions, no overflow).
- `flagsdp` CLI: S-expression spec in, exact text problem out (`flagcheck/cppbackend.py`).
- `test_cpp_backend.py`: graphs, complete flag lists, every Y entry and every constraint
  vector equal the reference **exactly** on 13 problems (2 and 3 colours, empty type,
  Ramsey zeros, degree polynomials via C++ multiply/average). Mutation-checked.

## 2026-11-20: independent soundness review (reviewer agent) and fixes

The reviewer found no way for the mathematics or the exact checker to produce a false
proof. It checked: the certificate implication, the off-diagonal factor 2 against the
storage convention, the exact integer scaling, Y validity (one labelled representative
per type class), localising blocks, every Ramsey constraint family (blow-up zeros, green
powers, degree window with correct red/blue roles, degree polynomials, multipliers),
generation completeness, nauty canonical invariance with fixed roots, absence of
overflow, and the loader/format. Fixes for the gaps it raised:

1. `verify` was not independent by default. It now checks, against the reference: the
   full graph list (completeness; a missing graph could hide a false proof), every
   flag list, one plain block per type class, and for K sampled graphs every Y matrix,
   objective/constraint coefficient and localising-block entry.
2. C++ now enforces n <= kMaxVertices (10) instead of undefined behaviour.
3. Results record the git revision and the flagsdp sha256; verify notes mismatches.
4. RAMSEY_UPPER entries cited (DS1.18 Table Ia/Ib).
5. GMP rationals are canonicalised; degree options raise an error instead of silently
   skipping; an empty-type option was added; the loader rejects bad/duplicate Y keys.

Results so far (exact certificates): R(4,5) <= 29 on n = 7, in both the objective form
(1/d = 28.0856) and the feasibility form (N = 31, 30, 29; margins 7.4e-7, 3.6e-7, 1.4e-7).
N = 28 is not provable on n = 7 without extra constraints (margin -5e-9).
