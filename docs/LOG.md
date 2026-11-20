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
