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
