# Legacy (thesis-era) code

This is the original implementation used for the thesis (C. Mattes, 2023). It is kept
**unchanged** for reference and for regression testing. New work should use `cpp/` and
`python/`.

`harness/` builds `flags.h` (truncated before the Mosek-dependent copositive functions)
against Ubuntu's nauty, and runs standard problems through the original
`plainFlagAlgebra`:

```bash
bash legacy/harness/build.sh ~/build/legacy
cd python && ~/venvs/flag/bin/python -m pytest -q tests/test_legacy_backend.py
```

## Known issues

Found during the 2026 review. See `docs/LOG.md` for details.

1. `mainRamsey.cpp` (R(4,6) ≤ 40 attempt): `Equation deg24(deg23Vec, ...)` should use
   `deg24Vec`. As written, the added constraint ⟦(d-22)(d-23)²⟧ = 0 is false for real
   (4,6)-graphs with a vertex of red degree 24. It is an **unsound** extra constraint.
2. `general.h` `gcd`: `int s = b;` truncates 64-bit values, so it returns a wrong gcd for
   operands above 2^31. `fractions.h` `Frac::reduce` has the same truncation
   (`int c = gcd(...)`). Fractions with large denominators (e.g. 1/N^k from exact-N
   constraints) can be silently wrong.
3. `plainFlagAlgebra`: constraint coefficients are written via `static_cast<int>(...)` and
   can overflow.
4. Each SDP row is scaled by the integer LCM of its denominators, and the objective by the
   LCM of the right-hand sides (up to N^(n-1)). The resulting programs are badly
   conditioned.
5. Infeasibility was read from the solver status of an unnormalised problem. Posed as a
   normalised phase-1 problem, the same SDP data proves R(4,5) ≤ 29 on 7 vertices with
   CSDP (numerically), versus the ≤ 31 on 8 vertices reported in the thesis.

For n ≤ 7 the generated SDP data agrees with the independent Python reference: every
objective value and feasibility margin matches, on all tested problems.
