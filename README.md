# Plain Flag Algebra

Tools for the plain flag algebra method (Razborov; see C. Mattes, *Optimization
Techniques in Extremal Graph Theory*, 2023): bounds on subgraph densities via
semidefinite programming, and proofs of non-existence of finite graphs (e.g. Ramsey
bounds) via the blow-up encoding of Lidický and Pfender. Every result comes with an
**exact rational certificate** that is checked independently.

## Layout

| path | what |
|------|------|
| `cpp/` | Fast exact core: admissible graph/flag generation, flag algebra in GMP rationals, SDP data (thesis Algorithm 2) as integer counts. CLI: `flagsdp`. |
| `python/flagcheck/` | Independent reference implementation (brute force, from definitions), problem specs, CSDP driver, rounding, **exact certificate checker**. |
| `python/tests/` | High-level tests (Goodman, Mantel, inducibility, Ramsey bounds), exact equivalence of C++ and reference, and soundness tests on real Ramsey graphs. |
| `legacy/` | The original thesis code, unchanged, plus a harness to test it. Known issues are listed in `legacy/README.md`. |
| `docs/PLAN.md`, `docs/LOG.md` | Work plan and running log of findings. |

## Building (Ubuntu / WSL)

```bash
sudo apt install build-essential cmake ninja-build libgmp-dev nauty libnauty-dev coinor-csdp python3-venv
python3 -m venv ~/venvs/flag && ~/venvs/flag/bin/pip install pytest numpy scipy
cmake -S cpp -B ~/build/flagalg -G Ninja && cmake --build ~/build/flagalg
~/build/flagalg/unit_tests
bash legacy/harness/build.sh ~/build/legacy      # optional: legacy regression tests
cd python && ~/venvs/flag/bin/python -m pytest -q
```

## How a proof works

1. A problem (`flagcheck/problems.py`) states what we want to prove: `objective >= d`,
   or infeasibility of a set of constraints.
2. `flagsdp` (or the slow reference) builds the exact SDP data.
3. CSDP solves the SDP numerically.
4. The solution is rounded to rational PSD-by-construction matrices and multipliers.
5. `check_certificate` computes, in exact arithmetic, the bound the certificate proves.
   Rounding can only weaken that bound, never make it wrong.
