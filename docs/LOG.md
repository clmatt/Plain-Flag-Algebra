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

## 2026-11-20: scaling to n = 8; n = 7 independently confirmed

- Exact cross-check at research scale: C++ and the Python reference agree on every graph
  (1360), flag list, Y entry and constraint vector for R(4,5), N = 29, n = 7 (23-minute
  run). So the C++ is independently confirmed at the size the results come from.
- Degree constraints at n = 7, N = 28 (localising PSD, localised equalities, both): none
  proves N = 28 (margins -4e-9, -3e-8, -1.1e-8). So n = 7 is exhausted for N = 28 with
  these families.
- n = 8 probe for R(4,5): 10681 graphs, 287 type blocks, built in < 1 min with 126 MB peak.
  Python cannot hold ~10^7 Y entries as Fractions in 3.5 GB, so:
  * `flagsdp --sdpa F --normalize` writes the CSDP input directly. It is identical, entry
    for entry, to Python's write_sdpa (tests/test_large_pipeline.py).
  * `check_flagsdp_file` checks certificates by streaming the exact file. It gives slacks
    identical to check_certificate on random certificates.
  * `ramsey_run` now always uses this path.

## 2026-11-20: plain R(4,6) n = 7 sweep; neighbourhood constraints

Plain sweep (n = 7, degree options none/avg/psd/local/all): N = 42 and N = 41 are not
proved with any option (margins between -6e-10 and -3.4e-8). Every failing run is a tiny
negative, which is what a genuinely feasible relaxation looks like, so solver precision
is not the bottleneck. The relaxation needs facts that 7-vertex samples cannot see.
The sweep was stopped before N = 40 and 39 in favour of the following.

New constraint families (`flagcheck/neighbourhoods.py`, `--extras`). Each is valid
pointwise on the balanced blow-up of every (k,l)-graph on N vertices and is localised
([[g F]] >= 0 or = 0 for every rooted F; `lpsd` adds [[g v v^T]] >= 0):

- green: every clone class has measure exactly 1/N, at 1- and 2-vertex roots.
- codeg: at a red or blue pair of roots, the vertices with a fixed colour pattern to the
  pair form a (k-a, l-b)-graph, so there are at most R(k-a, l-b) - 1 of them.
- edges: the red edge counts inside the red and blue neighbourhoods of every vertex lie
  in [min, max] over (k-1,l)- and (k,l-1)-graphs of that order. Imposed through the exact
  Lagrange interpolant over the degree window (when it fits in n) and the convex-hull
  lines (always).

Edge data (`flagcheck/edge_tables.py`):
- computed from McKay's complete lists of (3,4)-, (3,5)-, (3,6)- and (4,4)-graphs, his
  min/max-edge (4,5) sets for 4-23 vertices, and all 352366 (4,5)-graphs on 24 vertices;
- graph counts agree with his page, and 12937 sampled graphs were re-checked to be
  genuine Ramsey graphs;
- every (3,k) minimum agrees with Goedgebeur-Radziszowski (2013) Table 3;
- orders <= 6 agree with exhaustive brute force.

The user approved these published computational results as hypotheses. `verify` prints
every hypothesis a certificate relies on.

`verify` spot checks now use `flagcheck/direct.py`. It evaluates every constraint
coefficient and localising entry in the sampled graph straight from the definitions
(disjoint random subsets), which makes it a third implementation. It agrees with the
reference on random expressions.

Tests:
- every multiplier is >= 0 (green: = 0) at every root of blow-ups of C5, Wagner, C13,
  Paley(17), and C5 as a (3,4)-graph (cubic interpolant);
- Paley(17) is tight for the red-pair codegree bound and the red-neighbourhood maximum,
  and a bound one tighter is violated;
- C++ matches the reference exactly with all extras;
- R(3,4) <= 9 is proved and <= 8 never;
- no option refutes C13 (N = 13) or Paley(17) (N = 17, n = 5 and 6).

Benchmark, R(4,6), N = 42, n = 7:

| degrees | extras            | margin    | proved |
|---------|-------------------|-----------|--------|
| none    | green,codeg,edges | -7.4e-9   | no     |
| local   | green,codeg,edges | +1.623e-3 | YES    |

So there is no (4,6)-graph on 42 vertices, i.e. R(4,6) <= 42, one better than the
recursive bound 43 that the plain method could not improve. (Classical counting: with
degree 17 everywhere, the red triangle count from red neighbourhoods, at most 42*42/3,
contradicts the count of red edges in blue neighbourhoods, which is at least 42*116.) The
degree information has to be pointwise: the identity uses sum_v d_v^2.
`python -m flagcheck.verify` re-checked the certificate from scratch (margin 1.623e-3).
Convention reminder: infeasible at N means R(k,l) <= N. N = 40 matches the best known
bound (R(4,6) <= 40, DS1.18); N = 39 would be new.

R(4,6) with localised degrees + green, codeg, edges (n = 7): N = 40 margin -7.1e-9 and
N = 39 margin -4.0e-9, not proved. So this level reaches R(4,6) <= 41 (Mackey's bound),
but not the best known 40.

Second level (`--extras pedges`): at a red or blue pair of roots, the common
neighbourhood S with a fixed colour pattern is a (k-a, l-b)-graph. Its red edge count
lies in the edge-table range for |S|. In densities, D = |S|/N is a 3-vertex 2-flag and
T = 2 e(S)/N^2 is a 4-vertex 2-flag. |S| ranges over 0..R(k-a,l-b)-1, so the constraints
are the lower/upper convex-hull lines T >= a D + b (each re-checked exactly against every
table point). Tests:
- on Paley(17), D and T at concrete pairs equal |S|/N and 2e(S)/N^2 computed from the
  graph;
- every multiplier is >= 0 at every pair;
- C++ matches the reference;
- no option refutes C13 or Paley(17).

Test-helper fix: pointwise checks now only use root placements that induce the
multiplier's type. Elsewhere every value is trivially 0, which could make a "tight
somewhere" assertion pass vacuously. The earlier tightness results (Paley(17) red-pair
codegree bound, red-neighbourhood maximum) still hold with the stricter helper.

## 2026-11-20: SDPA vs CSDP

`--solver sdpa` runs SDPA 7.3.20 (Ubuntu package) on the same input file. SDPA's
standard form is CSDP's with primal and dual swapped: xVec = y (densities),
yMat = X (certificate), objValDual = CSDP's primal objective. param.sdpa is written to
print 17 significant digits (the default %+8.3e would ruin exact rounding) and
tolerances of 1e-8. On small problems the two solvers give the same exact-checker
verdicts and objectives within 1e-6 (tests/test_ramsey_run.py).

Benchmark on the same file, same background load, 2 threads each:

| problem (n = 7)                                | SDPA  | CSDP  | margin SDPA / CSDP  |
|------------------------------------------------|-------|-------|---------------------|
| R(4,6) N=41, local + green,codeg,edges (LP 5280) | 567 s | 508 s | +6.979e-7 / +6.998e-7 |
| R(4,6) N=42, plain (LP 1441)                   | 237 s | 169 s | -3.8e-9 / -4.1e-9   |

No speed-up, so CSDP stays the default. SDPA remains as an independently implemented
second solver. Both certificates go through the same exact checker anyway.

## 2026-11-20: R(4,5) calibration with extras; N = 40 with pedges

R(4,5) = 25, n = 7, degrees=all + green,codeg,edges:

| N  | margin    | proved |
|----|-----------|--------|
| 28 | +4.317e-4 | yes: R(4,5) <= 28 (plain n = 7 reached 29) |
| 27 | -1.6e-8   | no     |
| 26 | -1.6e-8   | no     |
| 25 | -5.1e-9   | no     |
| 24 | -4.2e-9   | no, as it must be (352366 graphs exist): soundness check passed |

`verify` re-checked R46_N41 (R(4,6) <= 41): passes.
R(4,6) N = 40, local + green,codeg,edges,pedges: -5.4e-9, not proved.

## 2026-11-20: N = 40 with localising PSD blocks; n = 8 stopped

R(4,6) N = 40, n = 7, degrees=all:

| extras                           | blocks | margin   | build / solve |
|----------------------------------|--------|----------|---------------|
| green,codeg,edges,lpsd           | 93     | -6.4e-9  | 1480 s / 382 s |
| green,codeg,edges,pedges,lpsd    | 132    | -1.0e-8  | 991 s / 643 s  |

Neither proves anything. Every n = 7 variant tried so far is feasible at N = 40, so N = 39
was not run (per the stopping rule). `verify` re-checked R45_N28 (R(4,5) <= 28): passes.

The n = 8 run of plain R(4,5) at N = 28 was stopped after about 3 h 45 m, at the user's
request: it was using about 2 GB of WSL memory and Windows was down to 337 MB free.

## 2026-11-20: option 2 (triple codegrees, larger localisation) at n = 7

New family `tcodeg`: |S| <= R(k-a, l-b) - 1 at red/blue triangles of roots, keeping only
the bounds no pair implies. For R(4,6), N = 40 there are 7, e.g. the common blue
neighbourhood of a blue triangle is a (4,3)-graph, so it has at most 8 vertices.

| problem (n = 7)  | degrees | extras                        | loc_size | margin   | proved |
|------------------|---------|-------------------------------|----------|----------|--------|
| R(4,6) N = 40    | local   | green,codeg,tcodeg,edges      | 6        | -1.0e-8  | no     |
| R(4,6) N = 40    | local   | green,codeg,edges             | 7        | -4.9e-9  | no     |
| R(4,6) N = 40    | local   | green,codeg,tcodeg,edges      | 7        | -9.5e-9  | no     |
| R(4,5) N = 27    | all     | green,codeg,tcodeg,edges      | 6        | -1.0e-8  | no     |

Nothing moves, so n = 7 looks exhausted for these constraint families. Next: an
overnight n = 8 run of R(4,6) N = 40 (local + green,codeg,edges, loc_size 6), with CSDP's
log streamed to ~/runs/night_n8/csdp.log and a memory guard.

## 2026-11-20 (late): n = 8 with extras does not fit this machine

The size probe for R(4,6) N = 40, n = 8, local + green,codeg,edges (loc_size 6) never
finished building: the kernel OOM-killed `flagsdp` at 3.3 GB resident (WSL's whole cap)
while it was generating the problem. So the night run was not started. For comparison,
plain R(4,5) at n = 8 builds in 126 MB, so the cost is in the localised constraints lifted
to n = 8 (each is held as an exact vector over every 8-vertex graph). Options: stream
those constraints out of flagsdp instead of holding them, run plain/local n = 8 here
(fits), or a bigger machine.

## 2026-11-21: streaming flagsdp; n = 8 with extras now builds in 166 MB

The OOM was flagsdp's dense (constraints x graphs) table of GMP rationals for the SDPA
LP rows, plus holding every lifted constraint and the Y body in memory. flagsdp now lifts
each constraint, writes it to both files and drops it; Y goes through a temporary file.
Exact output is byte-identical to the old binary on 9 problems (every extras family, up to
R(4,6) N = 40 n = 7); the SDPA files hold the same entries in a different order; all
pipeline/solver tests pass.

Probe, R(4,6) N = 40, n = 8, local + green,codeg,edges, loc_size 6 (under ulimit -v 2.5 GB):
11667 graphs, 293 blocks, 961 eq + 2096 le constraints, 18.8M LP entries, peak RSS 166 MB,
38 min (almost all of it lifting constraints to 8 vertices). Files: 1.4 GB.
The solve was started 2026-11-21 08:17 (~/jobs/night_n8.log, CSDP log in ~/runs/night_n8).

## 2026-11-21: n = 8 solve timing; indexed lift

R(4,6) N = 40, n = 8, local + green,codeg,edges (m = 11667): CSDP started 08:56, iterations
0-2 done by 09:37, iteration 3 not reached by 09:57, so 20+ min per iteration (10-20 h for
a full solve). WSL memory during the solve: ~2.25 GB. The memory guard stopped it at 09:57
when a comparison test (running the old, dense-table flagsdp) was started next to it;
the solve itself was not the cause.

Indexed unrooted lift (theory.cpp): p(F, H) read from an index of all k-subsets of the
n-vertex graphs instead of scanning every n-vertex graph per F. Byte-identical flagsdp
output on 9 problems; unit test against density() (mutation-checked); 42 s -> 16 s for
R(4,6) N = 40 n = 7.

## 2026-11-21: nstats (neighbourhood structure from complete lists) at n = 7

| problem (n = 7) | degrees | extras                        | nstats rows | margin   | proved |
|-----------------|---------|-------------------------------|-------------|----------|--------|
| R(4,5) N = 27   | all     | green,codeg,edges,nstats      | 23          | -1.9e-8  | no     |
| R(4,6) N = 40   | local   | green,codeg,edges,nstats      | 96          | -1.2e-8  | no     |

For R(4,5) N = 27 both neighbourhoods are fully listed at every degree, so the 3-vertex
statistics are used at full strength (16 more rows need 8 vertices). No effect: the
relaxation's slack is not in 3-vertex neighbourhood statistics.

Green-free 8-vertex graphs (for eliminating the clone variables): R(4,5) 5588 of 10681,
R(4,6) 6343 of 11667.

## 2026-11-21: clone elimination (--blowup)

For the balanced blow-up W of a red/blue graph H on N >= n vertices, the n-vertex
densities satisfy y_G = sum_Q M[G][Q] z_Q exactly, where z_Q is the probability that n
distinct uniform vertices of H induce the green-free graph Q, and
M[G][Q] = n!/aut(G) (N)_r/N^n aut(Q_r)/r! d(Q_r, Q) (r = number of green classes of G,
Q_r its quotient). z lies in the simplex over the green-free admissible graphs, so the SDP
over z (vectors c -> c M, Y_Q = sum_G M[G][Q] Y_G) is a valid relaxation, with fewer
graph rows. It is also at least as strong: y is confined to the image M(simplex) instead
of the whole simplex.

Implementation: flagcheck/blowup.py (reference), flagsdp spec directive `blowup <N>`
(integer counts over den_b N^n n!, so the file format and exact checker are unchanged),
ramsey_run --blowup, verify (graph list = green-free admissible graphs; sampled entries
recomputed from columns of M counted without lift()).

Checks: y = M z exactly on random blow-ups (N = 5..8, n = 3..5) and Paley17; columns of M
sum to 1; nauty automorphism counts equal brute force; C++ reduced output equals the
reference reduction exactly (4 problems, incl. all extras and localising blocks); reduced
runs prove R(3,4) <= 9 and never <= 8, never refute C13 / Paley17 at n = 6; verify and
--reference pass, tampering (blowup flag, N) is rejected. Dropping the aut(Q_r) factor in
flagsdp is caught by all of these (it produced false "proofs" for C13 / Paley17).
Unreduced flagsdp output is byte-identical to before (18 comparisons).

Reduced runs at n = 7 (graphs = number of green-free graph rows):

| problem (n = 7) | degrees | extras                   | graphs | margin (exact) | proved |
|-----------------|---------|--------------------------|--------|----------------|--------|
| R(4,5) N = 27   | all     | green,codeg,edges        | 627    | -6.9e-9        | no     |
| R(4,6) N = 40   | local   | green,codeg,edges        | 678    | -3.6e-9        | no     |
| R(4,6) N = 40   | local   | green,codeg,edges,nstats | 678    | -6.9e-9        | no     |

The tighter relaxation does not close the n = 7 gap (margins are solver-level zeros).
The n = 8 reduced run (6343 rows instead of 11667) was started 2026-11-21 11:00.
