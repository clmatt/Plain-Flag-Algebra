# Results

Each JSON file holds the parameters of a Ramsey feasibility problem, the solver output, and an
exact certificate. Re-check any file from scratch with

    cd python && python -m flagcheck.verify ../results/<file>.json [--reference] [--spot-check K]

A file proves R(k,l) <= N only if the verifier reports a positive proven bound.

## Committed certificates

Only certificates that prove something are committed; failed runs are summarised in
docs/LOG.md. "Hypotheses" are printed by `verify` (Ramsey upper bounds from DS1.18 and
the edge tables in python/flagcheck/edge_tables.py).

| file | proves | margin | verified |
|------|--------|--------|----------|
| R46_N42_n7_local_x-codeg-edges-green.json | R(4,6) <= 42 | 1.623e-3 | yes |
| R46_N41_n7_local_x-codeg-edges-green.json | R(4,6) <= 41 | 6.998e-7 | yes |
| R45_N28_n7_all_x-codeg-edges-green.json | R(4,5) <= 28 | 4.317e-4 | yes |
