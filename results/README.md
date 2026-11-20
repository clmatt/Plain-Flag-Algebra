# Results

Each JSON file holds the parameters of a Ramsey feasibility problem, the solver output, and an
exact certificate. Re-check any file from scratch with

    cd python && python -m flagcheck.verify ../results/<file>.json [--reference] [--spot-check K]

A file proves R(k,l) <= N only if the verifier reports a positive proven bound.
