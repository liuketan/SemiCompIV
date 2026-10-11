# SemiCompIV

Simulation code for *Instrumental variables estimation with semi-competing risks*.
It implements the counting-process IV recursions for the three cumulative
transition effects of an illness–death model (onset, death before onset, death
after onset), their influence-function standard errors, and the full Monte Carlo
study reported in the manuscript.

## Requirements

- R ≥ 4.1 with the packages `Rcpp` and `parallel` (base).
- A C++ toolchain that R can use to compile `cp_core.cpp`
  (Xcode Command Line Tools on macOS, Rtools on Windows). No Fortran or
  Armadillo is needed.
- Parallel runs use `parallel::mclapply`, which forks on macOS and Linux. On
  Windows set `NCORE=1`.

## Quick example

```bash
Rscript sample_code.R
```

Generates one data set of size 2,000 under Design A, fits the three curves,
and prints estimates, influence-function standard errors, pointwise 95%
intervals, and the denominator diagnostics that should be reported with an
application. It runs in a few seconds.

## Files

| File | Purpose |
|---|---|
| `cp_core.cpp` | Event-time recursions and influence functions (C++ via Rcpp); naive Aalen comparator |
| `cp_lib.R` | Data-generating designs A–E, nuisance fits, estimator wrapper, population-limit solver |
| `cp_sim.R` | Monte Carlo driver for one scenario and performance summaries |
| `run_single.R` | Runs the scenario named by the `SCENARIO` environment variable |
| `run_manifest.sh` | Runs all 16 main scenarios |
| `run_bootstrap.R` | Refitted subject-bootstrap check (100 data sets × 99 resamples) |
| `compute_poplimits.R` | Probability limits of the post-onset recursion when entry balance fails |
| `aggregate_results.R` | Bias, ESE, ASE, coverage and Monte Carlo standard errors for every scenario |
| `make_tables_figures.R` | LaTeX tables and the B12 figure, written to `output/` |

## Reproducing the simulation study

```bash
NSIM=2000 NCORE=7 ./run_manifest.sh   # 16 scenarios; about 30 min on 7 cores
Rscript run_bootstrap.R               # bootstrap check, single core
Rscript compute_poplimits.R
Rscript aggregate_results.R           # writes performance_summary.csv
Rscript make_tables_figures.R         # writes output/
```

Completed scenarios are skipped when `run_manifest.sh` is re-run. Each
replicate uses a fixed seed, so results are reproducible for a given `NSIM`.

| Manuscript item | Scenarios |
|---|---|
| Table of identification-valid designs | `S1e_rho50_n2000` (A), `S2_ufree_n2000` (B) |
| Table of entry-balance violations | `S3a`, `S3b` (C), `S4a`, `S4b` (D); limits from `compute_poplimits.R` |
| Figure of B12 curves | `S1e`, `S3a`, `S4a` |
| Bootstrap check | `run_bootstrap.R` |
| Supplementary scenarios | `S1a`–`S1d`, `S1f` (sample size, Corr(G,X) = 0.3), `S5a`, `S5b` (weak instrument), `S6a`–`S6c` (censoring) |

## Scope

This repository contains the simulation code only. The electronic health
record application uses patient-level data that cannot be shared.
