# truncProxy_weighting

This repository contains code for proximal weighting estimation under dependent
left truncation, including both the original research implementation and the
`truncProxy` R package.

Associated paper:

- Yuyao Wang, Andrew Ying, and Ronghui Xu, "Proximal Survival Analysis for
  Dependent Left Truncation"
- <https://arxiv.org/pdf/2512.21283>

## Repository structure

- `pkg/truncProxy/`: R package source for `truncProxy`
- `src/`: research implementation of the estimators
- `main_simu3.R`: main simulation script
- `results/`: saved outputs from experiments and simulations
- `dev/`: scripts for verifying that optimized implementations reproduce the
  original results (see `dev/README.md`)

The implementation used for the paper is preserved at the `v0.1.0-paper` git
tag; later commits contain performance improvements and bug fixes that
reproduce the same results (verified to numerical tolerance by the scripts in
`dev/`).

## R package

The `truncProxy` package is developed in `pkg/truncProxy/`.

For package documentation, installation instructions, and a simulation-based
example, see:

- [`pkg/truncProxy/README.md`](pkg/truncProxy/README.md)

Local installation from the repository root:

```r
devtools::install("pkg/truncProxy")
```

## Research code

The `src/` directory contains the implementation used for method development
and simulation studies. The package in `pkg/truncProxy/` is the user-facing
distribution of the current `PQB_estimator()` and `PQB_IPCW_estimator()`
implementations, and also provides the compiled bridge-process routines used
by the research code.

## Reproducing the simulation study (Table 1)

Table 1 of the paper uses the `simu3` setting of `main_simu3.R`. Dependencies:

```r
install.packages(c("Rcpp", "RcppArmadillo", "survival"))
```

Run from the repository root:

```bash
Rscript main_simu3.R
```

Settings, as in the paper: 500 simulation replicates, `B = 200` bootstrap draws
per replicate, seed 213; `N = 1000` and time-varying IPCW weights by default.
Each replicate re-seeds its own random-number stream, so the results are
identical for any number of worker cores (`n_cores`, set near the top of the
script). One run at `N = 1000` takes roughly 1.5 hours on a laptop (about a
third of that at `N = 500`); for a quick smoke test use
`SIMU_REPS=5 Rscript main_simu3.R`. The true value of the estimand is read from
`results/simu3_truth.rda` (the commented block in the script shows how it was
computed).

Table 1 combines four runs — two sample sizes (`SIMU_N`) by two IPCW weighting
schemes (`SIMU_TV=TRUE` time-varying weights, the main rows; `SIMU_TV=FALSE`
case weights, the "-cw" rows). `SIMU_OUT` saves each run's summary:

```bash
SIMU_N=500  SIMU_TV=TRUE  SIMU_OUT=results/simu3_table1_n500_tvTRUE   Rscript main_simu3.R
SIMU_N=500  SIMU_TV=FALSE SIMU_OUT=results/simu3_table1_n500_tvFALSE  Rscript main_simu3.R
SIMU_N=1000 SIMU_TV=TRUE  SIMU_OUT=results/simu3_table1_n1000_tvTRUE  Rscript main_simu3.R
SIMU_N=1000 SIMU_TV=FALSE SIMU_OUT=results/simu3_table1_n1000_tvFALSE Rscript main_simu3.R
Rscript dev/make_table1.R > results/table1.tex   # assembles the LaTeX table
```

Provenance: the exact numbers in the published Table 1 were produced with the
code at the `v0.1.0-paper` tag (`git checkout v0.1.0-paper`). Later commits fix
a bug in the data generator (censoring times were recycled across simulated
subjects), so re-running on `main` draws different simulated data for the same
seed. The estimators themselves are unchanged: `dev/check_equivalence.R`
verifies that the current implementation reproduces the tagged one on fixed data
to numerical tolerance. Table 1 of the revised manuscript was regenerated with
the current code (tag `v0.2.1-revision1`) from the four runs listed above; the
summaries are in `results/simu3_table1_n*_tv*.csv` and the assembled table in
`results/table1.tex`.
