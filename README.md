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
- `src/`: original research implementation of the estimators
- `main_simu3.R`: main simulation script
- `results/`: saved outputs from experiments and simulations

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

The `src/` directory contains the original implementation used for method
development and simulation studies. The package in `pkg/truncProxy/` is the
user-facing distribution of the current `PQB_estimator()` and
`PQB_IPCW_estimator()` implementations.

## Reproducing simulations

The main simulation workflow is in:

- `main_simu3.R`

Before running the simulation code, set the working directory to the repository
root and make sure the required R package dependencies are installed.
