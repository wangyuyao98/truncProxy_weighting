# Development helper scripts

Scripts used to verify that optimized implementations reproduce the
original (`v0.1.0-paper`) results.

- `ref_capture.R` — run on a baseline checkout (e.g. the `v0.1.0-paper`
  tag) to save a reference dataset plus all estimator outputs:

      Rscript dev/ref_capture.R ref_outputs_v0.rda

- `check_equivalence.R` — run on the current working tree to recompute
  everything on the saved dataset and compare at tolerance 1e-8:

      Rscript dev/check_equivalence.R ref_outputs_v0.rda

Both must be run from the repository root.
