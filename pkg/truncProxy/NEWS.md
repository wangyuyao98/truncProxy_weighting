# truncProxy 0.2.0

- Vectorized the compiled bridge-process recursion with Armadillo matrix
  operations; roughly 1.5x faster per call, with results identical up to
  floating-point summation order.
- Fixed the handling of pseudo-inverse failures inside the bridge
  recursion: the integrated bridge process is now carried forward across
  a failed time step instead of being reset to zero, and the number of
  failed steps is returned to R (`n_pinv_fail`), triggering a warning in
  `PQB_estimator()` and `PQB_IPCW_estimator()`.
- Vectorized `B_int_eval()` and clamped its interval lookup so evaluation
  times outside `[0, max(qtime)]` use the boundary column.
- Added validation of the `IPCW_weights` dimensions in `qbridge_IPCW()`
  and NA checks on the IPCW weights in `PQB_IPCW_estimator()`.
- Relaxed the regression-test tolerance from 1e-12 to 1e-8 to allow for
  BLAS summation-order differences.

# truncProxy 0.1.0

- Initial public package release.
- Added `PQB_estimator()` for proximal weighting estimation under dependent left truncation.
- Added `PQB_IPCW_estimator()` for IPCW-adjusted proximal weighting estimation under dependent left truncation and right censoring.
- Added compiled bridge-process routines based on `Rcpp` and
  `RcppArmadillo`.
- Added simulation-based examples, regression tests, and GitHub Actions
  `R CMD check` workflow.
