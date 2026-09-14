#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>
#include <exception>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

// Pseudo-inverse from a single SVD; the same decomposition supplies the
// rank-deficiency count, so the matrix is factorized only once. Mirrors
// arma::pinv's driver strategy: try the fast divide-and-conquer SVD (dgesdd)
// first and fall back to the robust standard driver (dgesvd) when it fails
// to converge -- dgesdd alone fails routinely on the rank-deficient designs
// that binary covariates produce. rank_deficient is incremented (not reset)
// so callers can accumulate a per-recursion count; messages are the caller's
// responsibility to keep per-step output from flooding the console.
arma::mat safe_pinv(const arma::mat& X, double tol, bool& success,
                    int& rank_deficient) {
  arma::mat pinv_X;
  try {
    if (X.n_rows == 0 || X.n_cols == 0 || X.has_nonfinite()) {
      success = false;
      return pinv_X;
    }

    arma::mat U, V;
    arma::vec s;
    bool ok = arma::svd_econ(U, s, V, X, "both", "dc");
    if (!ok) {
      ok = arma::svd_econ(U, s, V, X, "both", "std");
    }
    if (!ok) {
      success = false;
      return pinv_X;
    }

    // Same convention as arma::pinv(X, tol): singular values less than tol
    // are treated as zero.
    const arma::uvec keep = arma::find(s >= tol);
    if (keep.n_elem < std::min(X.n_rows, X.n_cols)) {
      rank_deficient++;
    }

    if (keep.n_elem == 0) {
      pinv_X.zeros(X.n_cols, X.n_rows);
    } else {
      pinv_X = V.cols(keep) * arma::diagmat(1.0 / s(keep)) * U.cols(keep).t();
    }
    success = true;
  } catch (std::exception& ex) {
    Rcpp::Rcerr << "Exception in pseudo-inverse calculation: " << ex.what() << "\n";
    success = false;
  } catch (...) {
    Rcpp::Rcerr << "Unknown error occurred in pseudo-inverse calculation.\n";
    success = false;
  }
  return pinv_X;
}

// Backward recursion shared by qbridge() and qbridge_IPCW().
// When IPCW_weights has zero size the plain (uncensored) recursion is run.
static List qbridge_core(const arma::vec& TT,
                         const arma::vec& Q,
                         const arma::vec& qtime,
                         const arma::mat& Z,
                         const arma::mat& W,
                         const arma::mat& IPCW_weights,
                         const arma::vec& weights,
                         double tol) {
  const arma::uword N = TT.n_elem;
  const arma::uword K_q = qtime.n_elem;
  const arma::uword dim_Z = Z.n_cols;
  const bool use_ipcw = (IPCW_weights.n_elem > 0);
  int n_pinv_fail = 0;
  int n_rank_deficient = 0;

  // Counting-process increments and at-risk indicators; neither depends on
  // the recursion, so both are computed in one vectorized pass.
  arma::mat dN(N, K_q), risk_q(N, K_q);
  for (arma::uword j = 0; j < K_q; j++) {
    dN.col(j) = -arma::conv_to<arma::vec>::from(Q == qtime[j]);
    risk_q.col(j) = arma::conv_to<arma::vec>::from((Q <= qtime[j]) % (TT > qtime[j]));
  }

  arma::mat B(K_q + 1, dim_Z, arma::fill::zeros);
  arma::mat dB(K_q, dim_Z, arma::fill::zeros);
  arma::mat B_int(N, K_q + 1, arma::fill::zeros);
  arma::mat dB_int(N, K_q, arma::fill::zeros);
  arma::vec B_numer;   // dim_W; sized by its assignment below
  arma::mat B_denom;   // dim_W x dim_Z
  arma::mat Zv(N, dim_Z);

  for (arma::uword jj = K_q; jj-- > 0;) {
    const arma::uword j = jj;

    // exp of the integrated bridge at t+, combined with all case weights
    arma::vec v = weights % arma::exp(B_int.col(j + 1));
    if (use_ipcw) {
      v %= IPCW_weights.col(j);
    }

    // B_numer = W' (v . dN_j),  B_denom = W' diag(v . risk_j) Z
    B_numer = W.t() * (v % dN.col(j));
    Zv = Z.each_col() % (v % risk_q.col(j));
    B_denom = W.t() * Zv;

    bool pinv_success = false;
    arma::mat pinv_B_denom = safe_pinv(B_denom, tol, pinv_success, n_rank_deficient);

    if (!pinv_success) {
      n_pinv_fail++;
      // Carry the recursion forward so earlier time steps still start from
      // the accumulated integral rather than from zero.
      B_int.col(j) = B_int.col(j + 1);
      B.row(j) = B.row(j + 1);
      continue;
    }

    arma::vec dB_j = pinv_B_denom * B_numer;
    dB.row(j) = dB_j.t();
    dB_int.col(j) = Z * dB_j;
    B_int.col(j) = B_int.col(j + 1) - dB_int.col(j);
    B.row(j) = B.row(j + 1) - dB_j.t();
  }

  // One summary line instead of a message per time step; the counts are also
  // returned so R callers can warn programmatically.
  if (n_pinv_fail > 0) {
    Rcpp::Rcerr << "qbridge: pseudo-inverse failed at " << n_pinv_fail << " of "
                << K_q << " time steps (increment skipped, integral carried forward).\n";
  }

  return List::create(
    Named("dN") = dN,
    Named("risk_q") = risk_q,
    Named("dB") = dB,
    Named("B") = B,
    Named("B_int") = B_int,
    Named("dB_int") = dB_int,
    Named("B_numer") = B_numer,
    Named("B_denom") = B_denom,
    Named("qtime") = qtime,
    Named("n_pinv_fail") = n_pinv_fail,
    Named("n_rank_deficient") = n_rank_deficient
  );
}

static void check_inputs(const arma::vec& TT, const arma::vec& Q,
                         const arma::mat& Z, const arma::mat& W,
                         const arma::vec& weights) {
  const arma::uword N = TT.n_elem;
  if (Q.n_elem != N) {
    stop("Error: Length of Q does not match length of TT.");
  }
  if (Z.n_rows != N) {
    stop("Error: Number of rows in Z does not match length of TT.");
  }
  if (W.n_rows != N) {
    stop("Error: Number of rows in W does not match length of TT.");
  }
  if (weights.n_elem != N) {
    stop("Error: Length of weights does not match length of TT.");
  }
}

// [[Rcpp::export]]
SEXP qbridge(const arma::vec& TT,
             const arma::vec& Q,
             const arma::vec& qtime,
             const arma::mat& Z,
             const arma::mat& W,
             const arma::vec& weights,
             double tol = 0.001) {
  check_inputs(TT, Q, Z, W, weights);
  arma::mat no_ipcw;
  return qbridge_core(TT, Q, qtime, Z, W, no_ipcw, weights, tol);
}

// [[Rcpp::export]]
SEXP qbridge_IPCW(const arma::vec& TT,
                  const arma::vec& Q,
                  const arma::vec& qtime,
                  const arma::mat& Z,
                  const arma::mat& W,
                  const arma::mat& IPCW_weights,
                  const arma::vec& weights,
                  double tol = 0.001) {
  check_inputs(TT, Q, Z, W, weights);
  if (IPCW_weights.n_rows != TT.n_elem || IPCW_weights.n_cols != qtime.n_elem) {
    stop("Error: IPCW_weights must be N x length(qtime).");
  }
  return qbridge_core(TT, Q, qtime, Z, W, IPCW_weights, weights, tol);
}
