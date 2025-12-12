#include <RcppArmadillo.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

// This is a example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//


// Helper: Robust pseudo-inverse function
arma::mat safe_pinv(const arma::mat& X, double tol, bool& success) {
  arma::mat pinv_X;
  try {
    if (X.n_rows == 0 || X.n_cols == 0) {
      success = false;
      Rcpp::Rcerr << "Matrix is empty, skipping pseudo-inverse.\n";
      return pinv_X;
    }
    
    // Check rank before attempting pseudo-inverse
    if (arma::rank(X, tol) < std::min(X.n_rows, X.n_cols)) {
      Rcpp::Rcerr << "Warning: Matrix is rank-deficient, pseudo-inverse may be unstable.\n";
    }
    
    pinv_X = arma::pinv(X, tol);
    success = true;
  } catch (std::exception &ex) {
    Rcpp::Rcerr << "Exception in pseudo-inverse calculation: " << ex.what() << "\n";
    success = false;
  } catch (...) {
    Rcpp::Rcerr << "Unknown error occurred in pseudo-inverse calculation.\n";
    success = false;
  }
  return pinv_X;
}



// [[Rcpp::export]]
SEXP qbridge(NumericVector TT,
             NumericVector Q,
             NumericVector qtime,
             NumericMatrix Z,
             NumericMatrix W,
             NumericVector weights,
             double tol = 0.001) {
  
  int N = TT.size();
  
  // Check that all inputs have consistent row counts
  if (Q.size() != N) {
    stop("Error: Length of Q does not match length of TT.");
  }
  if (Z.nrow() != N) {
    stop("Error: Number of rows in Z does not match length of TT.");
  }
  if (W.nrow() != N) {
    stop("Error: Number of rows in W does not match length of TT.");
  }
  if (weights.size() != N) {
    stop("Error: Length of weights does not match length of TT.");
  }
  
  int K_q = qtime.size();
  int dim_Z = Z.ncol();
  int dim_W = W.ncol();
  
  // Initialize all matrices and vectors
  NumericVector B_numer(dim_W);
  NumericMatrix B(K_q+1, dim_Z), dB(K_q, dim_Z), B_int(N, K_q+1), dB_int(N, K_q);  // Changed B_int and B to K_q+1 and with the following index changed to j+1
  NumericMatrix dN(N, K_q), risk_q(N, K_q);
  NumericMatrix B_denom(dim_W, dim_Z);
  
  // Initialize to zero
  std::fill(B_numer.begin(), B_numer.end(), 0.0);
  std::fill(B.begin(), B.end(), 0.0);
  std::fill(dB.begin(), dB.end(), 0.0);
  std::fill(B_int.begin(), B_int.end(), 0.0);
  std::fill(dB_int.begin(), dB_int.end(), 0.0);
  std::fill(dN.begin(), dN.end(), 0.0);
  std::fill(risk_q.begin(), risk_q.end(), 0.0);
  std::fill(B_denom.begin(), B_denom.end(), 0.0);
  
  // Main backward recursion loop
  for (int j = K_q - 1; j >= 0; j--) {
    
    // Reset for each time point
    std::fill(B_numer.begin(), B_numer.end(), 0.0);
    std::fill(B_denom.begin(), B_denom.end(), 0.0);
    
    // Compute increments and risks
    for (int i = 0; i < N; i++) {
      dN(i, j) = (Q[i] == qtime[j]) ? -1 : 0;
      risk_q(i, j) = (Q[i] <= qtime[j] && qtime[j] < TT[i]) ? 1 : 0;
      
      double exp_Bint = std::exp(B_int(i, j+1)); // this is at t+
      
      for (int dim_w = 0; dim_w < dim_W; dim_w++) {
        B_numer[dim_w] += weights[i] * W(i, dim_w) * exp_Bint * dN(i, j);
        for (int dim_z = 0; dim_z < dim_Z; dim_z++) {
          B_denom(dim_w, dim_z) += weights[i] * W(i, dim_w) * Z(i, dim_z) * risk_q(i, j) * exp_Bint;
        }
      }
    }
    
    // Convert B_denom to arma matrix for pseudo-inverse
    arma::mat B_denom_arma(B_denom.begin(), dim_W, dim_Z, false);
    
    // Calculate pseudo-inverse safely
    bool pinv_success = false;
    arma::mat pinv_B_denom = safe_pinv(B_denom_arma, tol, pinv_success);
    
    if (!pinv_success) {
      Rcpp::Rcerr << "Pseudo-inverse failed at time index j = " << j << "\n";
      continue; // Skip this time step if pseudo-inverse failed
    }
    
    // Compute dB = pinv(B_denom) * B_numer
    arma::colvec B_numer_arma(B_numer.begin(), dim_W, false);
    arma::colvec dB_arma = pinv_B_denom * B_numer_arma;
    
    // Copy back dB result to Rcpp matrix
    for (int dim_z = 0; dim_z < dim_Z; dim_z++) {
      dB(j, dim_z) = dB_arma(dim_z);
    }
    
    // Update B_int and B
    for (int i = 0; i < N; i++) {
      double dB_int_sum = 0.0;
      for (int dim_z = 0; dim_z < dim_Z; dim_z++) {
        dB_int_sum += Z(i, dim_z) * dB(j, dim_z);
      }
      dB_int(i, j) = dB_int_sum;
      B_int(i, j) = B_int(i, j + 1) - dB_int(i, j);
      // B_int(i, j) = (j + 1 < K_q) ? (B_int(i, j + 1) - dB_int(i, j)) : -dB_int(i, j);
    }
    
    for (int dim_z = 0; dim_z < dim_Z; dim_z++) {
      B(j, dim_z) = B(j + 1, dim_z) - dB(j, dim_z);
      // if (j == K_q - 1) {
      //   B(j, dim_z) = -dB(j, dim_z);
      // } else {
      //   B(j, dim_z) = B(j + 1, dim_z) - dB(j, dim_z);
      // }
    }
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
    Named("qtime") = qtime
  );
}


// [[Rcpp::export]]
SEXP qbridge_IPCW(NumericVector TT,
                  NumericVector Q,
                  NumericVector qtime,
                  NumericMatrix Z,
                  NumericMatrix W,
                  NumericMatrix IPCW_weights,
                  NumericVector weights,
                  double tol = 0.001) {
  
  int N = TT.size();
  
  // Check that all inputs have consistent row counts
  if (Q.size() != N) {
    stop("Error: Length of Q does not match length of TT.");
  }
  if (Z.nrow() != N) {
    stop("Error: Number of rows in Z does not match length of TT.");
  }
  if (W.nrow() != N) {
    stop("Error: Number of rows in W does not match length of TT.");
  }
  if (weights.size() != N) {
    stop("Error: Length of weights does not match length of TT.");
  }
  
  int K_q = qtime.size();
  int dim_Z = Z.ncol();
  int dim_W = W.ncol();
  
  // Initialize all matrices and vectors
  NumericVector B_numer(dim_W);
  NumericMatrix B(K_q+1, dim_Z), dB(K_q, dim_Z), B_int(N, K_q+1), dB_int(N, K_q);  // Changed B_int and B to K_q+1 and with the following index changed to j+1
  NumericMatrix dN(N, K_q), risk_q(N, K_q);
  NumericMatrix B_denom(dim_W, dim_Z);
  
  // Initialize to zero
  std::fill(B_numer.begin(), B_numer.end(), 0.0);
  std::fill(B.begin(), B.end(), 0.0);
  std::fill(dB.begin(), dB.end(), 0.0);
  std::fill(B_int.begin(), B_int.end(), 0.0);
  std::fill(dB_int.begin(), dB_int.end(), 0.0);
  std::fill(dN.begin(), dN.end(), 0.0);
  std::fill(risk_q.begin(), risk_q.end(), 0.0);
  std::fill(B_denom.begin(), B_denom.end(), 0.0);
  
  // Main backward recursion loop
  for (int j = K_q - 1; j >= 0; j--) {
    
    // Reset for each time point
    std::fill(B_numer.begin(), B_numer.end(), 0.0);
    std::fill(B_denom.begin(), B_denom.end(), 0.0);
    
    // Compute increments and risks
    for (int i = 0; i < N; i++) {
      dN(i, j) = (Q[i] == qtime[j]) ? -1 : 0;
      risk_q(i, j) = (Q[i] <= qtime[j] && qtime[j] < TT[i]) ? 1 : 0;
      
      double exp_Bint = std::exp(B_int(i, j+1)); // this is at t+
      
      for (int dim_w = 0; dim_w < dim_W; dim_w++) {
        B_numer[dim_w] += weights[i] * IPCW_weights(i,j) * W(i, dim_w) * exp_Bint * dN(i, j);
        for (int dim_z = 0; dim_z < dim_Z; dim_z++) {
          B_denom(dim_w, dim_z) += weights[i] * IPCW_weights(i,j) * W(i, dim_w) * Z(i, dim_z) * risk_q(i, j) * exp_Bint;
        }
      }
    }
    
    // Convert B_denom to arma matrix for pseudo-inverse
    arma::mat B_denom_arma(B_denom.begin(), dim_W, dim_Z, false);
    
    // Calculate pseudo-inverse safely
    bool pinv_success = false;
    arma::mat pinv_B_denom = safe_pinv(B_denom_arma, tol, pinv_success);
    
    if (!pinv_success) {
      Rcpp::Rcerr << "Pseudo-inverse failed at time index j = " << j << "\n";
      continue; // Skip this time step if pseudo-inverse failed
    }
    
    // Compute dB = pinv(B_denom) * B_numer
    arma::colvec B_numer_arma(B_numer.begin(), dim_W, false);
    arma::colvec dB_arma = pinv_B_denom * B_numer_arma;
    
    // Copy back dB result to Rcpp matrix
    for (int dim_z = 0; dim_z < dim_Z; dim_z++) {
      dB(j, dim_z) = dB_arma(dim_z);
    }
    
    // Update B_int and B
    for (int i = 0; i < N; i++) {
      double dB_int_sum = 0.0;
      for (int dim_z = 0; dim_z < dim_Z; dim_z++) {
        dB_int_sum += Z(i, dim_z) * dB(j, dim_z);
      }
      dB_int(i, j) = dB_int_sum;
      B_int(i, j) = B_int(i, j + 1) - dB_int(i, j);
      // B_int(i, j) = (j + 1 < K_q) ? (B_int(i, j + 1) - dB_int(i, j)) : -dB_int(i, j);
    }
    
    for (int dim_z = 0; dim_z < dim_Z; dim_z++) {
      B(j, dim_z) = B(j + 1, dim_z) - dB(j, dim_z);
      // if (j == K_q - 1) {
      //   B(j, dim_z) = -dB(j, dim_z);
      // } else {
      //   B(j, dim_z) = B(j + 1, dim_z) - dB(j, dim_z);
      // }
    }
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
    Named("qtime") = qtime
  );
}

