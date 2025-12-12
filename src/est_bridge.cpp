# include <RcppArmadillo.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//




// [[Rcpp::export]]
SEXP qbridge(NumericVector TT,
              NumericVector Q, 
              NumericVector qtime, // sorted left truncation time -- the jumps, !! including tau2 at the end, see if we still need the input tau2?
              NumericMatrix Z, // (1,W1,Z)
              NumericMatrix W, // (1,W2,Z)
              NumericVector weights){
  
  // variables for estimation
  int N = TT.size(), K_q = qtime.size(), dim_Z = Z.cols(), dim_W = W.cols();
  NumericVector B_numer(dim_W);
  NumericMatrix B(K_q, dim_Z), dB(K_q, dim_Z), B_int(N, K_q), dB_int(N, K_q), B_denom(dim_W, dim_Z);  // B_int = B(t)(1,W1,Z)
  NumericMatrix dN(N, K_q),  risk_q(N, K_q);
  double del = 0.001;

  
  for (int j = K_q-1; j >= 0; j--) {
    // initiate the vector B_numer and the matrix B_denom to be all zeros
    for (int dim_w = 0; dim_w < dim_W; dim_w++) {
      B_numer[dim_w] = 0;
      for (int dim_z = 0; dim_z < dim_Z; dim_z++) {
        B_denom(dim_w, dim_z) = 0;
      }
    }

    for (int i = 0; i < N; i++) {
      // estimation
      dN(i, j) = (Q[i] == qtime[j]) ? 1:0;
      dN(i, j) = - dN(i, j);  // jumps -1 if Q[i] = qtime[j]
      risk_q(i, j) = (Q[i] <= qtime[j] && qtime[j] < TT[i]) ? 1:0;
      for (int dim_w = 0; dim_w < dim_W; dim_w++) {
        B_numer[dim_w] += weights[i] * W(i, dim_w) * exp(B_int(i, j)) * dN(i, j);
        for (int dim_z = 0; dim_z < dim_Z; dim_z++) {
          B_denom(dim_w, dim_z) += weights[i] * W(i, dim_w) * Z(i, dim_z) * risk_q(i, j) * exp(B_int(i, j));
        }
      }
    }
    


    // solve linear system using RcppArmadillo
    try
    {
      arma::mat X(B_denom.begin(), dim_Z, dim_W, false);
      if (arma::rank(X) < std::min(dim_W, dim_Z)) {
        std::cerr << "Warning: X is rank-deficient, pseudo-inverse may be unstable." << std::endl;
      }
      arma::mat tmp_pinv = arma::pinv(X, del);
      arma::mat tmp_res = tmp_pinv * as<arma::vec>(B_numer);
      for (int dim_w = 0; dim_w < dim_W; dim_w++) {
        dB(j, dim_w) = tmp_res(dim_w, 0);
      }
    } catch (...)  
    {
       
    }
      
      

    
    // estimation
    for (int i = 0; i < N; i++) {
      for (int dim_z = 0; dim_z < dim_Z; dim_z++) {
        dB_int(i, j) += Z(i, dim_z) * dB(j, dim_z);
      }
      B_int(i, j) = B_int(i, j+1) - dB_int(i, j);
    }
    

    for (int dim_z = 0; dim_z < dim_Z; dim_z++) {
      if (j == K_q - 1) {
        B(j, dim_z) = 0 - dB(j, dim_z);
      }
      else{
        B(j, dim_z) = B(j + 1, dim_z) - dB(j, dim_z);
      }
    }
    
    
  }

  
  
  
  return List::create(
    Named("dN_Q") = dN,
    Named("risk_q") = risk_q,
    Named("dB") = dB,
    Named("B") = B,
    Named("B_int") = B_int,
    Named("dB_int") = dB_int,
    Named("B_numer") = B_numer,
    Named("B_denom") = B_denom,
    Named("qtime") = qtime);
}