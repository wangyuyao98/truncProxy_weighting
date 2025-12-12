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
SEXP pinv(NumericMatrix B_denom, double del) {
  
  // Get matrix dimensions
  int dim_row = B_denom.nrow();  // Number of rows
  int dim_col = B_denom.ncol();  // Number of columns
  arma::mat tmp_pinv;
  arma::mat X;
  
  // Solve linear system using RcppArmadillo
  try {
    arma::mat X(B_denom.begin(), dim_row, dim_col, false);  // Correct order
    tmp_pinv = arma::pinv(X, del);
    
    return List::create(
      Named("pinv") = tmp_pinv,
      Named("X") = X
    );
    
  } catch (std::exception &ex) {
    Rcpp::Rcerr << "Error: " << ex.what() << std::endl;
    return List::create(Named("error") = "Matrix inversion failed.");
  } catch (...) {
    Rcpp::Rcerr << "Unknown error in pinv()." << std::endl;
    return List::create(Named("error") = "Unknown error occurred.");
  }
  
}






