// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath> 


//' Matrix Square Root
//' @param S (p x p) numeric matrix. 
//' @param eps Scalar minimum eigenvalue.
//' @export
// [[Rcpp::export]]
SEXP MatrixSqrt(
    const arma::mat& S,
    double eps = 1e-8
) {
  
  // Check that S is symmetric.
  if (!S.is_symmetric()) {
    Rcpp::Rcout << "Input matrix is not symmetric."<< std::endl;
    return R_NilValue;
  }

  arma::mat out;
  
  if (S.is_sympd()) {
    
    out = sqrtmat_sympd(S);
    
  } else {
    
    arma::vec eigenval;
    arma::mat eigenvec;
    arma::eig_sym(eigenval, eigenvec, S);
    
    int p = eigenval.n_elem;
    arma::vec sqrt_eigenval = arma::zeros(p);
    for (int i = 0; i < p; i++) {
      double d = eigenval(i);
      if (d < eps) {
        if (d >= 0) {
          d = 0;
        } else {
          Rcpp::Rcout << "Eigenvalue: " << d << " is negative."<< std::endl;
          return R_NilValue;
        }
      }
      sqrt_eigenval(i) = std::sqrt(d);
    }
    arma::mat D = arma::diagmat(sqrt_eigenval);
    out = eigenvec * D * eigenvec.t();
  }
  return Rcpp::wrap(out);
}

