// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>

//' Matrix Square Root
//' @param S (p x p) symmetric numeric matrix.
//' @param eps Minimum eigenvalue threshold (smaller non-negative eigenvalues set to 0).
// [[Rcpp::export]]
SEXP MatrixSqrt(const arma::mat& S, double eps) {
  if (!S.is_symmetric()) {
    Rcpp::stop("Input matrix is not symmetric.");
  }
  arma::mat out;
  if (S.is_sympd()) {
    out = sqrtmat_sympd(S);
  } else {
    arma::vec eigenval;
    arma::mat eigenvec;
    arma::eig_sym(eigenval, eigenvec, S);
    const arma::uword p = eigenval.n_elem;
    arma::vec sqrt_eigenval(p, arma::fill::zeros);
    for (arma::uword i = 0; i < p; ++i) {
      double d = eigenval(i);
      if (d < 0) Rcpp::stop("Matrix has negative eigenvalue.");
      if (d >= eps) sqrt_eigenval(i) = std::sqrt(d);
    }
    out = eigenvec * arma::diagmat(sqrt_eigenval) * eigenvec.t();
  }
  return Rcpp::wrap(out);
}

