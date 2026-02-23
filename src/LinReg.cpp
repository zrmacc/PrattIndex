// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>

//' Ordinary Least Squares
//' @param y n x 1 numeric vector.
//' @param X n x p numeric matrix.
//' @return List with \code{beta}, \code{resid_var}, \code{var_beta}.
// [[Rcpp::export]]
SEXP FitOLS(const arma::colvec& y, const arma::mat& X) {
  const arma::uword n = y.n_elem;
  const arma::uword p = X.n_cols;
  const arma::mat A = X.t() * X;
  const arma::vec beta = arma::solve(A, X.t() * y, arma::solve_opts::likely_sympd);
  const arma::vec eps = y - X * beta;
  const double resid_var = arma::as_scalar(eps.t() * eps / double(n - p));
  arma::mat Ainv;
  if (!arma::inv_sympd(Ainv, A)) Ainv = arma::pinv(A);
  return Rcpp::List::create(
    Rcpp::Named("beta") = beta,
    Rcpp::Named("resid_var") = resid_var,
    Rcpp::Named("var_beta") = resid_var * Ainv
  );
}

