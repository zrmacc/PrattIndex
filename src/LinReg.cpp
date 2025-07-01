// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath> 


//' Ordinary Least Squares
//' 
//' Fits the standard OLS model.
//' 
//' @param y n x 1 Numeric vector.
//' @param X n x p Numeric matrix.
//' 
//' @return List containing the following:
//' \item{Beta}{Regression coefficient.}
//' \item{V}{Outcome variance.}
//' \item{Ibb}{Information matrix for beta.}
//' \item{Resid}{Outcome residuals.}
//' @export
// [[Rcpp::export]]
SEXP FitOLS(const arma::colvec y, const arma::mat X){
  
  // Observation and features.
  const int n = y.size();
  const int p = X.n_cols;
  
  // Estimate beta.
  const arma::mat A = X.t()*X;
  const arma::vec beta = arma::solve(A, X.t() * y, arma::solve_opts::likely_sympd);
  
  // Residual variance.
  const arma::vec eps = (y - X * beta);
  const double resid_var = arma::as_scalar(eps.t()*eps/(n-p));
  
  // Variance of beta.
  const arma::mat var_beta = resid_var * arma::pinv(A);
  
  // Output.
  return Rcpp::List::create(
    Rcpp::Named("beta")=beta,
    Rcpp::Named("resid_var")=resid_var,
    Rcpp::Named("var_beta")=var_beta
  );
}

