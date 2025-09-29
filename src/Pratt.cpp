// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath> 

// ----------------------------------------------------------------------------
// Pratt index.
// ----------------------------------------------------------------------------

// Output of OLS
struct OLSFit {
  arma::colvec beta;
  arma::mat var_beta;
  double var_resid;
};

// Fit Ordinary Least Squares
// @param X Design matrix.
// @param y Response vector.
OLSFit OLSCpp(const arma::mat& X, const arma::colvec& y) {
  const arma::uword n = X.n_rows;
  const arma::uword p = X.n_cols;
  const arma::mat XtX  = X.t() * X;
  const arma::mat XtXi = arma::pinv(XtX);
  const arma::colvec beta =
    arma::solve(XtX, X.t() * y, arma::solve_opts::likely_sympd);
  const arma::colvec resid = y - X * beta;
  const arma::uword df = n - p;
  const double var_resid = arma::as_scalar((resid.t() * resid) / double(df));
  const arma::mat var_beta = var_resid * XtXi;
  OLSFit out;
  out.beta = beta;
  out.var_beta = var_beta;
  out.var_resid = var_resid;
  return out;
}


//' Pratt Index
//' 
//' @param y (n x 1) numeric
//' @param g (n x 1) numeric
//' @param e (n x 1) numeric
//' @return List containing:
//'   * beta: Regression coefficients (b_G, b_E, b_H) from the joint model.
//'   * beta_null: Regression coefficient (b_G, b_E) from the model that omits H.
//'   * cov_xx: Covariance matrix of (G, E, H).
//'   * kappa: Pratt index
//'   * grad_kappa: Gradient of kappa with respect to beta.
//'   * var_beta: Covariance matrix of beta.
//'   * var_beta_null: covariance matrix of beta from the model that omits H.
//'   * var_kappa: Variance of Pratt indices.
//'   * var_resid: Residual variance from the model that includes H.
//'   * var_resid_null: Residual variance from the model the omits H.
// [[Rcpp::export]]
SEXP PrattIndexCpp(arma::colvec y, arma::colvec g, arma::colvec e) {
  
  // Design matrix ------------------------------------------------------------
  const arma::uword n = y.n_elem;
  const arma::vec intercept = arma::ones(n);
  const arma::vec h = g % e;
  const arma::mat X0 = arma::join_rows(intercept, g, e);
  const arma::mat X1 = arma::join_rows(intercept, g, e, h);
  const arma::mat X2 = arma::join_rows(g, e, h);
  
  // OLS fits -----------------------------------------------------------------
  OLSFit ols_without_h = OLSCpp(X0, y);
  OLSFit ols_with_h = OLSCpp(X1, y);

  const arma::colvec beta = ols_with_h.beta.subvec(1, 3);
  const arma::mat var_beta = ols_with_h.var_beta.submat(1, 1, 3, 3);
  const double var_resid = ols_with_h.var_resid;
  
  // Covariance ---------------------------------------------------------------
  const arma::mat cov_xx = arma::cov(X2);
  const arma::vec cov_yx = cov_xx * beta; 
  const double var_y = arma::as_scalar(beta.t() * cov_yx) + var_resid;
  
  // Pratt index --------------------------------------------------------------
  const arma::vec kappa = (beta % cov_yx) / var_y;
  
  // Gradients ----------------------------------------------------------------
  const arma::uword n_beta = beta.n_elem;
  arma::mat grad_kappa(n_beta, n_beta, arma::fill::zeros);
  
  for (arma::uword j = 0; j < n_beta; ++j) {
    arma::vec term(n_beta, arma::fill::zeros);
    term(j) = cov_yx(j);
    term += beta(j) * cov_xx.row(j).t();
    grad_kappa.col(j) = (term - 2.0 * kappa(j) * cov_yx) / var_y;
  }
  
  // Delta-method -------------------------------------------------------------
  arma::vec var_kappa(n_beta);
  for (arma::uword j = 0; j < n_beta; j++) {
    var_kappa(j) = arma::as_scalar(grad_kappa.col(j).t() * var_beta * grad_kappa.col(j));
  }
  
  // Results ------------------------------------------------------------------
  return Rcpp::List::create(
   Rcpp::Named("beta") = beta,
   Rcpp::Named("beta_null") = arma::colvec(ols_without_h.beta.subvec(1, 2)),
   Rcpp::Named("cov_xx") = cov_xx,
   Rcpp::Named("kappa") = kappa,
   Rcpp::Named("grad_kappa") = grad_kappa,
   Rcpp::Named("n") = n,
   Rcpp::Named("var_beta") = var_beta,
   Rcpp::Named("var_beta_null") = ols_without_h.var_beta.submat(1, 1, 2, 2),
   Rcpp::Named("var_kappa") = var_kappa,
   Rcpp::Named("var_resid") = var_resid,
   Rcpp::Named("var_resid_null") = ols_without_h.var_resid
  );
 }
