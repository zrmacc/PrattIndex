// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>

// ----------------------------------------------------------------------------
// Pratt index.
// ----------------------------------------------------------------------------

struct OLSFit {
  arma::colvec beta;
  arma::mat var_beta;
  double var_resid;
};

// Fit OLS: beta and var_beta from (X'X)^{-1}, residual variance from (n - p).
OLSFit OLSCpp(const arma::mat& X, const arma::colvec& y) {
  const arma::uword n = X.n_rows;
  const arma::uword p = X.n_cols;
  const arma::mat XtX = X.t() * X;
  const arma::colvec beta =
    arma::solve(XtX, X.t() * y, arma::solve_opts::likely_sympd);
  const arma::colvec resid = y - X * beta;
  const arma::uword df = n - p;
  const double var_resid = arma::as_scalar((resid.t() * resid) / double(df));
  arma::mat XtXi;
  if (!arma::inv_sympd(XtXi, XtX)) XtXi = arma::pinv(XtX);
  OLSFit out;
  out.beta = beta;
  out.var_beta = var_resid * XtXi;
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
SEXP PrattIndexCpp(const arma::colvec& y, const arma::colvec& g, const arma::colvec& e) {
  
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
  
  const arma::mat cov_xx = arma::cov(X2);
  arma::mat cov_xx_inv;
  if (!arma::inv_sympd(cov_xx_inv, cov_xx)) cov_xx_inv = arma::pinv(cov_xx);
  const arma::vec cov_yx = cov_xx * beta; 
  
  // Variance of y ------------------------------------------------------------
  // Model-based version:
  // const double var_y = arma::as_scalar(beta.t() * cov_yx) + var_resid;
  
  const arma::vec y_c = y - arma::mean(y);
  const double var_y = arma::as_scalar((y_c.t() * y_c) / double(n - 1));
  
  // Pratt index --------------------------------------------------------------
  const arma::vec kappa = (beta % cov_yx) / var_y;
  
  // Gradients ----------------------------------------------------------------
  const arma::uword n_beta = beta.n_elem;
  arma::mat grad_kappa(n_beta, n_beta, arma::fill::zeros);
  
  for (arma::uword j = 0; j < n_beta; ++j) {
    arma::vec term(n_beta, arma::fill::zeros);
    term(j) = cov_yx(j);
    term += beta(j) * cov_xx.row(j).t();
    grad_kappa.col(j) = term / var_y;
  }
  
  // Delta-method -------------------------------------------------------------
  arma::vec var_kappa(n_beta);
  for (arma::uword j = 0; j < n_beta; ++j) {
    var_kappa(j) = arma::as_scalar(grad_kappa.col(j).t() * var_beta * grad_kappa.col(j));
  }
  
  // Results ------------------------------------------------------------------
  return Rcpp::List::create(
   Rcpp::Named("beta") = beta,
   Rcpp::Named("beta_null") = arma::colvec(ols_without_h.beta.subvec(1, 2)),
   Rcpp::Named("cov_xx") = cov_xx,
   Rcpp::Named("cov_xx_inv") = cov_xx_inv,
   Rcpp::Named("kappa") = kappa,
   Rcpp::Named("grad_kappa") = grad_kappa,
   Rcpp::Named("n") = n,
   Rcpp::Named("var_beta") = var_beta,
   Rcpp::Named("var_beta_null") = ols_without_h.var_beta.submat(1, 1, 2, 2),
   Rcpp::Named("var_kappa") = var_kappa,
   Rcpp::Named("var_resid") = var_resid,
   Rcpp::Named("var_resid_null") = ols_without_h.var_resid,
   Rcpp::Named("var_y") = var_y
  );
}

//' Pratt Index Influence Function
//' 
//' Computes the plug-in influence function for the Pratt index components
//' corresponding to X = (G, E, H) with H = G * E.
//'
//' @param y (n x 1) numeric
//' @param g (n x 1) numeric
//' @param e (n x 1) numeric
//' @return (n x 3) numeric matrix psi, where psi.row(i) is IF_kappa for obs i.
// [[Rcpp::export]]
SEXP PrattInfluenceCpp(const arma::colvec& y, const arma::colvec& g, const arma::colvec& e) {
  const arma::uword n = y.n_elem;
  if (g.n_elem != n || e.n_elem != n) {
    Rcpp::stop("Input vectors y, g, e must have the same length.");
  }
  if (n < 3) {
    Rcpp::stop("Need n >= 3 to form sample covariances with (n - 1) denominator.");
  }
  
  // Design matrix (no intercept; centered moments handle means) ---------------
  const arma::vec h = g % e;                         // H = G * E
  const arma::mat X = arma::join_rows(g, e, h);      // n x 3
  
  // Centered variables --------------------------------------------------------
  const arma::rowvec mean_x = arma::mean(X, 0);      // 1 x 3
  const double mean_y = arma::mean(y);               // scalar
  
  const arma::mat Xc = X.each_row() - mean_x;        // n x 3
  const arma::vec yc = y - mean_y;                   // n x 1
  
  // Sample moments (unbiased, denominator n-1) --------------------------------
  const double denom = double(n - 1);
  const arma::mat SigmaXX = (Xc.t() * Xc) / denom;            // 3 x 3
  const arma::vec SigmaXY = (Xc.t() * yc) / denom;            // 3 x 1
  const double SigmaYY = arma::as_scalar((yc.t() * yc) / denom); // scalar
  
  // Inverses / coefficients ---------------------------------------------------
  const arma::mat SigmaXX_inv = arma::pinv(SigmaXX);          // 3 x 3
  const arma::vec beta = SigmaXX_inv * SigmaXY;               // 3 x 1
  const arma::vec kappa = (beta % SigmaXY) / SigmaYY;         // 3 x 1
  
  // Precompute diagonals ------------------------------------------------------
  const arma::mat Dxy = arma::diagmat(SigmaXY);               // 3 x 3
  const arma::mat Db  = arma::diagmat(beta);                  // 3 x 3
  
  // Output --------------------------------------------------------------------
  arma::mat psi(n, 3, arma::fill::zeros);
  
  for (arma::uword i = 0; i < n; ++i) {
    
    const arma::vec u = Xc.row(i).t();          // 3 x 1  (tilde X_i)
    const double v = yc(i);                     // scalar (tilde Y_i)
    
    const arma::mat DeltaXX = (u * u.t()) - SigmaXX;    // 3 x 3
    const arma::vec DeltaXY = (u * v) - SigmaXY;        // 3 x 1
    const double DeltaYY = (v * v) - SigmaYY;           // scalar
    
    const arma::vec term1 = (Dxy * SigmaXX_inv + Db) * DeltaXY; // 3 x 1
    const arma::vec term2 = Dxy * SigmaXX_inv * DeltaXX * beta; // 3 x 1
    
    const arma::vec if_kappa =
      (term1 - term2) / SigmaYY - (kappa * (DeltaYY / SigmaYY));
    
    psi.row(i) = if_kappa.t();
  }
  
  return Rcpp::wrap(psi);
}
