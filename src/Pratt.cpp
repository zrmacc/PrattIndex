// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// ----------------------------------------------------------------------------
// Pratt index.
// ----------------------------------------------------------------------------

// XY structure.
// @param y (n x 1) Numeric vector.
// @param X (n x p) Numeric matrix.
struct XY {
  arma::colvec y;
  arma::mat x;
};


// Standard Deviation
//
// @param x Numeric vector.
// @return double.
double Stdev(const arma::colvec x) {
  const int n = x.size();
  const double sigma2 = arma::accu(arma::square(x - arma::mean(x))) / n;
  const double sigma = std::sqrt(sigma2);
  return sigma;
};


// Standardize
// 
// @param y (n x 1) Numeric vector.
// @param x (n x p) Numeric matrix.
// @return XY
XY Standardize(
    const arma::colvec y,
    const arma::mat x
) {
  arma::colvec y_std = y;
  arma::mat x_std = x;
  
  // Number of columns in x.
  int p = x.n_cols;
  
  y_std = (y - arma::mean(y)) / Stdev(y);
  for (int i = 0; i < p; i++) {
    arma::colvec w = x.col(i);
    x_std.col(i) = (w - arma::mean(w)) / Stdev(w);
  }
  
  // Output structure.
  XY out;
  out.y = y_std;
  out.x = x_std;
  return out;
};


//' Pratt Index
//'
//' Calculate the Pratt index and intermediates.
//'
//' @param y (n x 1) Numeric vector.
//' @param x (n x p) Numeric matrix. 
//' @param Standardize y and x?
//' @return List containing the following:
//' \itemize{
//' \item{beta: (p x 1) vector of regression coefficient.}
//' \item{var_beta: (p x p) sampling variance of the regression coefficients.}
//' \item{r: (p x 1) vector of correlations between y and each column of x..}
//' \item{var_r: (p x p) sampling variance of r.}
//' \item{R: (p x p) correlation marix among the columns of x.}
//' \item{pratt: (p x 1) Pratt indices, defined as the element-wise product of beta and r.}
//' \item{var_pratt: (p x 1) sampling variance of Pratt indices.}
//' \item{resid_var: scalar residual variance of y conditional on x.}
//' \item{var_exp: scalar variance in y explained by x.}
//' }
//' @export
// [[Rcpp::export]]
SEXP PrattIndex(
    arma::colvec y, 
    arma::mat x,
    const bool standardize = true
) {

  // Dimensions.
  const int n = x.n_rows;
  const int p = x.n_cols;
  
  // Standardize.
  if (standardize) {
    XY xy_std = Standardize(y, x);
    y = xy_std.y;
    x = xy_std.x;
  }
  
  // Estimate regression coefficients.
  const arma::mat R = x.t() * x / double(n);
  const arma::colvec r = x.t() * y / double(n);
  const arma::colvec beta = arma::solve(R, r, arma::solve_opts::likely_sympd);
  
  // Estimate residual variance.
  const arma::colvec eps = (y - x * beta);
  const double v = arma::as_scalar((eps.t() * eps) / double(n));
  
  // Estimate marginal correlations.
  const arma::mat Rinv = arma::pinv(R);
  
  // Variance of beta.
  const arma::mat var_beta = v * Rinv / double(n);
  
  // Variance of r.
  const arma::mat var_r = v * R / double(n);
  
  // Pratt Index.
  const arma::colvec pratt = beta % r; 
  
  // Variation explained.
  const double var_exp = arma::as_scalar(y.t() * x * beta / n);
  
  // Calculate variance.
  arma::colvec var_pratt = arma::zeros(p);
  for(int i = 0; i < p; i++) {
    double sigma = r(i) * Rinv(i, i) * r(i) + beta(i) * R(i, i) * beta(i) + 2 * beta(i) * r(i) + Rinv(i, i) + 1;
    var_pratt(i) = v * sigma / double(n); 
  }
  
  // Output.
  return Rcpp::List::create(
    Rcpp::Named("beta") = beta,
    Rcpp::Named("var_beta") = var_beta,
    Rcpp::Named("r") = r,
    Rcpp::Named("var_r") = var_r,
    Rcpp::Named("R") = R,
    Rcpp::Named("pratt") = pratt,
    Rcpp::Named("var_pratt") = var_pratt,
    Rcpp::Named("resid_var") = v,
    Rcpp::Named("var_exp") = var_exp
  );
};
