// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

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
  
  arma::mat y_std = y;
  arma::mat x_std = x;
  
  // Dimensions.
  int p = x.n_cols;
  
  y_std = (y - arma::mean(y)) / Stdev(y);
  for (int i = 0; i < p; i++) {
    arma::colvec w = x.col(i);
    x_std.col(i) = (w - arma::mean(w)) / Stdev(w);
  };
  
  // Output structure.
  return {y_std, x_std};
};


// ----------------------------------------------------------------------------
// Marginal correlations.
// ----------------------------------------------------------------------------

//' Marginal correlations
//'
//' Fits the standard OLS model.
//'
//' @param y (n x 1) Numeric vector.
//' @param x (n x p) Numeric matrix.
//' @return List containing the following:
//' \itemize{
//' \item{r: Correlations between y and the columns of x.}
//' \item{R: Correlations among the columns of x.}
//' }
//' @export
// [[Rcpp::export]]
SEXP Corr(
   const arma::colvec y, 
   const arma::mat x
) {
  
  // Standardize the data.
  XY std = Standardize(y, x);
  
  // Get standardized versions.
  const arma::colvec& y_std = std.y;
  const arma::mat& x_std = std.x;
  
  // Dimensions.
  const int n = y.size();
  
  // y-x correlations.
  const arma::colvec r = x_std.t() * y_std / n;
  
  // x-x correlations.
  const arma::mat R = x_std.t() * x_std / n;
  
 return Rcpp::List::create(
   Rcpp::Named("r") = r,
   Rcpp::Named("R") = R
 );
};


// ----------------------------------------------------------------------------
// Ordinary least squares.
// ----------------------------------------------------------------------------

//' Ordinary Least Squares
//'
//' Fits the standard OLS model.
//'
//' @param y (n x 1) Numeric vector.
//' @param x (n x p) Numeric matrix.
//' @param standard Standardize y and x?
//' @return List containing the following:
//' \itemize{
//' \item{beta: Regression coefficients.}
//' \item{v: Residual variance.}
//' \item{se: Standard errors.}
//' \item{z: Z-scores.}
//' \item{pval: P-values based on the chi2 distribution.}
//' }
//' @export
// [[Rcpp::export]]
SEXP OLS(
   const arma::colvec y, 
   const arma::mat x,
   const bool standard=true
) {
  
  // Dimensions.
  const int n = x.n_rows;
  const int p = x.n_cols;
  
  // Standardize.
  arma::colvec y_std = y;
  arma::mat x_std = x; 
  if (standard) {
    XY std = Standardize(y, x);
    y_std = std.y;
    x_std = std.x;
  }
  
  // Information.
  const arma::mat A = x_std.t() * x_std;
  
  // Estimate beta.
  const arma::colvec b = arma::solve(A, x_std.t() * y_std, arma::solve_opts::likely_sympd);
  
  // Calculate residuals.
  const arma::colvec eps = (y_std - x_std * b);
  
  // Scale.
  const double v = arma::as_scalar((eps.t() * eps) / (n - p));
  
  // Information.
  const arma::mat Ibb = A / v;
  
  // Standard errors.
  const arma::vec se = arma::sqrt(arma::diagvec(arma::pinv(Ibb)));
  
  // Z-scores.
  const arma::vec z = b / se;
  
  // P-values.
  Rcpp::Environment base("package:stats");
  Rcpp::Function pchisq = base["pchisq"]; 
  SEXP pval = pchisq(
    Rcpp::_["q"]=arma::pow(z, 2), Rcpp::_["df"]=1, Rcpp::_["lower.tail"]=false);
  
  return Rcpp::List::create(
    Rcpp::Named("beta") = b,
    Rcpp::Named("v") = v,
    Rcpp::Named("se") = se,
    Rcpp::Named("z") = z,
    Rcpp::Named("pval") = pval,
    Rcpp::Named("resid") = eps
  );
 };


// ----------------------------------------------------------------------------
// Pratt index.
// ----------------------------------------------------------------------------


//' Pratt Index
//'
//' Calculate the Pratt index and intermediates.
//'
//' @param y (n x 1) Numeric vector.
//' @param x (n x p) Numeric matrix.
//' @param standard Standardize y and x?
//' @return List containing the following:
//' \itemize{
//' \item{beta: Regression coefficients.}
//' \item{r: Correlations between y and the columns of x.}
//' \item{R: Correlations among the columns of x.}
//' \item{pratt: Pratt indices for the columns of x.}
//' \item{var_exp: Variance in y explained by x.}
//' }
//' @export
// [[Rcpp::export]]
SEXP PrattIndex(
    const arma::colvec y, 
    const arma::mat x,
    const bool standard=true
) {

  // Standardize.
  arma::colvec y_std = y;
  arma::mat x_std = x; 
  if (standard) {
    XY std = Standardize(y, x);
    y_std = std.y;
    x_std = std.x;
  }
  
  // Correlations.
  const int n = x.n_rows;
  const int p = x.n_cols;
  const arma::colvec r = x_std.t() * y_std / n;
  const arma::mat R = x_std.t() * x_std / n;
  
  // Regression coefficients.
  // Note that beta = R^{-1} %*% r.
  const arma::colvec beta = arma::solve(R, r, arma::solve_opts::likely_sympd);
  
  // Pratt Index.
  const arma::colvec pratt = beta % r; 
  
  // Variation explained.
  const double var_exp = arma::as_scalar(y_std.t() * x_std * beta / n);
  
  // Calculate variance.
  arma::colvec pratt_var = arma::zeros(p);
  const arma::mat Rinv = arma::pinv(R);
  for(int i = 0; i < p; i++) {
    double sigma = r(i) * Rinv(i, i) * r(i) + beta(i) * R(i, i) * beta(i) + 2 * beta(i) * r(i);
    pratt_var(i) = sigma / n; 
  }
  
  // Output.
  return Rcpp::List::create(
    Rcpp::Named("beta") = beta,
    Rcpp::Named("pratt") = pratt,
    Rcpp::Named("r") = r,
    Rcpp::Named("R") = R,
    Rcpp::Named("pratt_var") = pratt_var,
    Rcpp::Named("var_exp") = var_exp
  );
};
