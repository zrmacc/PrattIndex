// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath> 

// ----------------------------------------------------------------------------
// Pratt index.
// ----------------------------------------------------------------------------

//' Pratt Index
//'
//' Calculate the Pratt index and intermediates.
//'
//' @param y (n x 1) Numeric vector.
//' @param g (n x 1) Genotype vector.
//' @param e (n x 1) Environment vector.
//' @return List containing the following:
//' \itemize{
//' \item{beta: (3 x 1) vector of regression coefficient.}
//' \item{pratt: Scalar Pratt index of the interaction term.}
//' \item{sigma2: Scalar residual variance from the joint model.}
//' \item{var_beta: (3 x 3) sampling variance of the regression coefficients.}
//' \item{var_exp: Scalar variance in y explained by x.}
//' \item{var_pratt: Scalar sampling variance of the pratt index.}
//' }
//' @export
// [[Rcpp::export]]
SEXP PrattIndex(
    arma::colvec y, 
    arma::colvec g,
    arma::colvec e 
) {

  // Dimensions.
  const int n = y.n_elem;
  const arma::vec intercept = arma::ones(n); 
  
  const arma::vec h = g % e; 
  const arma::mat x = arma::join_rows(intercept, g, e, h);
  const int p = x.n_cols;
  
  // Estimate beta.
  const arma::mat XtX = x.t() * x;
  const arma::mat XtXi = arma::pinv(XtX);
  
  // Estimate regression coefficients.
  const arma::colvec beta = arma::solve(XtX, x.t() * y, arma::solve_opts::likely_sympd);
  
  // Estimate residual variance.
  const arma::colvec eps = (y - x * beta);
  const double v = arma::as_scalar((eps.t() * eps) / double(n - p));

  // Variance of beta.
  const arma::mat var_beta = v * XtXi;
  
  // Pratt index.
  const double beta_g = beta(1);
  const double beta_e = beta(2);
  const double beta_h = beta(3);
  
  const double mu_g = arma::mean(g);
  const double var_g = arma::var(g);
  
  const double mu_e = arma::mean(e);
  const double var_e = arma::var(e);
  
  const double var_h = arma::var(h);
  const double var_y = arma::var(y);
  
  const double gamma_h = var_g * mu_e * beta_g + var_e * mu_g * beta_e + var_h * beta_h;
  const double pratt = beta_h * gamma_h / std::sqrt(var_h * var_y);
  
  const double var_pratt = var_beta(3, 3) * std::pow(
    var_g * mu_e * beta_g + var_e * mu_g * beta_e , 2) / (var_h * var_y);
  
  // Variation explained.
  const double var_exp = arma::as_scalar(y.t() * x * beta / n);
  
  // Output.
  return Rcpp::List::create(
    Rcpp::Named("beta") = beta,
    Rcpp::Named("pratt") = pratt,
    Rcpp::Named("sigma2") = v,
    Rcpp::Named("var_beta") = var_beta,
    Rcpp::Named("var_exp") = var_exp,
    Rcpp::Named("var_pratt") = var_pratt
  );
}
