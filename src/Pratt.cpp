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
//' \item{beta: (4 x 1) vector of regression coefficients (intercept, G, E, H).}
//' \item{kappa: (3 x 1) vector of Pratt indices (G, E, H).}
//' \item{grad_kappa: (3 x 3) matrix of Pratt index gradients with respect to beta. 
//' Columns corresponds to kappa (G, E, H), rows to beta (G, E, H).}
//' \item{rho: (3 x 1) vector of correlations (G, E, H).}
//' \item{var_beta: (4 x 4) covariance matrix of beta.}
//' \item{var_kappa: (3 x 1) vector of Pratt index variances (G, E, H).}
//' \item{var_resid: Residual variance from regression of y ~ (intercept, G, E, H).}
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
  
  const double beta_g = beta(1);
  const double beta_e = beta(2);
  const double beta_h = beta(3);
  
  // Estimate residual variance.
  const arma::colvec eps = (y - x * beta);
  const double v = arma::as_scalar((eps.t() * eps) / double(n - p));

  // Variance of beta.
  const arma::mat var_beta = v * XtXi;
  
  // Moments.
  const double mu_g = arma::mean(g);
  const double var_g = arma::var(g);
  const double std_g = std::sqrt(var_g);
  
  const double mu_e = arma::mean(e);
  const double var_e = arma::var(e);
  const double std_e = std::sqrt(var_e);
  
  // const double mu_h = arma::mean(h);
  const double var_h = arma::var(h);
  const double std_h = std::sqrt(var_h);
  
  // Model-based variance of Y.
  const double var_y = beta_g * beta_g * var_g + 
    beta_e * beta_e * var_e +
    beta_h * beta_h * var_h + 
    2.0 * beta_g * beta_h * var_g * mu_e + 
    2.0 * beta_e * beta_h * var_e * mu_g + 
    v;
  const double std_y = std::sqrt(var_y);
  
  // Correlations.
  arma::vec rho(3);
  rho(0) = (std_g / std_y) * (beta_g + mu_e * beta_h);
  rho(1) = (std_e / std_y) * (beta_e + mu_g * beta_h);
  rho(2) = (var_g * mu_e * beta_g + var_e * mu_g * beta_e + var_h * beta_h) / (std_h * std_y);
  
  // Pratt index.
  arma::vec kappa(3);
  kappa(0) = beta_g * rho(0);
  kappa(1) = beta_e * rho(1);
  kappa(2) = beta_h * rho(2);
  
  // Gradients of std_y wrt beta.
  arma::vec grad_std_y(3);
  grad_std_y(0) = std_g * rho(0);
  grad_std_y(1) = std_e * rho(1);
  grad_std_y(2) = std_h * rho(2);
  
  // -------- Gradients of lambdas --------
  // lambda_G = beta_G * sigma_G * (beta_G + mu_E * beta_H)
  arma::colvec grad_lambda_g(3, arma::fill::zeros);
  grad_lambda_g(0) = std_g * (2.0 * beta_g + mu_e * beta_h);
  grad_lambda_g(2) = std_g * mu_e * beta_g;

  // lambda_E = beta_E * sigma_E * (beta_E + mu_G * beta_H)
  arma::colvec grad_lambda_e(3, arma::fill::zeros);
  grad_lambda_e(1) = std_e * (2.0 * beta_e + mu_g * beta_h);
  grad_lambda_e(2) = std_e * mu_g * beta_e;
  
  // lambda_H = (beta_H / sigma_H) * (var_g * mu_E * beta_G + var_e * mu_G * beta_E + var_h * beta_H)
  arma::colvec grad_lambda_h(3, arma::fill::zeros);
  grad_lambda_h(0) = (var_g * mu_e * beta_h) / std_h;
  grad_lambda_h(1) = (var_e * mu_g * beta_h) / std_h;
  grad_lambda_h(2) = (var_g * mu_e * beta_g + var_e * mu_g * beta_e + 2.0 * var_h * beta_h) / std_h;
  
  // -------- Gradients of kappas (cols) --------
  // ∂κ_j/∂β = (1/σ_Y) * (∂λ_j/∂β - κ_j * ∂σ_Y/∂β)
  arma::mat grad_kappa(3, 3);
  grad_kappa.col(0) = (grad_lambda_g - kappa(0) * grad_std_y) / std_y; // dκ_G/dβ
  grad_kappa.col(1) = (grad_lambda_e - kappa(1) * grad_std_y) / std_y; // dκ_E/dβ
  grad_kappa.col(2) = (grad_lambda_h - kappa(2) * grad_std_y) / std_y; // dκ_H/dβ
  
  // Variance of Pratt index.
  const arma::mat sigma = var_beta.submat(1, 1, 3, 3);
  arma::vec var_kappa(3);
  var_kappa(0) = arma::as_scalar( grad_kappa.col(0).t() * sigma * grad_kappa.col(0) ); // Var(κ_G)
  var_kappa(1) = arma::as_scalar( grad_kappa.col(1).t() * sigma * grad_kappa.col(1) ); // Var(κ_E)
  var_kappa(2) = arma::as_scalar( grad_kappa.col(2).t() * sigma * grad_kappa.col(2) ); // Var(κ_H)
  
  // Output.
  return Rcpp::List::create(
    Rcpp::Named("beta") = beta,
    Rcpp::Named("grad_kappa") = grad_kappa,
    Rcpp::Named("kappa") = kappa,
    Rcpp::Named("rho") = rho,
    Rcpp::Named("var_beta") = var_beta,
    Rcpp::Named("var_kappa") = var_kappa,
    Rcpp::Named("var_resid") = v
  );
}
