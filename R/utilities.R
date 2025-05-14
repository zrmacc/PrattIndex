# Purpose: Utility functions for Pratt Index
# Updated: 25-05-14

#' Expected values
#' 
#' Calculates expected values under the generative model:
#' \deqn{Y = G\beta_{G} + E\beta_{E} + H\beta_{H} + \epsilon}
#' 
#' @param beta_g Joint coefficient for G.
#' @param mu_g Mean of G.
#' @param var_g Variance of G.
#' @param beta_e Joint coefficient for E.
#' @param mu_e Mean of E.
#' @param var_e Variance of E.
#' @param beta_h Joint coefficient for H = G x E.
#' @param var_resid Residual variance.
#' @return List of named expected values.
#' @export
Expected <- function(
  beta_g,
  mu_g,
  var_g,
  beta_e,
  mu_e,
  var_e,
  beta_h,
  mu_h,
  var_h,
  var_resid 
) {
  
  # Moments of H.
  mu_h <- mu_g * mu_e
  var_h <- var_g * var_e + (mu_g^2) * var_e + (mu_e^2) * var_g
  
  # Moments of Y.
  mu_y <- mu_g * beta_g + mu_e * beta_e + mu_h * beta_h
  var_y <- var_g * (beta_g^2) + var_e * (beta_e^2) + var_e * (beta_h^2) + 
    2 * (var_g * mu_e) * beta_g * beta_h + 
    2 * (var_e * mu_g) * beta_e * beta_h + 
    var_resid

  # Correlations.
  r_g <- sqrt(var_g) * (beta_g + mu_e * beta_h)
  r_e <- sqrt(var_e) * (beta_e + mu_g * beta_h)
  r_h <- (var_g * mu_e * beta_g + var_e * mu_g * beta_e + var_h * beta_h) / sqrt(var_h * var_y)
  
  # Pratt indices.
  pratt_g <- beta_g * r_g
  pratt_e <- beta_e * r_e
  pratt_h <- beta_h * r_h
  
  # Output.
  out <- list(
    mu_h = mu_h,
    var_h = var_h,
    r_g = r_g,
    r_e = r_e,
    r_h = r_h,
    pratt_g = pratt_g,
    pratt_e = pratt_e,
    pratt_h = pratt_h,
    mu_y = mu_y,
    var_y = var_y
  )
  return(out)
}
