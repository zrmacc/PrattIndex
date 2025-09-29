#' Pratt Index
#' 
#' @param y Phenotype.
#' @param g Genotype.
#' @param e Environment.
#' @return List of results.
#' @export 
PrattIndex <- function(y, g, e) {
  
  # Calculate in C++.
  results <- PrattIndexCpp(y = y, g = g, e = e)
  
  # Add labeling.
  rownames(results$beta) <- c("G", "E", "H")
  rownames(results$beta_null) <- c("G", "E")
  
  dimnames(results$cov_xx) <- list(c("G", "E", "H"), c("G", "E", "H"))

  dimnames(results$var_beta) <- list(c("bG", "bE", "bH"), c("bG", "bE", "bH"))
  dimnames(results$var_beta_null) <- list(c("bG", "bE"), c("bG", "bE"))
  
  rownames(results$kappa) <- c("G", "E", "H")
  dimnames(results$grad_kappa) <- list(c("G", "E", "H"), c("dkG", "dkE", "dkH"))
  rownames(results$var_kappa) <- c("kG", "kE", "kH")
  
  # Output.
  return(results)
}


#' Score Test
#' 
#' @param beta_null Estimate of (b_G, b_E) from the model that omits H.
#' @param cov_xx Covariance matrix of (G, E, H).
#' @param kappa_h Pratt index for H.
#' @param n Sample size.
#' @param var_resid_null Estimate of residual variance from the model that omits H.
ScoreTest <- function(
    beta_null, 
    cov_xx,
    kappa_h, 
    n,
    var_resid_null
) {
  
  # Variance of Y under the null.
  var_y_null <- as.numeric(
    t(beta_null) %*% cov_xx[1:2, 1:2] %*% beta_null + var_resid_null
  )
  
  # Schur complement.
  cov_h_given_ge <- as.numeric(
    cov_xx[3, 3] - cov_xx[3, 1:2] %*% solve(cov_xx[1:2, 1:2], cov_xx[1:2, 3])
  )
  
  # Pratt variance.
  c0 <- as.numeric(cov_xx[3, 1:2] %*% beta_null)
  var_beta_h_null <- var_resid_null / (n * cov_h_given_ge)
  var_kappa_null <- (c0 / var_y_null)^2 * var_beta_h_null
    
  # Output.
  out <- data.frame(
    term = "H",
    method = "Score",
    kappa = kappa_h,
    se = sqrt(var_kappa_null)
  )
  return(out)
}


#' Pratt Test
#' 
#' Performs the Pratt test for the interaction between G and E.
#' 
#' @param y Phenotype.
#' @param g Genotype.
#' @param e Environment.
#' @param use_score_test Logical.
#' @return Data.frame of results.
#' @export 
PrattTest <- function(y, g, e, use_score_test = TRUE) {
  
  # Pratt components.
  results <- PrattIndex(y = y, g = g, e = e)
  
  if (use_score_test) {
    out <- ScoreTest(
      beta_null = results$beta_null,
      cov_xx = results$cov_xx,
      kappa_h = results$kappa[3],
      n = results$n,
      var_resid_null = results$var_resid_null
    )
  } else {
    out <- data.frame(
      term = "H",
      method = "Wald",
      kappa = results$kappa[3],
      se = sqrt(results$var_kappa[3])
    )
  }
  
  # Finalize test results.
  out$chisq <- (out$kappa / out$se)^2
  out$pval <- stats::pchisq(q = out$chisq, df = 1, lower.tail = FALSE)
  return(out)
}

