# Purpose: Functions related to quadratic forms.
# Updated: 2025-09-09

#' Calculate P-value for Quadratic Form
#' 
#' Calculate the p-value for a statistic \deqn{T = z^{\top}Az}, where
#' \deqn{z ~ N(0, cov_z)} and A is any symmetric positive-definite matrix.
#'
#' @param cov_z Covariance matrix of the random vector.
#' @param kernel Central matrix of the quadratic form.
#' @param test_stat Test statistic.
#' @param eps Minimum eigenvalue magnitude.
#' @param method "davies" or "liu". If "davies", "liu" will provide a fallback.
#' @return Numeric p-value.
#' @export
CalcQuadFormP <- function(
    cov_z,
    kernel, 
    test_stat, 
    eps = 1e-8,
    method = "davies"
) {
  
  # Ensure symmetric.
  cov_z <- 0.5 * (cov_z + t(cov_z))
  kernel <- 0.5 * (kernel + t(kernel))
  
  # Calculate eigenvalues.
  s <- kernel %*% cov_z
  lambda <- eigen(s, symmetric = FALSE, only.values = TRUE)$values
  lambda <- lambda[abs(lambda) > eps]
  
  # Calculate p-value
  if (method == "liu") {
    pval <- CompQuadForm::liu(q = test_stat, lambda = lambda)
    method <- "liu"
  } else {
    
    # Try Davies.
    davies_out <- tryCatch(
      {CompQuadForm::davies(q = test_stat, lambda = lambda)},
      error = function(e) {return(list(ifault = 1))}
    )
    
    if (davies_out$ifault != 0) {
      pval <- CompQuadForm::liu(q = test_stat, lambda = lambda)
      method <- "liu"
    } else {
      pval <- davies_out$Qq
      method <- "davies"
    }
    
  }
  
  # Output.
  out <- data.frame(
    stat = test_stat,
    pval = pval,
    method = method
  )
  return(out)
}

