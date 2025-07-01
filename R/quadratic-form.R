# Purpose: Functions related to quadratic forms.
# Updated: 2025-06-15

#' Calculate P-value for Quadratic Form
#' 
#' Calculate the p-value for a statistic \deqn{T = y^{\top}Ay}, where
#' \deqn{y ~ N(0, v)} and A is any symmetric positive-definite matrix.
#'
#' @param a Central matrix of the quadratic form.
#' @param test_stat Test statistic.
#' @param v Variance of the random vector.
#' @param eps Minimum eigenvalue magnitude.
#' @param method "davies" or "liu". If "davies", "liu" will provide a fallback.
#' @return Numeric p-value.
#' @export
CalcQuadFormP <- function(
    a, 
    test_stat, 
    v,
    eps = 1e-8,
    method = "davies"
) {
  
  # Calculate eigenvalues.
  s <- a %*% v
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

