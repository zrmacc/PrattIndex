# Purpose: Pratt index test.
# Updated: 25-05-14

#' Pratt Test
#' 
#' @param y Phenotype.
#' @param g Genotype vector.
#' @param e Environmental exposure.
#' @param alpha Alpha for confidence intervals.
#' @return Data.frame.
#' @export
PrattTest <- function(y, g, e, alpha = 0.05) {
  
  # Calculate Pratt Index.
  pratt <- PrattIndex(
    y = y,
    g = g,
    e = e
  )
  
  # Output.
  z <- stats::qnorm(p = 1 - alpha / 2)
  out <- data.frame(
    pratt = pratt$pratt,
    se_pratt = sqrt(pratt$var_pratt)
  )
  out$lower <- out$pratt - z * out$se_pratt
  out$upper <- out$pratt + z * out$se_pratt
  out$pval <- stats::pnorm(out$pratt / out$se_pratt, lower.tail = FALSE)
  return(out)
}
