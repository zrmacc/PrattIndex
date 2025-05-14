# Purpose: Pratt index test.
# Updated: 25-05-12

#' Pratt Test
#' 
#' @param g Genotype vector.
#' @param e Environmental exposure.
#' @param y Phenotype.
#' @param standardize Standardize?
#' @return Data.frame.
PrattTest <- function(g, e, y, standardize = TRUE) {
  
  # Standardize.
  if (standardize) {
    g <- Standardize(g)
    e <- Standardize(e)
  }
  
  # Form interaction.
  h <- g * e
  if (standardize) {
    h <- Standardize(h)
  }
  
  x <- cbind(g, e, h)
  
  # Pratt index.
  pratt <- PrattIndex(y = y, x = x)
  kappa <- pratt$pratt[3]
  
  # P-value calculation.
  lambda <- pratt$eigenval
  lambda <- lambda[abs(lambda) > 1e-12]
  davies_output <- CompQuadForm::davies(q = kappa, lambda = pratt$eigenval)
  if (davies_output$ifault != 0) {
    pval <- CompQuadForm::liu(q = kappa, lambda = pratt$eigenval)
  } else {
    pval <- davies_output$Qq
  }
  
  # Output.
  out <- data.frame(
    pratt = kappa,
    mean = sum(lambda),
    var = 2 * sum(lambda^2),
    pval = pval
  )
  return(out)
}



