# Purpose: Simulate data.
# Updated: 2025-03-08


#' Simulate Genotype
#' 
#' @param n Sample size.
#' @param maf Minor allele frequency.
#' @return Standardized genotype.
#' @noRd
.GenGeno <- function(n, maf) {
  g <- stats::rbinom(n = n, size = 2, prob = maf)
  out <- (g - mean(g)) / sd(g)
  return(out)
}


#' Simulate Environment
#' 
#' @param n Sample size.
#' @return Standardized environment.
#' @noRd
.GenEnv <- function(n) {
  e <- stats::rnorm(n = n)
  out <- (e - mean(e)) / sd(e)
  return(out)
}


#' Generate Interaction
#' 
#' @param g Standardized genotype.
#' @param e Standardized environment.
#' @return Standardized interaction.
#' @noRd
.GenInt <- function(g, e) {
  h <- g * e
  out <- (h - mean(g)) / sd(h)
  return(out)
}


#' Generate Covariates
#' 
#' @param n Sample size.
#' @param maf Minor allele frequency.
#' @return n x 3 covariate matrix.
GenCovar <- function(n, maf) {
  
  # Genotype.
  g <- .GenGeno(n = n, maf = maf)
  
  # Environment.
  e <- .GenEnv(n = n)
  
  # Interaction.
  h <- .GenInt(g, e)
  
  # Output.
  out <- data.frame(
    g = g,
    e = e,
    h = h
  )
  return(out)
}


#' Generate Phenotype
#' 
#' @param beta 3 x 1 beta vector for the joint model.
#' @param x n x 3 covariate data.frame. 
#' @return n x 1 phenotype vector.
GenPheno <- function(beta, x) {
  
  # Check coefficients.
  resid_var <- (1 - sum(beta^2))
  if (resid_var < 0) {
    stop(glue::glue("Implied residual variance {sprintf('%.1e', resid_var)} is negative."))
  }
  
  # Linear predictor.
  eta <- as.numeric(data.matrix(x) %*% beta)
  
  # Residuals.
  eps <- stats::rnorm(n) 
  eps <- (eps - mean(eps)) / sd(eps)
  eps <- sqrt(resid_var) * eps
  
  # Phenotype.
  y <- eta + eps
  return(y)
}


#' Generate Data
#' 
#' @param n Sample size. 
#' @param beta_g Genetic beta.
#' @param beta_e Environment beta.
#' @param beta_h Interaction beta.
#' @param maf Genetic minor allele frequency. 
#' @param 
#' @export 
GenData <- function(
  n = 1000,  
  beta_g = 0.1,
  beta_e = 0.1,
  beta_h = 0.0,
  maf = 0.25
) {
  
  # Covariates.
  x <- GenCovar(n = n, maf = maf)
  
  # Phenotype.
  beta <- c(beta_g, beta_e, beta_h)
  y <- GenPheno(beta, x)
  
  # Output.
  out <- x
  out$y <- (y - mean(y)) / sd(y)
  return(out)
}

