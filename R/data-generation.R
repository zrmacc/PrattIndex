# Purpose: Simulate data.
# Updated: 2025-05-14


#' Simulate Genotype
#' 
#' @param n Sample size.
#' @param maf Minor allele frequency.
#' @return Genotype.
#' @noRd
.GenGeno <- function(n, maf) {
  g <- stats::rbinom(n = n, size = 2, prob = maf)
  return(g)
}


#' Simulate Environment
#' 
#' @param n Sample size.
#' @return Environment.
#' @noRd
.GenEnv <- function(n) {
  e <- stats::rnorm(n = n)
  return(e)
}


#' Generate Interaction
#' 
#' @param g Genotype.
#' @param e Environment.
#' @return Interaction.
#' @noRd
.GenInt <- function(g, e) {
  h <- g * e
  return(h)
}


#' Standardize
#' 
#' @param x Numeric vector.
#' @return Standardized vector.
#' @noRd
Standardize <- function(x) {
  out <- x - mean(x)
  out <- out / sqrt(mean(out^2))
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
#' Uses the centered g and e to form the interaction.
#' 
#' @param beta 3 x 1 beta vector for the joint model.
#' @param x n x 3 covariate data.frame. 
#' @param var_exp Proportion of variance explained by main effects.
#' @param var_resid Residual variance explained. Set to null to specify `var_exp` instead.
#' @return n x 1 phenotype vector.
GenPheno <- function(
    beta, 
    x,
    var_exp,
    var_resid = NULL
) {
  
  if (is.null(var_resid)) {
    # Generate covariance matrix.
    S <- stats::cov(x)
    
    # Calculate variance of main effect.
    main_var <- as.numeric(t(beta) %*% S %*% beta)
    
    # Calculate residual variance.
    sigma2 <- (main_var - var_exp * main_var) / var_exp
    
  } else {
    
    sigma2 <- var_resid
    
  }

  # Linear predictor.
  eta <- as.numeric(data.matrix(x) %*% beta)
  
  # Residuals.
  n <- length(eta)
  eps <- sqrt(sigma2) * stats::rnorm(n) 
  
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
#' @param var_exp Proportion of variance explained by main effects.
#' @param var_resid Residual variance. Set to null to specify `var_exp` instead.
#' @export 
GenData <- function(
  n = 1000,  
  beta_g = 0.1,
  beta_e = 0.1,
  beta_h = 0.0,
  maf = 0.25,
  var_exp = 0.25,
  var_resid = NULL
) {
  
  # Covariates.
  x <- GenCovar(n = n, maf = maf)

  # Phenotype.
  beta <- c(beta_g, beta_e, beta_h)
  y <- GenPheno(
    beta = beta, 
    x = x,
    var_exp = var_exp,
    var_resid = var_resid
  )
  
  # Output.
  out <- x
  out$y <- y
  return(out)
}

