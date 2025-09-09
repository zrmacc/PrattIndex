# Purpose: Simulate data.
# Updated: 2025-09-09


#' Simulate Genotype
#' 
#' @param n Sample size.
#' @param maf Minor allele frequency.
#' @param mu_g Mean of G. Only used if type_g = "normal".
#' @param type_g Either "binom" or "normal".
#' @param var_g Variance of G. Only used if type_g = "normal".
#' @return Genotype.
#' @noRd
.GenGeno <- function(
    n, 
    maf = 0.25,
    mu_g = NULL,
    type_g = "binom",
    var_g = NULL
  ) {
  
  if (type_g == "binom") {
    g <- stats::rbinom(n = n, size = 2, prob = maf)
  } else if (type_g == "normal") {
    g <- stats::rnorm(n = n, mean = mu_g, sd = sqrt(var_g))
  } else {
    stop(glue::glue("G of type {type_g} not implemented."))
  }
  
  return(g)
}


#' Simulate Environment
#' 
#' @param n Sample size.
#' @param mu_e Mean of E.
#' @param var_e Variance of E.
#' @return Environment.
#' @noRd
.GenEnv <- function(n, mu_e = 0, var_e = 1) {
  e <- stats::rnorm(n = n, mean = mu_e, sd = sqrt(var_e))
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


#' Generate Covariates
#' 
#' @param n Sample size.
#' @param maf Minor allele frequency.
#' @param mu_e Mean of E.
#' @param mu_g Mean of G, only used if type_g = "normal".
#' @param type_g Either "binom" or "normal".
#' @param var_e Variance of E.
#' @param var_g Variance of G, only used if type_g = "normal".
#' @return n x 3 covariate matrix.
GenCovar <- function(
    n, 
    maf = 0.25,
    mu_e = 0,
    mu_g = NULL,
    type_g = "binom",
    var_e = 1,
    var_g = NULL
) {
  
  # Genotype.
  g <- .GenGeno(
    n = n, 
    maf = maf,
    mu_g = mu_g,
    type_g = type_g,
    var_g = var_g
  )
  
  # Environment.
  e <- .GenEnv(
    n = n,
    mu_e = mu_e,
    var_e = var_e
  )
  
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
#' @param beta_g Genetic beta.
#' @param beta_e Environment beta.
#' @param beta_h Interaction beta.
#' @param maf Genetic minor allele frequency. 
#' @param mu_e Mean of E.
#' @param mu_g Mean of G, only used if type_g = "normal".
#' @param n Sample size. 
#' @param type_g Either "binom" or "normal".
#' @param var_e Variance of E.
#' @param var_exp Proportion of variance in Y explained by main effects.
#' @param var_g Variance of G, only used if type_g = "normal".
#' @param var_resid Residual variance. Set to null to specify `var_exp` instead.
#' @export 
GenData <- function(
  beta_g = 0.1,
  beta_e = 0.1,
  beta_h = 0.0,
  maf = 0.25,
  mu_g = NULL,
  mu_e = 0,
  n = 1000,  
  type_g = "binom",
  var_e = 1,
  var_exp = 0.25,
  var_g = NULL,
  var_resid = NULL
) {
  
  # Covariates.
  x <- GenCovar(
    maf = maf,
    n = n, 
    mu_e = mu_e,
    mu_g = mu_g,
    type_g = type_g,
    var_e = var_e,
    var_g = var_g
  )

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

