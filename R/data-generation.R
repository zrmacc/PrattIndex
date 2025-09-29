# Purpose: Simulate data.
# Updated: 2025-09-27


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
    assertthat::assert_that(
      !is.null(maf),
      msg = "If type_g is binom, specify maf."
    )
    g <- stats::rbinom(n = n, size = 2, prob = maf)
  } else if (type_g == "normal") {
    assertthat::assert_that(
      !is.null(mu_g) && !is.null(var_g),
      msg = "If type_g is normal, specify mu_g and var_g."
    )
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
#' @param n_e Number of trials for E, if type_e = "binom".
#' @param p_e Probability of success for E, if type_e = "binom".
#' @param type_e Either "binom" or "normal".
#' @param var_e Variance of E.
#' @return Environment.
#' @noRd
.GenEnv <- function(
    n, 
    mu_e = 0,
    n_e = NULL,
    p_e = NULL,
    type_e = "normal",
    var_e = 1
) {

  if (type_e == "binom") {
    assertthat::assert_that(
      !is.null(n_e) && !is.null(p_e),
      msg = "If type_e is binom, specify n_e and p_e."
    )
    e <- stats::rbinom(n = n, size = n_e, prob = p_e)
  } else if (type_e == "normal") {
    assertthat::assert_that(
      !is.null(mu_e) && !is.null(var_e),
      msg = "If type_e is normal, specify mu_e and var_e."
    )
    e <- stats::rnorm(n = n, mean = mu_e, sd = sqrt(var_e))
  } else {
    stop(glue::glue("E of type {type_e} not implemented."))
  }
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
#' @param maf Minor allele frequency for G, if type_g = "binom".
#' @param mu_e Mean of E, if type_e = "normal".
#' @param mu_g Mean of G, if type_g = "normal".
#' @param n_e Number of trials for E, if type_e = "binom".
#' @param p_e Probability of success for E, if type_e = "binom".
#' @param type_e Either "binom" or "normal".
#' @param type_g Either "binom" or "normal".
#' @param var_e Variance of E, if type_e = "normal".
#' @param var_g Variance of G, if type_g = "normal".
#' @return n x 3 covariate matrix.
GenCovar <- function(
    n, 
    maf = 0.25,
    mu_e = 0,
    mu_g = NULL,
    n_e = NULL,
    p_e = NULL,
    type_e = "normal",
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
    n_e = n_e,
    p_e = p_e,
    type_e = type_e,
    var_e = var_e
  )
  
  # Interaction.
  h <- .GenInt(g, e)
  
  # Output.
  out <- data.frame(g = g, e = e, h = h)
  return(out)
}


#' Generate Phenotype
#' 
#' @param beta 3 x 1 beta vector for the joint model.
#' @param x n x 3 covariate data.frame. 
#' @param var_exp Proportion of variance explained by (G, E, H).
#' @param var_resid Residual variance Set to null to specify `var_exp` instead.
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
#' Simulates a data set containing the phenotype, genotype, environment, and 
#' interaction. Genotype can be either binomial or normal. If binomial, specify the
#' minor allele frequency. If normal, specify mu_g and var_g. Likewise, environment
#' can be binomial or normal. If binomial, specify n_e and p_e. If normal, specify
#' mu_e and var_e. The residual phenotypic variance can either be specified directly,
#' or it can be inferred by specifying the variance explained by (G, E, H).
#' 
#' @param beta_g Genetic beta.
#' @param beta_e Environment beta.
#' @param beta_h Interaction beta.
#' @param maf Minor allele frequency for G, if type_g = "binom". 
#' @param mu_e Mean of E, if type_e = "normal".
#' @param mu_g Mean of G, if type_g = "normal".
#' @param n Sample size. 
#' @param n_e Number of trials for E, if type_e = "binom".
#' @param p_e Probability of success for E, if type_e = "binom".
#' @param type_e Either "binom" or "normal".
#' @param type_g Either "binom" or "normal".
#' @param var_e Variance of E, if type_e = "normal".
#' @param var_exp Proportion of variance in Y explained by (G, E, H).
#' @param var_g Variance of G, if type_g = "normal".
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
  n_e = NULL,
  p_e = NULL,
  type_e = "normal",
  type_g = "binom",
  var_e = 1,
  var_exp = 0.25,
  var_g = NULL,
  var_resid = NULL
) {
  
  # Covariates.
  x <- GenCovar(
    n = n,
    maf = maf,
    mu_e = mu_e,
    mu_g = mu_g,
    n_e = n_e,
    p_e = p_e,
    type_e = type_e,
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

