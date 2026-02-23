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
#' @param var_y Marginal variance of Y.
ScoreTest <- function(
    beta_null, 
    cov_xx,
    kappa_h, 
    n,
    var_resid_null,
    var_y
) {

  # Schur complement.
  cov_h_given_ge <- as.numeric(
    cov_xx[3, 3] - cov_xx[3, 1:2] %*% solve(cov_xx[1:2, 1:2], cov_xx[1:2, 3])
  )
  
  # Pratt variance.
  c0 <- as.numeric(cov_xx[3, 1:2] %*% beta_null)
  var_beta_h_null <- var_resid_null / (n * cov_h_given_ge)
  var_kappa_null <- (c0 / var_y)^2 * var_beta_h_null
    
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
      var_resid_null = results$var_resid_null,
      var_y = results$var_y
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


# ------------------------------------------------------------------------------


#' Second Order Pratt Test
#' 
#' Performs the Pratt test for the interaction between G and E.
#' 
#' @param y Phenotype.
#' @param g Genotype.
#' @param e Environment.
#' @return Data.frame of results.
#' @export 
SecPrattTest <- function(y, g, e) {
  
  # Pratt components.
  results <- PrattIndex(y = y, g = g, e = e)

  # Test statistic.
  n <- length(y)
  kappa_h <- results$kappa[3]
  test_stat <- n * kappa_h
  
  # Calculate eigenvalues.
  cov_xx <- results$cov_xx
  cov_xx_inv <- solve(cov_xx)
  lambda1 <- 0.5 * (1 + sqrt(cov_xx[3, 3] * cov_xx_inv[3, 3]))
  lambda2 <- 0.5 * (1 - sqrt(cov_xx[3, 3] * cov_xx_inv[3, 3]))
  
  # Calculate p-value.
  pval <- CompQuadForm::davies(
    q = test_stat,
    lambda = c(lambda1, lambda2)
  )$Qq
  
  # Output.
  out <- data.frame(
    term = "H",
    method = "MixChi2",
    kappa = kappa_h,
    test_stat = test_stat,
    lambda1 = lambda1,
    lambda2 = lambda2,
    pval = pval
  )
  return(out)
}


# ------------------------------------------------------------------------------


#' Pratt Influence Function Test
#'
#' Computes Wald tests for the Pratt indices corresponding to G, E, and
#' H = G * E with inference via the empirical influence function.
#'
#' @param y (n x 1) numeric vector. Phenotype.
#' @param g (n x 1) numeric vector. Genotype.
#' @param e (n x 1) numeric vector. Environment.
#' @return Data.frame of results.
#' @export
PrattIFTest <- function(y, g, e) {
  
  # Pratt components.
  results <- PrattIndex(y = y, g = g, e = e)
  
  # Influence function.
  psi <- PrattInfluenceCpp(y = y, g = g, e = e)

  # Variance estimate.
  n <- nrow(psi)
  var_kappa_if <- colSums(psi^2) / (n^2)
  se <- sqrt(var_kappa_if)
  
  # Results.
  out <- data.frame(
    term = c("G", "E", "H"),
    method = "Wald",
    kappa = as.numeric(results$kappa),
    se = as.numeric(se),
    stringsAsFactors = FALSE
  )
  
  # Calculate p-values.
  out$chisq <- (out$kappa / out$se)^2
  out$pval <- stats::pchisq(q = out$chisq, df = 1, lower.tail = FALSE)
  return(out)
}


# ------------------------------------------------------------------------------


#' Pratt influence-function-style test from summary statistics
#'
#' Computes Wald tests for the Pratt indices for G, E, and H using the same
#' inputs as \code{PrattTestSS}. The standard errors use the delta-method
#' variance (gradient of kappa with respect to beta times variance of beta),
#' which is asymptotically equivalent to the empirical influence-function
#' variance used by \code{PrattIFTest} when the model holds. The exact
#' influence-function variance (mean of IF^2) would require individual-level
#' data or an analytic expectation under the full distribution of (G, E, Y).
#'
#' For \eqn{\kappa_H}{kappa_H}, the point estimate agrees with
#' \code{PrattTestSS}, but the standard error (and hence test statistic and
#' p-value) is not numerically equal: \code{PrattTestSS} uses the
#' null-calibrated (score) variance under \eqn{H_0: \beta_H = 0}{H0: beta_H = 0},
#' while \code{PrattIFTestSS} uses the full-model delta-method variance. They
#' are asymptotically equivalent under the null.
#'
#' @param n Sample size.
#' @param bg Coefficient of G from the joint model Y ~ G + E + H.
#' @param be Coefficient of E from the joint model Y ~ G + E + H.
#' @param bh Coefficient of H (interaction) from the joint model Y ~ G + E + H.
#' @param var_y Marginal variance of the phenotype Y.
#' @param maf Minor allele frequency of G (additive coding 0, 1, 2 assumed).
#' @param mean_e Marginal mean of E.
#' @param var_e Marginal variance of E.
#' @return Data.frame with \code{term} (G, E, H), \code{method} (\dQuote{Wald}),
#'   \code{kappa}, \code{se}, \code{chisq}, and \code{pval}.
#' @export
PrattIFTestSS <- function(n, bg, be, bh, var_y, maf, mean_e, var_e) {

  mu_g <- 2 * maf
  sigma_gg <- 2 * maf * (1 - maf)
  sigma_ge <- 0
  sigma_ee <- var_e
  mu_e_val <- mean_e
  sigma_gh <- sigma_gg * mu_e_val
  sigma_eh <- sigma_ee * mu_g
  sigma_hh <- sigma_gg * sigma_ee + sigma_ee * (mu_g^2) + sigma_gg * (mu_e_val^2)

  sigma_xx <- matrix(c(
    sigma_gg, sigma_ge, sigma_gh,
    sigma_ge, sigma_ee, sigma_eh,
    sigma_gh, sigma_eh, sigma_hh
  ), 3L, 3L)
  beta <- c(bg, be, bh)
  cov_yx <- as.numeric(sigma_xx %*% beta)

  # Pratt indices for G, E, H
  kappa <- (beta * cov_yx) / var_y

  # Residual variance from joint model: var(Y) = Var(X'beta) + var_resid
  var_resid <- var_y - as.numeric(t(beta) %*% sigma_xx %*% beta)
  var_resid <- max(0, var_resid)

  # var(beta_hat) = var_resid * (X'X)^{-1}; with centered X, (1/n) X'X ~ SigmaXX
  sigma_xx_inv <- solve(sigma_xx)
  var_beta <- (var_resid / n) * sigma_xx_inv

  # Gradient of kappa w.r.t. beta: (d/d beta_k) kappa_j = (1/var_y)[ cov_yx(j)*1(k=j) + beta(j)*Sigma_XX[j,k] ]
  grad_kappa <- matrix(0, 3L, 3L)
  for (j in 1:3) {
    term <- rep(0, 3)
    term[j] <- cov_yx[j]
    term <- term + beta[j] * sigma_xx[j, ]
    grad_kappa[, j] <- term / var_y
  }

  # Delta-method variance of kappa
  var_kappa <- numeric(3L)
  for (j in 1:3) {
    var_kappa[j] <- as.numeric(t(grad_kappa[, j]) %*% var_beta %*% grad_kappa[, j])
  }
  se <- sqrt(var_kappa)

  out <- data.frame(
    term = c("G", "E", "H"),
    method = "Wald",
    kappa = as.numeric(kappa),
    se = se,
    stringsAsFactors = FALSE
  )
  out$chisq <- (out$kappa / out$se)^2
  out$pval <- stats::pchisq(q = out$chisq, df = 1L, lower.tail = FALSE)
  return(out)
}


#' Null-calibrated Pratt index test from summary statistics
#'
#' Performs the null-calibrated (score) test for the interaction Pratt index
#' using only summary statistics, under the assumption G is independent of E.
#' Uses the joint model coefficients (G, E, H = G*E), marginal variance of Y,
#' MAF, and the mean and variance of E to reconstruct the null variance.
#'
#' For \eqn{\kappa_H}{kappa_H}, the point estimate agrees with
#' \code{PrattIFTestSS}, but the standard error (and hence test statistic and
#' p-value) is not numerically equal: this function uses the null-calibrated
#' (score) variance under \eqn{H_0: \beta_H = 0}{H0: beta_H = 0}, while
#' \code{PrattIFTestSS} uses the full-model delta-method variance. They are
#' asymptotically equivalent under the null.
#'
#' @param n Sample size.
#' @param bg Coefficient of G from the joint model Y ~ G + E + H.
#' @param be Coefficient of E from the joint model Y ~ G + E + H.
#' @param bh Coefficient of H (interaction) from the joint model Y ~ G + E + H.
#' @param var_y Marginal variance of the phenotype Y.
#' @param maf Minor allele frequency of G (additive coding 0, 1, 2 assumed).
#' @param mean_e Marginal mean of E.
#' @param var_e Marginal variance of E.
#' @return Data.frame with \code{term}, \code{method}, \code{kappa} (Pratt index
#'   for H), \code{se}, \code{chisq}, and \code{pval}.
#' @export
PrattTestSS <- function(n, bg, be, bh, var_y, maf, mean_e, var_e) {

  # G: additive coding 0,1,2 => E[G] = 2*maf, Var(G) = 2*maf*(1-maf)
  mu_g <- 2 * maf
  sigma_gg <- 2 * maf * (1 - maf)
  sigma_ge <- 0

  sigma_ee <- var_e
  mu_e <- mean_e

  # Under G ⊥ E: Σ_GH = Σ_GG*μ_E, Σ_EH = Σ_EE*μ_G,
  # Σ_HH = Σ_GG*Σ_EE + Σ_EE*μ_G² + Σ_GG*μ_E²
  sigma_gh <- sigma_gg * mu_e
  sigma_eh <- sigma_ee * mu_g
  sigma_hh <- sigma_gg * sigma_ee + sigma_ee * (mu_g^2) + sigma_gg * (mu_e^2)

  # Pratt index for H: κ_H = β_H (Σ_HG*β_G + Σ_HE*β_E + Σ_HH*β_H) / σ²_Y
  numer_kappa <- bh * (sigma_gh * bg + sigma_eh * be + sigma_hh * bh)
  kappa_h <- numer_kappa / var_y

  # Null (reduced) model coefficients: under independence
  # Σ_GY = Σ_GG(β_G + μ_E*β_H), Σ_EY = Σ_EE(β_E + μ_G*β_H)
  # (β_G,0, β_E,0)' = Σ^{-1}_{[G,E]} Σ_{[G,E],Y} => β_G,0 = β_G + μ_E*β_H, β_E,0 = β_E + μ_G*β_H
  bg_null <- bg + mu_e * bh
  be_null <- be + mu_g * bh

  sigma_gy <- sigma_gg * bg_null
  sigma_ey <- sigma_ee * be_null

  # Residual variance under null: σ²_ε,0 = σ²_Y - Σ'_{Y,[G,E]} Σ^{-1}_{[G,E]} Σ_{[G,E],Y}
  # With diagonal Σ_{[G,E]}: = σ²_Y - (Σ_GY²/Σ_GG + Σ_EY²/Σ_EE)
  var_resid_null <- var_y - (sigma_gy^2 / sigma_gg + sigma_ey^2 / sigma_ee)
  var_resid_null <- max(0, var_resid_null)

  # σ²_H|GE = Σ_HH - Σ_{H,[G,E]} Σ^{-1}_{[G,E]} Σ_{[G,E],H}
  sigma_h_given_ge <- sigma_hh - (sigma_gh^2 / sigma_gg + sigma_eh^2 / sigma_ee)
  sigma_h_given_ge <- max(.Machine$double.eps, sigma_h_given_ge)

  # Null variance of κ̂_H: V0(κ̂_H) = (Σ_HG*β_G,0 + Σ_HE*β_E,0)²/σ⁴_Y * σ²_ε,0/(N*σ²_H|GE)
  c0 <- sigma_gh * bg_null + sigma_eh * be_null
  var_kappa_null <- (c0 / var_y)^2 * (var_resid_null / (n * sigma_h_given_ge))
  se <- sqrt(var_kappa_null)

  chisq <- (kappa_h / se)^2
  pval <- stats::pchisq(q = chisq, df = 1, lower.tail = FALSE)

  out <- data.frame(
    term = "H",
    method = "Score",
    kappa = kappa_h,
    se = se,
    chisq = chisq,
    pval = pval,
    stringsAsFactors = FALSE
  )
  return(out)
}

