library(testthat)

test_that("Test Pratt index calculation.", {
  
  # Compare.
  Compare <- function(x, y, tol = 1e-3) {
    x <- unname(as.numeric(x))
    y <- unname(as.numeric(y))
    expect_equal(x, y, tolerance = tol, ignore_attr = TRUE)
  }
  
  # Manual.
  withr::local_seed(101)
  n <- 10000
  df <- PrattIndex::GenData(
    n = n,
    beta_g = 0.5,
    beta_e = 0.5,
    beta_h = 0.5,
    var_resid = 1
  )
  y <- df$y
  g <- df$g
  e <- df$e
  x <- as.matrix(df[, c("g", "e", "h")])
  
  # Expected.
  fit <- lm(y ~ x)
  beta <- coef(fit)
  sigma2 <- sigma(fit)^2
  
  # Pratt.
  pratt <- beta[2:4] * cor(y, x)
  pratt_h <- pratt[3]
  
  # Variance explained.
  var_exp <- sum(pratt)
  
  # Observed.
  # Note: turning on standardization results in small numeric discrepancies.
  obs <- PrattIndex(y = y, g = g, e = e)
  
  # Compare.
  Compare(obs$beta, beta)
  Compare(obs$sigma2, sigma2)
  Compare(obs$pratt, pratt_h, tol = 5e-2)
})
