library(testthat)

test_that("Test Pratt index calculation.", {
  
  # Compare.
  Compare <- function(x, y) {
    x <- unname(as.numeric(x))
    y <- unname(as.numeric(y))
    expect_equal(x, y, tolerance = 1e-3, ignore_attr = TRUE)
  }
  
  # Manual.
  withr::local_seed(1010)
  n <- 1000
  df <- PrattIndex::GenData(n = n)
  y <- df$y
  x <- as.matrix(df[, c("g", "e", "h")])
  
  # Expected.
  fit <- lm(y ~ 0 + x)
  beta <- coef(fit)
  
  r <- t(x) %*% y / n
  pratt <- r * beta
  
  # Observed.
  # Note: turning on standardization results in small numeric discrepancies.
  obs <- PrattIndex(y = y, x = x, standardize = FALSE)
  
  # Compare.
  Compare(obs$beta, beta)
  Compare(obs$r, r)
  Compare(obs$pratt, pratt)
  expect_equal(sum(pratt), obs$var_exp)
  
})
