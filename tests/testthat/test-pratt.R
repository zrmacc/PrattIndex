# ------------------------------------------------------------------------------
# PrattIndex, PrattTest, ScoreTest
# ------------------------------------------------------------------------------

test_that("PrattIndex returns expected structure", {
  set.seed(1)
  df <- GenData(n = 200, beta_h = 0, var_resid = 1)
  out <- PrattIndex(y = df$y, g = df$g, e = df$e)
  expect_named(out, c(
    "beta", "beta_null", "cov_xx", "cov_xx_inv", "kappa", "grad_kappa",
    "n", "var_beta", "var_beta_null", "var_kappa", "var_resid",
    "var_resid_null", "var_y"
  ))
  expect_equal(nrow(out$beta), 3L)
  expect_equal(rownames(out$beta), c("G", "E", "H"))
  expect_equal(length(out$kappa), 3L)
  expect_equal(out$n, 200L)
  expect_true(all(dim(out$cov_xx) == c(3L, 3L)))
})

test_that("PrattTest returns data.frame with expected columns", {
  set.seed(2)
  df <- GenData(n = 300, beta_h = 0)
  out_score <- PrattTest(df$y, df$g, df$e, use_score_test = TRUE)
  out_wald <- PrattTest(df$y, df$g, df$e, use_score_test = FALSE)
  expect_s3_class(out_score, "data.frame")
  expect_named(out_score, c("term", "method", "kappa", "se", "chisq", "pval"))
  expect_equal(out_score$term, "H")
  expect_equal(out_score$method, "Score")
  expect_equal(out_wald$method, "Wald")
  expect_true(out_score$se > 0)
  expect_true(out_score$pval >= 0 && out_score$pval <= 1)
})

test_that("PrattTest under null gives kappa_H near zero", {
  set.seed(3)
  df <- GenData(n = 2000, beta_g = 0.05, beta_e = 0.1, beta_h = 0, var_resid = 1)
  out <- PrattTest(df$y, df$g, df$e)
  expect_equal(out$term, "H")
  expect_lt(abs(out$kappa), 0.05)
})

# ------------------------------------------------------------------------------
# PrattTestSS
# ------------------------------------------------------------------------------

test_that("PrattTestSS returns expected structure", {
  out <- PrattTestSS(
    n = 1000, bg = 0.1, be = 0.2, bh = 0,
    var_y = 2, maf = 0.3, mean_e = 0, var_e = 1
  )
  expect_s3_class(out, "data.frame")
  expect_named(out, c("term", "method", "kappa", "se", "chisq", "pval"))
  expect_equal(out$term, "H")
  expect_equal(out$method, "Score")
  expect_equal(nrow(out), 1L)
})

test_that("PrattTestSS with bh=0 gives kappa=0", {
  out <- PrattTestSS(
    n = 500, bg = 0.1, be = 0.2, bh = 0,
    var_y = 2, maf = 0.25, mean_e = 0, var_e = 1
  )
  expect_equal(out$kappa, 0)
  expect_true(out$se > 0)
  expect_equal(out$pval, 1)
})

test_that("PrattTestSS kappa_H matches PrattIndex from same data", {
  set.seed(4)
  df <- GenData(n = 500, beta_g = 0.1, beta_e = 0.2, beta_h = 0.15, var_resid = 1)
  fit <- lm(y ~ g + e + I(g * e), data = df)
  bg <- unname(coef(fit)["g"])
  be <- unname(coef(fit)["e"])
  bh <- unname(coef(fit)["I(g * e)"])
  pi_full <- PrattIndex(df$y, df$g, df$e)
  kappa_h_full <- as.numeric(pi_full$kappa[3])
  out_ss <- PrattTestSS(
    n = nrow(df), bg = bg, be = be, bh = bh,
    var_y = var(df$y), maf = mean(df$g) / 2,
    mean_e = mean(df$e), var_e = var(df$e)
  )
  # Summary-statistic kappa_H matches full-data (sample vs population moments => relax tolerance)
  expect_true(isTRUE(all.equal(out_ss$kappa, kappa_h_full, tolerance = 0.02)))
})

# ------------------------------------------------------------------------------
# PrattIFTest, PrattIFTestSS
# ------------------------------------------------------------------------------

test_that("PrattIFTest returns three rows (G, E, H)", {
  set.seed(5)
  df <- GenData(n = 250, beta_h = 0)
  out <- PrattIFTest(df$y, df$g, df$e)
  expect_equal(nrow(out), 3L)
  expect_equal(out$term, c("G", "E", "H"))
  expect_named(out, c("term", "method", "kappa", "se", "chisq", "pval"))
  expect_true(all(out$se > 0))
  expect_true(all(out$pval >= 0 & out$pval <= 1))
})

test_that("PrattIFTestSS returns expected structure", {
  out <- PrattIFTestSS(
    n = 1000, bg = 0.1, be = 0.2, bh = 0.1,
    var_y = 2, maf = 0.3, mean_e = 0, var_e = 1
  )
  expect_equal(nrow(out), 3L)
  expect_equal(out$term, c("G", "E", "H"))
  expect_equal(out$method, rep("Wald", 3L))
  expect_true(all(out$se > 0))
})

test_that("PrattIFTestSS kappa matches PrattIndex from same data", {
  set.seed(6)
  df <- GenData(n = 400, beta_g = 0.08, beta_e = 0.15, beta_h = 0.1, var_resid = 1)
  fit <- lm(y ~ g + e + I(g * e), data = df)
  bg <- unname(coef(fit)["g"])
  be <- unname(coef(fit)["e"])
  bh <- unname(coef(fit)["I(g * e)"])
  pi_full <- PrattIndex(df$y, df$g, df$e)
  out_ss <- PrattIFTestSS(
    n = nrow(df), bg = bg, be = be, bh = bh,
    var_y = var(df$y), maf = mean(df$g) / 2,
    mean_e = mean(df$e), var_e = var(df$e)
  )
  expect_true(isTRUE(all.equal(as.numeric(pi_full$kappa), out_ss$kappa, tolerance = 0.02)))
})

# ------------------------------------------------------------------------------
# SecPrattTest
# ------------------------------------------------------------------------------

test_that("SecPrattTest returns expected structure", {
  set.seed(7)
  df <- GenData(n = 200, beta_h = 0)
  out <- SecPrattTest(df$y, df$g, df$e)
  expect_named(out, c("term", "method", "kappa", "test_stat", "lambda1", "lambda2", "pval"))
  expect_equal(out$term, "H")
  expect_equal(out$method, "MixChi2")
  expect_true(out$pval >= 0 && out$pval <= 1)
})

# ------------------------------------------------------------------------------
# GenData
# ------------------------------------------------------------------------------

test_that("GenData returns data.frame with g, e, h, y", {
  set.seed(8)
  df <- GenData(n = 100, beta_h = 0, var_resid = 1)
  expect_s3_class(df, "data.frame")
  expect_named(df, c("g", "e", "h", "y"))
  expect_equal(nrow(df), 100L)
  expect_equal(df$h, df$g * df$e)
})

test_that("GenData is reproducible with set.seed", {
  set.seed(9)
  a <- GenData(n = 50, beta_h = 0)
  set.seed(9)
  b <- GenData(n = 50, beta_h = 0)
  expect_equal(a$y, b$y)
  expect_equal(a$g, b$g)
  expect_equal(a$e, b$e)
})

# ------------------------------------------------------------------------------
# FitOLS, MatrixSqrt
# ------------------------------------------------------------------------------

test_that("FitOLS recovers coefficients for simple design", {
  set.seed(10)
  n <- 100
  X <- cbind(1, runif(n), runif(n))
  beta_true <- c(1, 0.5, -0.3)
  y <- X %*% beta_true + rnorm(n, 0, 0.5)
  out <- FitOLS(y, X)
  expect_true(isTRUE(all.equal(as.numeric(out$beta), beta_true, tolerance = 0.25)))
  expect_true(out$resid_var > 0)
  expect_equal(nrow(out$var_beta), 3L)
})

test_that("MatrixSqrt returns square root of sympd matrix", {
  S <- crossprod(matrix(rnorm(9), 3, 3)) + diag(3)
  sqrtS <- MatrixSqrt(S, eps = 1e-8)
  expect_equal(dim(sqrtS), c(3L, 3L))
  expect_true(isTRUE(all.equal(c(S), c(sqrtS %*% sqrtS), tolerance = 1e-7)))
})

test_that("MatrixSqrt errors on non-symmetric input", {
  A <- matrix(1:9, 3, 3)
  expect_error(MatrixSqrt(A, eps = 1e-8), "not symmetric")
})
