## Description

This package performs inference on the Pratt index.

## Installation

    remotes::install_github("zrmacc/PrattIndex", build_vignettes = TRUE)

## Simulate data

Simulates data with a phenotype *Y* that depends on 3 factors: genotype
*G*, environment *E*, and the gene-by-environment interaction
*H* = *G* × *E*. The phenotype is generated from the model:
*Y* = *G**β*<sub>*G*</sub> + *E**β*<sub>*E*</sub> + *H**β*<sub>*H*</sub> + *ϵ*,
where *ϵ* is normally distributed with mean 0 and variance
*σ*<sub>*ϵ*</sub><sup>2</sup>.

    set.seed(101)
    data <- PrattIndex::GenData(
      n = 1e3,
      beta_g = 0.5,
      beta_e = 0.5,
      beta_h = 0.5,
      var_resid = 0.5
    )
    head(data)

    ##   g           e          h           y
    ## 1 0 -1.27235813  0.0000000 -0.60459127
    ## 2 0 -0.66307158  0.0000000  0.07030683
    ## 3 1 -0.29654313 -0.2965431 -0.23086574
    ## 4 1 -0.57396643 -0.5739664  0.32733142
    ## 5 0 -0.02617408  0.0000000  0.11470254
    ## 6 0  0.25972722  0.0000000 -2.18571720

## Pratt test from individual-level data

The function `PrattTest` estimates the Pratt index, calculates the
standard error, and computes the 2-sided p-value against the null
hypothesis *H*<sub>0</sub> : *κ*<sub>*H*</sub> = 0. Here, we expect a
significant result because `beta_h != 0`.

    # Run the Pratt index test for interaction.
    result <- PrattIndex::PrattTest(y = data$y, g = data$g, e = data$e, use_score_test = TRUE)

    # Here we expect a significant result because beta_h != 0.
    show(result)

    ##   term method     kappa       cyh         se    chisq         pval
    ## 1    H  Score 0.1695818 0.4002422 0.01197316 200.6047 1.541216e-45

In the output,

- `term` specifies the coefficient for which the Pratt index is being
  estimated; `kappa` is the estimate and `se` its standard error.

- `cyh` is the estimated covariance between the phenotype *Y* and the
  interaction term *H*. For the score test (`use_score_test = TRUE`),
  this is computed under the null model; for the Wald test, it is
  computed under the full model. When `|cyh|` is small (e.g., less than
  0.02), the p-value from the score test is unreliable and will be set
  to `NA`. This threshold can be modified via `tau_cyh`.

- `chisq = (kappa/se)^2` is the *χ*<sup>2</sup> statistic, and `pval`
  the corresponding 2-sided p-value.

## Pratt test from summary statistics

`PrattTestSS` performs Wald tests of *H*<sub>0</sub> : *κ* = 0 for the
Pratt indices of *G*, *E*, and *H* starting from summary statistics. The
required inputs are:

1.  The joint-model coefficients *β̂*<sub>*G*</sub>, *β̂*<sub>*E*</sub>,
    *β̂*<sub>*H*</sub>.
2.  The mean and variance of the variant *G*.
3.  The mean and variance of the environment *E*.
4.  The marginal variance of *Y*.

<!-- -->

    # 1. Fit joint model and extract coefficients
    fit <- lm(y ~ g + e + I(g * e), data = data)
    bg  <- unname(coef(fit)["g"])
    be  <- unname(coef(fit)["e"])
    bh  <- unname(coef(fit)["I(g * e)"])

    # 2. Summary statistics required by PrattTestSS
    var_y  <- var(data$y)
    mu_g   <- mean(data$g)
    var_g  <- var(data$g)
    mean_e <- mean(data$e)
    var_e  <- var(data$e)

    # 3. Run the test from summary statistics
    result_ss <- PrattIndex::PrattTestSS(
      n = nrow(data),
      bg = bg,
      be = be,
      bh = bh,
      mu_g = mu_g,
      var_g = var_g,
      mean_e = mean_e,
      var_e = var_e,
      var_y = var_y
    )
    show(result_ss)

    ##   term method     kappa       cyh         se   chisq         pval
    ## 1    H  Score 0.1639126 0.3860677 0.01191262 189.326 4.461253e-43

#### Notes

- In order to minimize the number of inputs required, `PrattTestSS`
  makes the assumption that *G* and *E* are independent. Note that the
  test starting from individual-level data does not assume *G* ⟂ *E*.

- Minor allele frequency can be provided in place of the mean and
  variance of *G*. If provided, the mean and variance are calculated as
  `mu_g <- 2 * maf` and `var_g <- 2 * maf * (1 - maf)`.
