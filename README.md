
## Description

This package performs inference on the Pratt index.

## Installation

``` r
remotes::install_github("zrmacc/PrattIndex", build_vignettes = TRUE)
```

## Simulate data

Simulates data with a phenotype $Y$ that depends on 3 factors: genotype
$G$, environment $E$, and the gene-by-environment interaction
$H = G \times E$. The phenotype is generated from the model: $$
Y = G\beta_{G} + E\beta_{E} + H\beta_{H} + \epsilon,
$$ where $\epsilon$ is normally distributed with mean 0 and variance
$\sigma_{\epsilon}^{2}$.

``` r
set.seed(101)
data <- PrattIndex::GenData(
  n = 1e3,
  beta_g = 0.5,
  beta_e = 0.5,
  beta_h = 0.5,
  var_resid = 0.5
)
head(data)
```

    ##   g           e          h           y
    ## 1 0 -1.27235813  0.0000000 -0.60459127
    ## 2 0 -0.66307158  0.0000000  0.07030683
    ## 3 1 -0.29654313 -0.2965431 -0.23086574
    ## 4 1 -0.57396643 -0.5739664  0.32733142
    ## 5 0 -0.02617408  0.0000000  0.11470254
    ## 6 0  0.25972722  0.0000000 -2.18571720

## Calculate the Pratt Index

The function `PrattTest` estimates the p-value, calculates the standard
error, and calculates the 1-sided p-value against the null hypothesis
$H_{0}: \kappa_{H} = 0$. Here, we expect a significant result because
`beta_h != 0`.

``` r
# Run the Pratt index test for interaction.
result <- PrattIndex::PrattTest(y = data$y, g = data$g, e = data$e, use_score_test = TRUE)

# Here we expect a significant result because beta_h != 0.
show(result)
```

    ##   term method     kappa         se    chisq         pval
    ## 1    H  Score 0.1693793 0.01196251 200.4823 1.639005e-45
