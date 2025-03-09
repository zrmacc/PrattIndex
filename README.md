
## Description

This package performs inference on the Pratt index.

## Installation

``` r
remotes::install_github("zrmacc/PrattIndex", build_vignettes = TRUE)
```

## Simulate data

Simulates data with a phenotype $y$ and 3 covariates $g$ denoting
genotype, $e$ denoting the environment, and $h = g \times e$ denoting
the interaction. All variables are standardized to have mean 0 and
variance 1. The phenotype is related to the covariates by: $$
y = g\beta_{G} + e\beta_{E} + h\beta_{H} + \epsilon,
$$ where $\epsilon$ is normally distributed with mean 0 and variance
$1 - \beta_{G}^{2} - \beta_{E}^{2} - \beta_{H}^{2}$.

``` r
data <- PrattIndex::GenData(
  n = 1e3,
  beta_g = 0.1,
  beta_e = 0.1,
  beta_h = 0.0
)
head(data)
```

    ##            g          e           h          y
    ## 1  0.8303156 -0.6290048 -0.52599458 -0.3364087
    ## 2 -0.7913946 -0.1200337  0.09567101 -0.4974543
    ## 3 -0.7913946  1.5456556 -1.23194107  3.0596695
    ## 4 -0.7913946  1.2509762 -0.99707138  1.2758921
    ## 5  0.8303156  0.5429623  0.45404307  1.7484262
    ## 6  0.8303156  2.6401242  2.20775917  2.5729902

## Calculate the Pratt Index

The Pratt index is the element-wise product of the marginal correlations
$\hat{r} = X^{\top}y$ with the regression coefficients from the joint
model $\hat{\beta} + (X^{\top}X)^{-1}X^{\top}y$. `PrattIndex` calculates
the Pratt index as well as various intermediates used in the
calculation.

``` r
y <- data$y
x <- as.matrix(data[, c("g", "e", "h")])
pratt_comps <- PrattIndex::PrattIndex(y = y, x = x)

cat("Joint model regression coefficients:")
show(pratt_comps$beta)
cat("\n")

cat("Marginal correlation coefficients:")
show(pratt_comps$r)
cat("\n")

cat("Pratt indices:")
show(pratt_comps$pratt)
cat("\n")
```

    ## Joint model regression coefficients:            [,1]
    ## [1,]  0.06568554
    ## [2,]  0.12561034
    ## [3,] -0.04462223
    ## 
    ## Marginal correlation coefficients:            [,1]
    ## [1,]  0.06776428
    ## [2,]  0.12749280
    ## [3,] -0.05589877
    ## 
    ## Pratt indices:            [,1]
    ## [1,] 0.004451133
    ## [2,] 0.016014414
    ## [3,] 0.002494328
