
# GHSGEM: Generalized expectation-maximization (GEM) algorithm with graphical horseshoe (GHS) prior for network estimation

This R package implements a generalized expectation-maximization
algorithm with agraphical horseshoe prior for network estimation and
covariance and precision (inverse of the covariance) matrices
estimation.

## Installation

The GHSGEM package can be installed using the following code:

``` r
if(!require("devtools", quietly = TRUE)) {
  install.packages("devtools")
}
devtools::install_github("THautamaki/GHSGEM")
```

## Example of usage

This is a minimum working example.

``` r
# Install R package huge if not yet installed.
if(!require("huge", quietly = TRUE)) {
  install.packages("huge")
}

library(GHSGEM)
library(huge)

n <- 200  # number of observations (sample size)
p <- 100  # number of variables

# Generate simulated data using huge.generator.
sim <- huge.generator(n = n, d = p, graph = "scale-free")
```

    ## Generating data from the multivariate normal distribution with the scale-free graph structure....done.

``` r
# Run GHS GEM algorithm.
map <- GHS_MAP_estimation(sim$data, verbose = 0)
```

    ## Total iterations: 80. Elapsed time: 1.39644 s. Final difference: 9.82423e-05

``` r
# Calculate and print confusion matrix.
(cm <- conf_matrix(sim$theta, map$Theta_est))
```

    ##        Estim. P Estim. N
    ## True P       61       38
    ## True N       14     4837

``` r
# Calculate and print some performace scores.
round(calculate_scores(cm)[, c("MCC", "F1", "TPR", "FDR")], 4)
```

    ##      MCC     F1    TPR    FDR
    ## 1 0.7029 0.7011 0.6162 0.1867
