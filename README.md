
<!-- README.md is generated from README.Rmd. Please edit that file -->

# stIsing

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![CRAN
status](https://www.r-pkg.org/badges/version/stIsing)](https://CRAN.R-project.org/package=stIsing)
<!-- badges: end -->

The stIsing package allows to the estimation of binary graph models

## Installation

You can install the development version of stIsing like so:

``` r
devtools::install_github("giuseppealfonzetti/stIsing")
```

## Example

Set-up the structure of the data generating graph:

``` r
library(stIsing)

p <- 10
d <- p + p*(p-1)/2

parNodes <- rep(c(-1,1), p/2)
parEdgesMat <- matrix(0,p, p)
counter <- 1
for (col in 1:(p-1)) {
    for (row in (col+1):p) {
        if(abs(col-row)==1 & (min(col,row)!=p/2)) parEdgesMat[row, col] <- .5
        if(abs(col-row)==p/2) parEdgesMat[row, col] <- -.5
    }
}
parEdgesMat
#>       [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10]
#>  [1,]  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0     0
#>  [2,]  0.5  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0     0
#>  [3,]  0.0  0.5  0.0  0.0  0.0  0.0  0.0  0.0  0.0     0
#>  [4,]  0.0  0.0  0.5  0.0  0.0  0.0  0.0  0.0  0.0     0
#>  [5,]  0.0  0.0  0.0  0.5  0.0  0.0  0.0  0.0  0.0     0
#>  [6,] -0.5  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0     0
#>  [7,]  0.0 -0.5  0.0  0.0  0.0  0.5  0.0  0.0  0.0     0
#>  [8,]  0.0  0.0 -0.5  0.0  0.0  0.0  0.5  0.0  0.0     0
#>  [9,]  0.0  0.0  0.0 -0.5  0.0  0.0  0.0  0.5  0.0     0
#> [10,]  0.0  0.0  0.0  0.0 -0.5  0.0  0.0  0.0  0.5     0

true_theta <- parNodes
for (col in 1:(p-1)) {
    for (row in (col+1):p) {
        true_theta <- c(true_theta, parEdgesMat[row, col])
    }
}
true_theta
#>  [1] -1.0  1.0 -1.0  1.0 -1.0  1.0 -1.0  1.0 -1.0  1.0  0.5  0.0  0.0  0.0 -0.5
#> [16]  0.0  0.0  0.0  0.0  0.5  0.0  0.0  0.0 -0.5  0.0  0.0  0.0  0.5  0.0  0.0
#> [31]  0.0 -0.5  0.0  0.0  0.5  0.0  0.0  0.0 -0.5  0.0  0.0  0.0  0.0  0.0 -0.5
#> [46]  0.5  0.0  0.0  0.0  0.5  0.0  0.0  0.5  0.0  0.5
length(true_theta)
#> [1] 55
```

Generate data from the IsingSampler package

``` r
n <- 2000
graph_mat <- ising_from_theta_to_emat(true_theta, p)
graph_mat <- graph_mat + t(graph_mat)
thr_vec <- true_theta[1:p]

seed <- 123
set.seed(seed)
data <- IsingSampler::IsingSampler(n, graph_mat, thr_vec, 1, method = "direct")
```

Fit the model:

``` r
Q <- rep(TRUE, length(true_theta))                                # No constraints
theta_init <- rep(0, length(true_theta))                          # Initialisation

fit <- fit_isingGraph(
    DATA_LIST = list(DATA = as.matrix(data), CONSTRAINTS = Q),
    METHOD = 'ucminf',
    CPP_CONTROL = list(),
    INIT = theta_init,
    ITERATIONS_SUBSET = NULL,
    VERBOSEFLAG = 0
)
#> 1. Initialising at init vector.
#> 2. Optimising with ucminf...
#> 3. Done! (2.54 secs)
mean((fit$theta-true_theta)^2)
#> [1] 0.02966388
mean((theta_init-true_theta)^2)
#> [1] 0.2409091
```
