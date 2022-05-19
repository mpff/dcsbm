
<!-- README.md is generated from README.Rmd. Please edit that file -->

# dcsbm

<!-- badges: start -->

[![R-CMD-check](https://github.com/mpff/dcsbm/workflows/R-CMD-check/badge.svg)](https://github.com/mpff/dcsbm/actions)
<!-- badges: end -->

The goal of dcsbm is to provide methods for estimating a two-way degree
corrected stochastic block model for directed, weighted graphs. Uses the
‘igraph’ library <https://igraph.org> for graph handling. See Peixoto
(2014) \<10.1103/PhysRevE.89.012804\> for details on the inference
algorithm.

## Installation

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("mpff/dcsbm")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(dcsbm)
#> Loading required package: igraph
#> 
#> Attaching package: 'igraph'
#> The following objects are masked from 'package:stats':
#> 
#>     decompose, spectrum
#> The following object is masked from 'package:base':
#> 
#>     union

# Generate graph by a planted partition model.
g <- sample_ppm(30, 0.3, 0.01, block.sizes=c(10,10,10))
V(g)$color <- c(rep(1,10), rep(2,10), rep(3,10))
plot(g)
```

<img src="man/figures/README-example-1.png" width="100%" />
