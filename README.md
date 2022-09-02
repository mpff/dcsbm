
<!-- README.md is generated from README.Rmd. Please edit that file -->

# dcsbm

<!-- badges: start -->

[![R-CMD-check](https://github.com/mpff/dcsbm/workflows/R-CMD-check/badge.svg)](https://github.com/mpff/dcsbm/actions)
[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![Codecov test
coverage](https://codecov.io/gh/mpff/dcsbm/branch/main/graph/badge.svg)](https://app.codecov.io/gh/mpff/dcsbm?branch=main)
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
#> 
#> Attaching package: 'dcsbm'
#> The following object is masked from 'package:igraph':
#> 
#>     sbm

# Generate graph by a planted partition model and estimate an sbm.
g <- sample_dcppm(60, c=0.95, k=20, B=3, k_coef=1.5, directed=TRUE, loops=TRUE)
m <- sbm(g, degree_correction=TRUE, n.sweeps=5, verbose=FALSE)

# Plot result for 3 blocks.
plot(g, vertex.color = m$partition[[8]], vertex.label = NA)
```

<img src="man/figures/README-example-1.png" width="66%" />

``` r

# Plot entropy delta by number of blocks.
plot(m$block_sequence, m$entropy_delta, type = "line")
#> Warning in plot.xy(xy, type, ...): plot type 'line' will be truncated to first
#> character
```

<img src="man/figures/README-example-2.png" width="66%" />
