
<!-- README.md is generated from README.Rmd. Please edit that file -->

# dcsbm

<!-- badges: start -->
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
library(igraph)
#> Warning: package 'igraph' was built under R version 4.1.3
#> 
#> Attaching package: 'igraph'
#> The following objects are masked from 'package:stats':
#> 
#>     decompose, spectrum
#> The following object is masked from 'package:base':
#> 
#>     union

# Generate graph by a planted partition model.
partition <- sample(1:3, 30, replace = T)
g <- sample_ppm(partition, 0.25, 0.01)
V(g)$color <- partition
plot(g)
```

<img src="man/figures/README-example-1.png" width="100%" />
