% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model.R
\name{block_sequence}
\alias{block_sequence}
\title{Build a sequence of block numbers}
\usage{
block_sequence(Bmax, Bmin = 1, sigma = 1.5)
}
\arguments{
\item{Bmax}{Start value}

\item{Bmin}{End value}

\item{sigma}{Division factor.}
}
\value{
A vector of block numbers in decreasing order.
}
\description{
So that B_{i+1} = B_i/sigma.
}
