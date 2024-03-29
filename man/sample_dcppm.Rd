% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sample.R
\name{sample_dcppm}
\alias{sample_dcppm}
\title{Sample planted partition model with degree variability (Piexoto 2020)}
\usage{
sample_dcppm(N, c, k, B, k_coef = 0, directed = FALSE, loops = FALSE)
}
\arguments{
\item{N}{Number of vertices in the graph.}

\item{c}{A number between 0 and 1 controlling the}

\item{k}{The average degree of a vertex.}

\item{B}{The number of blocks.}

\item{k_coef}{A number describing the degree variability within each block.
A value of 0 means that all vertices have the same expected degree.#'}

\item{directed}{Logical scalar, whether to generate a directed graph.}

\item{loops}{Logical scalar, whether self-loops are allowed in the graph.}
}
\value{
An igraph graph.
}
\description{
This function samples graphs from a stochastic block model by building a
block adjacency matrix from the Parameters N, k, B and c and accounting
for a exponential degree sequence inside each block controlled by \code{k_coef}.
}
\examples{
## Three groups with only a few connection between groups
G <- sample_dcppm(300, c=0.9, k=10, B=3, k_coef=2)
p <- c(rep(1,100), rep(2, 200), rep(3, 300))
plot(G, vertex.label=NA, vertex.color=p)
}
\keyword{graphs,}
\keyword{partition}
\keyword{planted}
\keyword{sample,}
