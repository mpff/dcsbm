% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sample.R
\name{sample_ppm}
\alias{sample_ppm}
\title{Sample planted partition model (Piexoto 2014)}
\usage{
sample_ppm(N, c, k, B, directed = FALSE, loops = FALSE)
}
\arguments{
\item{N}{Number of vertices in the graph.}

\item{c}{A number between 0 and 1 controlling the}

\item{k}{The average degree of a vertex.}

\item{B}{The number of blocks.}

\item{directed}{Logical scalar, whether to generate a directed graph.}

\item{loops}{Logical scalar, whether self-loops are allowed in the graph.}
}
\value{
An igraph graph.
}
\description{
Sampling from the planted partition model as in Piexoto 2014.
}
\details{
This function samples graphs from a stochastic block model by building a
block adjacency matrix from the Parameters N, k, B and c.
}
\examples{
## Three groups with only a few connection between groups
G <- sample_ppm(300, c=0.9, k=10, B=3)
p <- c(rep(1,100), rep(2, 200), rep(3, 300))
plot(G, vertex.label=NA, vertex.color=p)
}
\keyword{graphs,}
\keyword{partition}
\keyword{planted}
\keyword{sample,}
