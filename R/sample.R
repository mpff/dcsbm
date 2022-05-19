#' Sample planted partition model
#'
#' Sampling from the planted partition model
#'
#' This function samples graphs from a stochastic block model by (doing the
#' equivalent of) Bernoulli trials for each potential edge with the probabilities
#' given by \code{p} for in-group edges and \code{q} for between group edges.
#'
#' @param n Number of vertices in the graph.
#' @param p Probability of creating an edge between vertices of the same group.
#' @param q Probability of creating an edge between vertices of different groups.
#' @param block.sizes Numeric vector giving the number of vertices in each group.
#' The sum of the vector must match the number of vertices.
#' @param directed Logical scalar, whether to generate a directed graph.
#' @param loops Logical scalar, whether self-loops are allowed in the graph.
#' @return An igraph graph.
#' @keywords graphs, sample, planted partition
#' @examples
#'
#' ## Three groups with only a few connection between groups
#' g <- sample_ppm(1000, p=0.3, q=0.01, block.sizes=c(100,600,300))
#' g
#' @export

sample_ppm <- function (n, p, q, block.sizes, directed = FALSE, loops = FALSE)
{
  n <- as.integer(n)
  p <- as.double(p)
  q <- as.double(q)
  block.sizes <- as.integer(block.sizes)
  directed <- as.logical(directed)
  loops <- as.logical(loops)
  pm <- diag(p, length(block.sizes))
  pm[upper.tri(pm) | lower.tri(pm)] <- q
  res <- igraph::sample_sbm(n, pm, block.sizes, directed, loops)
  if (igraph::igraph_opt("add.params")) {
    res$name <- "Planted partition model"
    res$loops <- loops
  }
  res
}
