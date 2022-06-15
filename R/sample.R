#' Sample planted partition model
#'
#' Sampling from the planted partition model
#'
#' This function samples graphs from a stochastic block model by (doing the
#' equivalent of) \code{n.trial} Bernoulli trials for each potential edge with
#' probabilities given by \code{p} for in-group edges and \code{q} for between
#' group edges.
#'
#' @param n Number of vertices in the graph.
#' @param p Probability of creating an edge between vertices of the same group.
#' @param q Probability of creating an edge between vertices of different groups.
#' @param block.sizes Numeric vector giving the number of vertices in each group.
#' The sum of the vector must match the number of vertices.
#' @param n.trials Number of repeated Bernoulli trials per edge. Graphs are
#' weighted for \code{n.trials} greater than 1.
#' @param directed Logical scalar, whether to generate a directed graph.
#' @param loops Logical scalar, whether self-loops are allowed in the graph.
#' @return An igraph graph.
#' @keywords graphs, sample, planted partition
#' @examples
#'
#' ## Three groups with only a few connection between groups
#' g1 <- sample_ppm(1000, p=0.3, q=0.01, block.sizes=c(100,600,300))
#' g1
#' ## Three groups with weighted connections.
#' g2 <- sample_ppm(1000, p=0.1, q=0.003, block.sizes=c(100,600,300), n.trials=3)
#' g2
#' @export
#' @import igraph

sample_ppm <- function (n, p, q, block.sizes,
                        n.trials = 1, directed = FALSE, loops = FALSE)
{
  n <- as.integer(n)
  p <- as.double(p)
  q <- as.double(q)
  block.sizes <- as.integer(block.sizes)
  n.trials <- as.integer(n.trials)
  directed <- as.logical(directed)
  loops <- as.logical(loops)
  pm <- diag(p, length(block.sizes))
  pm[upper.tri(pm) | lower.tri(pm)] <- q
  res <- sample_sbm(n, pm, block.sizes, directed, loops)
  if(n.trials > 1){
    # See https://stackoverflow.com/a/27766128 (Combine two graphs with weights)
    E(res)$weight <- 1
    for (i in seq(n.trials - 1)){
      resnew <- sample_sbm(n, pm, block.sizes, directed, loops)
      E(resnew)$weight <- 1
      res <- union(res, resnew)
      E(res)$weight_1[is.na(E(res)$weight_1)] <- 0
      E(res)$weight_2[is.na(E(res)$weight_2)] <- 0
      E(res)$weight <- E(res)$weight_1 + E(res)$weight_2
      res <- remove.edge.attribute(res, "weight_1")
      res <- remove.edge.attribute(res, "weight_2")
      if (igraph_opt("add.params")) {
        res$name <- res$name_1
        res <- remove.graph.attribute(res, "name_1")
        res <- remove.graph.attribute(res, "name_2")
        res$loops <- res$loops_1
        res <- remove.graph.attribute(res, "loops_1")
        res <- remove.graph.attribute(res, "loops_2")
      }
    }
  }
  if (igraph_opt("add.params")) {
    res$name <- "Planted partition model"
    if(n.trials > 1) res$name <- "Weighted planted partition model"
    res$loops <- loops
  }
  res
}
