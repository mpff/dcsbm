#' Sample planted partition model
#'
#' Generate a random graph using the planted partition model
#'
#' This function samples graphs from a stochastic block model by sampling
#' each potential transmission with the probabilities given by \code{p} for
#' in-group transmissions and \code{q} for between group transmissions and
#' group membership given by the vector \code{partitions}.
#'
#' @param partition Vector of group memberships.
#' @param p Edge probability between nodes of same group.
#' @param q Edge probability between nodes of different groups.
#' @return An igraph graph.
#' @keywords graphs, sample, planted partition
#' @examples
#'
#' ## Three groups with only a few connection between groups
#' partition <- sample(1:3, 50, replace = TRUE)
#' g <- sample_ppm(partition, 0.3, 0.01)
#' g
#' @export
#' @importFrom stats runif

sample_ppm <- function (partition, p, q)
{
  partition <- structure(as.factor(partition))
  p <- as.double(p)
  q <- as.double(q)
  n <- length(partition)
  nodes <- seq_len(n)
  al <- lapply(nodes, function(i){
    # Adapted from https://bldavies.com/blog/generating-random-graphs-communities
    prob <- c(q, p)[1 + (partition[i] == partition)]
    nodes[which(runif(n) < prob)]
  })
  res <- igraph::graph_from_adj_list(al)
  res$name <- "Planted partition model"
  res
}
