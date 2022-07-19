#' Matrix of edge counts between blocks
#'
#' Entries r, s give nubmer of edges between nodes of blocks r and s
#' or, for convenience twice that number if r = s.
#'
#' @param graph An igraph graph.
#' @param partition A vector of group id's for each node of graph.
#' @import igraph

block_edge_counts <- function(graph, partition) {
  CG <- contract(graph, partition)
  E <- as_adjacency_matrix(CG, sparse=FALSE)
  if(!is.directed(graph)) diag(E) <- 2 * diag(E)
  E
}


#' Node counts for blocks
#'
#' @param partition A vector of group id's for each node of graph.
#' @import igraph

block_node_counts <- function(partition) {
  as.integer(table(partition))
}


#' Binary entropy function
#'
#' @param x Numeric, vector or matrix with numbers between 0 and 1.

H_binary <- function (x)
{
  stopifnot(0 <= x & x <= 1)
  H <- - x * log(x) - (1 - x) * log(1 - x)
  replace(H, is.na(H), 0)
}

