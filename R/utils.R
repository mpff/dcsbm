#' Matrix of edge counts between blocks
#'
#' Entries r, s give nubmer of edges between nodes of blocks r and s
#' or, for convenience twice that number if r = s.
#'
#' @param graph An igraph graph.
#' @param partition A vector of group id's for each node of graph.
#' @import igraph

block_edge_counts <- function(graph, partition) {
  c <- make_clusters(graph, partition, modularity = FALSE)
  # Todo: Can optimize for undirected graph (only consider triangle matrix).
  Ec <- lapply(seq_along(groups(c)), function(r){
    sapply(seq_along(groups(c)), function(s){
      es <- E(graph)[groups(c)[[r]] %->% groups(c)[[s]]]
      if(r == s & !is.directed(graph)){
        2 * length(es) # See SBM model in Piexoto (2014)
      } else {
        length(es)
      }
    })
  })
  do.call("rbind", Ec)
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
