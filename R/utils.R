#' Matrix of edge counts between blocks
#'
#' Entries r, s give nubmer of edges between nodes of blocks r and s
#' or, for convenience twice that number if r = s.
#'
#' @param graph An igraph graph.
#' @param partition A vector of group id's for each node of graph.
#' @import igraph

block_edge_counts <- function(graph, partition, n.blocks = NULL) {
  CG <- contract(graph, partition)
  E <- as_adjacency_matrix(CG, sparse=FALSE)
  if(is.null(n.blocks)) n.blocks <- max(partition)
  if(n.blocks != dim(E)[1]){
    counts <- table(partition)
    miss.blocks <- which(!1:n.blocks %in% dimnames(counts)$partition)
    Enew <- diag(0, nrow = n.blocks)
    Enew[-miss.blocks, -miss.blocks] <- E
    E <- Enew
  }
  if(!is.directed(graph)) diag(E) <- 2 * diag(E)
  E
}


#' Node counts for blocks
#'
#' @param partition A vector of group id's for each node of graph.
#' @param n.blocks The number of groups. Inferred from the
#' number of individual group id's in partition by default.
#' @import igraph

block_node_counts <- function(partition, n.blocks = NULL) {
  counts <- table(partition)
  if(is.null(n.blocks)){
    n.blocks <- max(names(counts))
  }
  if(n.blocks != dim(counts)){
    miss.blocks <- which(!1:n.blocks %in% dimnames(counts)$partition)
    for(block in miss.blocks) counts <- append(counts, 0, after = block - 1)
  }
  as.integer(counts)
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

#' Safer resample
#' @param x Numeric or vector from which to resample.
#' @param ... other arguments of \code{sample}.

resample <- function(x, ...) x[sample.int(length(x), ...)]
