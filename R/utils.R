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
  blocks <- sort(unique(partition))
  # Order bei number in partition behalten and insert empty rows/cols where needed
  if(max(blocks) != dim(E)[1]){
    miss.blocks <- which(!1:max(blocks) %in% blocks)
    Enew <- diag(0, nrow = max(blocks))
    Enew[-miss.blocks, -miss.blocks] <- E
    E <- Enew
  }
  if(!is.directed(graph)) diag(E) <- 2 * diag(E)
  E
}


#' Node counts for blocks
#'
#' @param partition A vector of group id's for each node of graph.
#' @import igraph

block_node_counts <- function(partition) {
  counts <- table(partition)
  blocks <- sort(unique(partition))
  # Order bei number in partition behalten and insert empty rows/cols where needed
  if(max(blocks) != length(blocks)){
    miss.blocks <- which(!1:max(blocks) %in% blocks)
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


#' Check partition
#'
#' Renumber partitions from 1 to the number of unique blocks (B).
#' @param partition A vector of group id's for each node of graph.
#' @return A vector of group id's going from 1 to B.

check_partition <- function(partition)
{
  partition <- as.integer(partition)
  blocks <- sort(unique(partition))
  in.partition <- 1:max(blocks) %in% partition
  if(!all(in.partition)){
    for(i in 1:length(blocks)) partition[which(partition == blocks[i])] <- i
  }
  partition
}
