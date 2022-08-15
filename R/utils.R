#' Matrix of edge counts between blocks
#'
#' Entries r, s give nubmer of edges between nodes of blocks r and s
#' or, for convenience twice that number if r = s.
#'
#' @param graph An igraph graph.
#' @param partition A vector of group id's for each node of graph.
#' @param n.blocks The number of groups. Inferred from the
#' number of individual group id's in partition by default.
#' @import igraph

block_edge_counts <- function(graph, partition, n.blocks = NULL) {

  CG <- contract(graph, partition)
  E <- as_adjacency_matrix(CG, sparse=FALSE)

  # Add blocks that are currently not represented in partition.
  if ( !is.null(n.blocks) ) {
    if( dim(E)[1] < n.blocks ){
      miss.blocks <- (dim(E)[1]+1):n.blocks
      Enew <- diag(0, nrow = n.blocks)
      Enew[-miss.blocks, -miss.blocks] <- E
      E <- Enew
    }
  }

  # Count incoming and outgoing edges on diagonal
  if(!is.directed(graph)) diag(E) <- 2 * diag(E)

  # Return
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

#' Sample at least once with replacement
#' See https://stackoverflow.com/a/26350427/10042003.
#' @param x a numeric or vector from which to sample.
#' @param size a non-negative integer giving the number of items to choose.
#' @param prob a vector of probability weights for obtaining the elements of the
#' vector being sampled.
#' @value A vector of length \code{size} with elements drawn from \code{x} where
#' every element of \code{x} occurs at least once if \code{size} is larger than
#' \code{length(x)}.

sample_at_least_once <- function(x, size, prob = NULL){
  # Only consider unique items.
  items <- unique(x)
  if(length(items) > size){
    stop("Not enough unique items in input to sample at least one of each.")
  }

  # Get values for vector - force each item in at least once
  vals <- c(items, sample(items, size - length(items), replace = TRUE, prob = prob))

  # Now shuffle them
  sample(vals)
}


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
