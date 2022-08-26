#' Entropy of a block partition
#'
#' Calculate the entropy associated with the current block partition and type
#' of degree correction. Can be used with directed and undirected graphs.
#'
#' @param graph An igraph graph.
#' @param partition Vector of integer values giving the block membership of each
#' vertex
#' @param degree_correction Type of degree correction to use. "oneway" for one
#' parameter per vertex, "twoway" for two parameters per vertex (input/output),
#' "none" for no degree correction.
#' @return Entropy value (numeric) for the given graph and partition.
#' @keywords graphs, inference, stochastic block model, degree correction
#' @export
#' @import igraph

get_entropy <- function (graph, partition, degree_correction = FALSE)
{
  # Initial checks
  stopifnot(is.igraph(graph))
  partition <- as.integer(partition)
  degree_correction <- as.logical(degree_correction)
  stopifnot(length(graph) == length(partition))

  # Add vertices with degree zero to first group.
  zero.vertices <- which(degree(graph) == 0)
  partition[zero.vertices] <- 1

  # Renumber blocks in partition.
  partition <- check_partition(partition)

  # Calculate Entropy
  if(!degree_correction) {
    E <- block_edge_counts(graph, partition)
    n <- block_node_counts(partition)
    S <- entropy_trad(E, n, directed = is.directed(graph))
  } else {
    E <- block_edge_counts(graph, partition)
    d <- degree(graph)
    S <- entropy_corrected(E, d, directed = is.directed(graph))
  }
  S
}


#' Entropy (no degree correction)
#'
#' Calculate the entropy associated with the current block partition.
#'
#' \deqn{S_t = \frac{1}{2} \sum_{rs} n_r n_s H_b\left(\frac{e_{rs}}{n_r n_s}\right)}{St = 0.5 * sum_(n[r] * n[s] * H_b(E[r,s] / (n[r]*n[s]))}
#'
#' @param E Integer matrix of edge counts associated with current partition
#' @param n Integer vector of node counts associated with current partition
#' @param directed Whether graph is directed
#' @import igraph

entropy_trad <- function (E, n, directed = FALSE)
{
  if(any(n == 0)){
    miss.blocks <- which(n == 0)
    E <- E[-miss.blocks, -miss.blocks]
    n <- n[-miss.blocks]
  }
  Enn <- (E / n) %*% diag(1/n, nrow = length(n))
  Hmat <- H_binary(Enn)
  S <- .5 * matrix(n, nrow = 1) %*% Hmat %*% matrix(n)
  if(directed) S <- 2 * S
  as.numeric(S)
}


#' Entropy (degree corrected)
#'
#' Calculate the entropy associated with the current block partition.
#'
#' \deqn{
#'   S_c = - E - \sum_k N_k \ln k! - \frac{1}{2} \sum_{rs} e_{rs} \ln\left(\frac{e_{rs}}{e_r e_s}\right)
#' }{
#'   St = - E - sum(N[k] * ln k!) - 0.5 * sum(e[r,s] * ln(e[r,s] / (e[r]*e[s])))
#' }
#'
#' @param E Integer matrix of edge counts associated with current partition.
#' @param d Integer vector of node degrees.
#' @param directed Whether graph is directed.
#' @import igraph

entropy_corrected <- function (E, d, directed = FALSE)
{
  # Calcualte part of entropy relating to total amount of edges
  edge_entropy <- sum(E)/2

  # Calculate part of entropy relating to different degrees
  degree_unique <- sort(unique(d))
  degree_counts <- table(d)
  degree_entropy <- sum(degree_counts * log(factorial(degree_unique)))

  # Calculate part of entropy relating to block edge counts
  e <- colSums(E)
  Eee <- (E / e) %*% diag(1/e, nrow = length(e))
  block_edges_entropy <- 0.5 * sum(E * log(Eee))

  # Add everything together
  S <- - edge_entropy - degree_entropy - block_edges_entropy
  if(directed) stop("Entropy for these parameters not yet implemented.")

  as.numeric(S)
}
