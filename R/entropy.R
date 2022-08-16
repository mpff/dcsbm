#' Entropy of a block partition
#'
#' Calculate the entropy associated with the current block partition and type
#' of degree correction. Can be used with directed and undirected graphs.
#'
#' @param graph An igraph graph.
#' @param partition Vector of integer values giving the block membership of each
#' vertex
#' @param degree_correction Type of degree correction to use. "oneway" for one
#' parameter per vertex, "twoways" for two parameters per vertex (input/output),
#' "none" for no degree correction.
#' @return Entropy value (numeric) for the given graph and partition.
#' @keywords graphs, inference, stochastic block model, degree correction
#' @examples
#' ## Three groups with weighted connections.
#' g1 <- sample_ppm(100, 0.3, 0.03, block.sizes = c(30, 50, 20))
#' ## random partition
#' p1 <- sample(c(1,2,3), 100, replace = TRUE)
#' get_entropy(g1, p1)
#' @export
#' @import igraph

get_entropy <- function (graph, partition,
                         degree_correction = c("none", "oneway", "twoways"))
{
  # Initial checks
  stopifnot(is.igraph(graph))
  partition <- as.integer(partition)
  degree_correction <- match.arg(degree_correction)
  stopifnot(length(graph) == length(partition))
  # Add vertices with degree zero to first group.
  zero.vertices <- which(degree(graph) == 0)
  partition[zero.vertices] <- 1
  # Renumber blocks in partition.
  partition <- check_partition(partition)
  # Calculate Entropy
  if(degree_correction == "none") {
    E <- block_edge_counts(graph, partition)
    n <- block_node_counts(partition)
    S <- entropy_trad(E, n, directed = is.directed(graph))
  } else stop("Entropy for these parameters not yet implemented.")
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
