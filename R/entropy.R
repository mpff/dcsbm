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
#' @param directed Logical scalar, Whether to use the directed or undirected
#' entropy. Ignored for undirected graphs.
#' @return Entropy value (numeric) for the given graph and partition.
#' @keywords graphs, inference, stochastic block model, degree correction
#' @examples
#' ## Three groups with weighted connections.
#' g1 <- sample_ppm(100, 0.3, 0.03, block.sizes = c(30, 50, 20))
#' ## random partition
#' p1 <- sample(c(1,2,3), 100, replace = TRUE)
#' get_entropy(g1, p1, degree_correction = "none", directed = FALSE)
#' @export
#' @import igraph

get_entropy <- function (graph, partition,
                         degree_correction = c("none", "oneway", "twoways"))
{
  stopifnot(is.igraph(graph))
  partition <- as.integer(partition)
  stopifnot(length(graph) == length(partition))
  degree_correction <- match.arg(degree_correction)
  if(!is.directed(graph)){
    if(degree_correction == "none") {
      E <- block_edge_counts(graph, partition)
      n <- block_node_counts(partition)
      S <- entropy_undirected_trad(E, n)
    }
    else stop("Entropy for these parameters not yet implemented.")
  }
  else stop("Entropy for directed graphs not yet implemented.")
  S
}


#' Entropy (undirected, no degree correction)
#'
#' Calculate the entropy associated with the current block partition.
#'
#' \deqn{S_t = \frac{1}{2} \sum_{rs} n_r n_s H_b\left(\frac{e_{rs}}{n_r n_s}\right)}{St = 0.5 * sum_(n[r] * n[s] * H_b(E[r,s] / (n[r]*n[s]))}
#'
#' @param E Integer matrix of edge counts associated with current partition
#' @param n Integer vector of node counts associated with current partition
#' @import igraph

entropy_undirected_trad <- function (E, n)
{
  Enn <- (E / n) %*% diag(1/n)
  Hmat <- H_binary(Enn)
  S <- .5 * matrix(n, nrow = 1) %*% Hmat %*% matrix(n)
  as.numeric(S)
}
