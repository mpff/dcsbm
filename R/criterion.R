#' Description length of a block partition
#'
#' Calculate the description length associated with the current block partition
#' and type of degree correction. Can be used with directed and undirected graphs.
#'
#' @param graph An igraph graph.
#' @param partition Vector of integer values giving the block membership of each
#' vertex
#' @param degree_correction Type of degree correction to use. "oneway" for one
#' parameter per vertex, "twoway" for two parameters per vertex (input/output),
#' "none" for no degree correction.
#' @return Description length value (numeric) for the given graph and partition.
#' @keywords graphs, inference, stochastic block model, degree correction
#' @export
#' @import igraph

calculate_dl <- function(graph, partition, degree_correction = FALSE) {
  # Initial checks
  stopifnot(is.igraph(graph))
  partition <- check_partition(partition)
  degree_correction <- as.logical(degree_correction)

  # Get params
  directed <- is.directed(graph)
  E <- length(E(graph))
  N <- length(V(graph))
  B <- length(unique(partition))

  # Get information necessary to describe the model

  # Traditional information
  if (directed) {
    information <- E * h_func(B * B / E) + N * log(B)
  } else {
    information <- E * h_func(B * (B + 1) / (2 * E)) + N * log(B)
  }

  # Adjust for degree corrected information
  if (degree_correction) {
    if (directed) {
      deg_in <- degree(graph, mode = "in")
      deg_in_counts <- table(deg_in)
      deg_ou <- degree(graph, mode = "out")
      deg_ou_counts <- table(deg_ou)
      information <- information - N * sum(deg_in_counts / N * log(deg_in_counts / N))
      information <- information - N * sum(deg_ou_counts / N * log(deg_ou_counts / N))
    } else {
      deg <- degree(graph)
      deg_counts <- table(deg)
      information <- information - N * sum(deg_counts / N * log(deg_counts / N))
    }
  }

  # Get entropy of the current partiton
  entropy <- get_entropy(graph, partition, degree_correction)

  # Calculate MDL
  entropy + information
}


h_func <- function(x) {
  (1 + x) * log(1 + x) - x * log(x)
}
