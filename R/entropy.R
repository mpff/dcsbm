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
  # Parameters.
  partition <- as.integer(partition)
  degree_correction <- as.logical(degree_correction)

  # Initial checks.
  stopifnot(is.igraph(graph))
  stopifnot(length(graph) == length(partition))

  # Add vertices with degree zero to first group. (TODO: Can this be removed?)
  zero.vertices <- which(degree(graph) == 0)
  partition[zero.vertices] <- 1

  # Renumber blocks in partition.
  partition <- check_partition(partition)

  # Calculate Entropy
  if(!degree_correction) {
    E <- block_edge_counts(graph, partition)
    n <- block_node_counts(partition)
    S <- entropy_trad(E, n, directed = is.directed(graph), simple = is.simple(graph))
  } else {
    E <- block_edge_counts(graph, partition)
    if (!is.directed(graph)) {
      d <- list("total" = degree(graph))
    } else {
      d <- list("in" = degree(graph, mode = "in"), "out" = degree(graph, mode = "out"))
    }
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
#' @param simple Whether graph is simple.
#' @import igraph

entropy_trad <- function (E, n, directed = FALSE, simple = TRUE)
{
  if(any(n == 0)){
    miss.blocks <- which(n == 0)
    E <- E[-miss.blocks, -miss.blocks]
    n <- n[-miss.blocks]
  }
  if (simple) {
    Enn <- (E / n) %*% diag(1/n, nrow = length(n))
    Hmat <- H_binary(Enn)
    S <- .5 * matrix(n, nrow = 1) %*% Hmat %*% matrix(n)
  } else {
    nn <- n %*% t(n)
    Hmat <- H_binary(nn/(nn + E))
    S <- .5 * sum((nn + E) * Hmat)
  }
  if(directed) S <- 2 * S
  as.numeric(S)
}


#' Entropy (degree corrected)
#'
#' Calculate the entropy associated with the current block partition.
#'
#' Undirected case:
#'
#' \deqn{
#'   S_c = - E - \sum_k N_k \ln k! - \frac{1}{2} \sum_{rs} e_{rs} \ln\left(\frac{e_{rs}}{e_r e_s}\right)
#' }{
#'   St = - E - sum(N[k] * ln k!) - 0.5 * sum(e[r,s] * ln(e[r,s] / (e[r]*e[s])))
#' }
#'
#' @param E Integer matrix of edge counts associated with current partition.
#' @param d Integer vector of node degrees.
#' @param directed Wether to calcualted directed entropy.
#' @import igraph

entropy_corrected <- function (E, d, directed = FALSE)
{

   if (!directed) {

    # Calcualte part of entropy relating to total amount of edges
    edge_entropy <- sum(E)/2  # /2 because of undirected edges

    d <- d$total

    # Calculate part of entropy relating to different degrees
    degree_unique <- sort(unique(d))
    degree_counts <- table(d)
    degree_entropy <- sum(degree_counts * log(factorial(degree_unique)))
    # TODO: Check if a factor "degree_unique" is missing here?!

    # Calculate part of entropy relating to block edge counts
    e <- colSums(E)
    Eee <- (E / e) %*% diag(1/e, nrow = length(e))
    block_edges_entropy <- 0.5 * sum(E * log(Eee), na.rm = TRUE)

    # Add everything together
    S <- - edge_entropy - degree_entropy - block_edges_entropy

  } else {

    # Calcualte part of entropy relating to total amount of edges
    edge_entropy <- sum(E)

    d_in <- d$`in`
    d_ou <- d$out

    # Calculate part of entropy relating to different degrees
    degree_in_unique <- sort(unique(d_in))
    degree_in_counts <- table(d_in)
    degree_in_entropy <- sum(degree_in_counts * log(factorial(degree_in_unique)))
    degree_ou_unique <- sort(unique(d_ou))
    degree_ou_counts <- table(d_ou)
    degree_ou_entropy <- sum(degree_ou_counts * log(factorial(degree_ou_unique)))

    # Calculate part of entropy relating to block edge counts
    e_in <- colSums(E)
    e_ou <- rowSums(E)
    Eee <- (E / e_in) %*% diag(1/e_ou, nrow = length(e_ou))
    block_edges_entropy <- sum(E * log(Eee), na.rm = TRUE)

    # Add everything together
    S <- - edge_entropy - degree_in_entropy - degree_ou_entropy - block_edges_entropy

  }

  as.numeric(S)
}
