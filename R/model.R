#' Degree corrected stochastic block model
#'
#' Estimate a ...
#'
#' @param graph An igraph graph.
#' @param n.blocks Number of blocks.
#' @param n.moves Number of merge trials per block.
#' @param n.sweeps Number of sweeps after block merge.
#' @param control List of parameters for the inference algorithm.
#' @return Todo
#' @keywords graphs, inference, stochastic block model, degree correction
#' @examples
#' ## Three groups with weighted connections.
#' g1 <- sample_ppm2(300, 0.9, 10, 3)
#' g1 <- simplify(g1)
#' g1 <- delete.vertices(g1, degree(g1) == 0)
#' model <- sbm(g1)
#' @export
#' @import igraph
#' @importFrom utils setTxtProgressBar txtProgressBar

sbm <- function (graph, n.blocks = c(1, Inf), n.moves = 10, n.sweeps = 0,
                 control = list(sigma = 1.5, eps = 0.1, beta = 1))
{
  # Initial graph checks
  stopifnot(is.igraph(graph))
  stopifnot(is.simple(graph))
  stopifnot(all(degree(graph) > 0))
  stopifnot(!is.directed(graph))

  # Initial parameter checks
  n.moves = as.integer(n.moves)
  n.sweeps = as.integer(n.sweeps)

  # Check n.blocks
  min.blocks <- as.integer(n.blocks[1])
  if(!is.infinite(n.blocks[2])){
    max.blocks <- as.integer(n.blocks[2])
  } else {
    max.blocks <- length(V(graph))
  }
  stopifnot(min.blocks <= max.blocks)

  # Check control parameters
  sigma <- control$sigma
  eps <- control$eps
  beta <- control$beta
  stopifnot(sigma > 1, eps > 0, beta > 0)

  # Graph parameters
  N <- length(V(graph))
  stopifnot(max.blocks <= N)
  E <- length(E(graph))

  # Start partition
  start_partition <- sample_at_least_once(1:max.blocks, N)
  start_partition <- check_partition(start_partition)
  start_entropy <- get_entropy(graph, start_partition)/E

  # Bookkeeping
  entropy_delta_iter <- 0
  partition_iter <- list(start_partition)

  # Merge sequence
  block_iter <- block_sequence(max.blocks, sigma, min.blocks)
  merges_iter <- -1 * diff(block_iter)

  # Loop trough all merges
  for(curr.iter in 1:length(merges_iter)) {
    collapse_result <- collapse_step(graph, partition_iter[[curr.iter]],
                                     n.merges = merges_iter[curr.iter],
                                     n.moves = n.moves, n.sweeps = n.sweeps,
                                     eps = eps, beta = beta)
    partition_iter[[curr.iter+1]] <- collapse_result$new_partition
    entropy_delta_iter <- append(entropy_delta_iter, collapse_result$entropy_delta/E)
  }


  list("block_sequence" = block_iter, "entropy_delta" = entropy_delta_iter,
       "partitions" = partition_iter)
}


block_sequence <- function(B, sigma, min.B = 1) {

  B <- as.integer(B)

  stopifnot(sigma > 1)
  stopifnot(B >= 1)

  B.seq <- c(B)
  next.B <- floor(B/sigma)
  while (next.B/sigma > min.B) {
    B.seq <- append(B.seq, next.B)
    next.B <- floor(next.B/sigma)
  }
  return(append(B.seq, min.B))
}
