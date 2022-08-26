#' Stochastic block model
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
#' g1 <- sample_ppm(60, 0.9, 10, 3)
#' g1 <- simplify(g1)
#' g1 <- delete.vertices(g1, degree(g1) == 0)
#' model <- sbm(g1)
#' @export
#' @import igraph
#' @importFrom utils setTxtProgressBar txtProgressBar

sbm <- function (graph, n.blocks = c(1, Inf), n.moves = 10, n.sweeps = 0,
                 control = list(sigma = 1.5, eps = 0.1, beta = 1, start_partition = NULL))
{
  # Initial graph checks
  stopifnot(is.igraph(graph))
  stopifnot(is.simple(graph))
  stopifnot(all(degree(graph) > 0))

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
  if(!is.null(control$start_partition)) {
    start_partition <- control$start_partition
  } else {
    start_partition <- sample_at_least_once(1:max.blocks, N)
  }
  start_partition <- check_partition(start_partition)
  start_entropy <- get_entropy(graph, start_partition)/E

  # Bookkeeping
  entropy_delta_iter <- 0
  partition_iter <- list(start_partition)

  # Merge sequence
  block_iter <- block_sequence(max.blocks, min.blocks, sigma)
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



block_sequence <- function(Bmax, Bmin = 1, sigma = 1.5) {

  Bmax <- as.integer(Bmax)
  Bmin <- as.integer(Bmin)

  stopifnot(Bmax >= 1, Bmax >= Bmin)
  stopifnot(Bmin >= 1)
  stopifnot(sigma > 1)

  effBmax <- Bmax - Bmin
  effBseq <- c(effBmax)
  effBnext <- floor(effBmax/sigma)
  while (effBnext/sigma > 0) {
    effBseq <- append(effBseq, effBnext)
    effBnext <- floor(effBnext/sigma)
  }
  return(append(effBseq + Bmin, Bmin))
}
