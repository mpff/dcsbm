#' Stochastic block model with degree correction
#'
#' Estimate a stochastic block model with or without degree correction. Works
#' for directed multigraphs with self-loops.
#'
#' @param graph An igraph graph.
#' @param degree_correction Whether to use degree correction.
#' @param n.blocks Number of blocks.
#' @param n.moves Number of merge trials per block.
#' @param n.sweeps Number of sweeps after block merge.
#' @param verbose Wether to print verbose output during estimation.
#' @param control List of parameters for the inference algorithm.
#' @return Todo
#' @keywords graphs, inference, stochastic block model, degree correction
#' @examples
#' ## Three groups with weighted connections.
#' g1 <- sample_ppm(60, 0.9, 10, 3)
#' model <- dcsbm(g1)
#' @export
#' @import igraph

dcsbm <- function (graph, degree_correction = FALSE, n.blocks = c(1, Inf),
                 n.moves = 10, n.sweeps = 0, verbose = TRUE,
                 control = list(sigma = 1.5, eps = 0.1, beta = 1, start_partition = NULL))
{
  # Initial graph checks
  stopifnot(is.igraph(graph))
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

  if (verbose) {
    # Graph Info
    graph_msg <- "\nGRAPH:"
    graph_msg <- paste(graph_msg, ifelse(is.directed(graph), "Directed -", "Undirected -"))
    graph_msg <- paste(graph_msg, ifelse(is.simple(graph), "Simple graph -", "Multi graph -"))
    graph_msg <- paste(graph_msg, ifelse(any(is.loop(graph)), "With loops", "Without loops"))
    graph_msg <- paste0(graph_msg, " (N = ", N, ",", " E = ", E, ")")
    cat(graph_msg)

    # SBM Info
    sbm_msg <- "\n\nEstimating a"
    sbm_msg <- paste(sbm_msg, ifelse(degree_correction, "degree corrected", "traditional"))
    if(degree_correction & is.directed(graph)) sbm_msg <- paste(sbm_msg, "(Input/Output)")
    sbm_msg <- paste(sbm_msg, "SBM via agglomerative merging.\n")
    cat(sbm_msg)
  }

  # Start partition
  part_msg <- "\nUsing"
  if(!is.null(control$start_partition)) {
    part_msg <- paste(part_msg, "the provided starting partition")
    start_partition <- control$start_partition
    start_partition <- check_partition(start_partition)
  } else {
    part_msg <- paste(part_msg, "a random starting partition with", max.blocks, "blocks.\n")
    start_partition <- sample_at_least_once(1:max.blocks, N)
  }
  if (verbose) cat(part_msg)

  # Start partition parameters
  start_description_length <- calculate_dl(graph, start_partition, degree_correction)/E
  start_entropy <- get_entropy(graph, start_partition, degree_correction)/E

  # Bookkeeping
  collapse_iterations <- list(list(
    "number_of_blocks" = max.blocks,
    "entropy" = start_entropy,
    "description_length" = start_description_length,
    "partition" = start_partition
  ))

  # Has converged?
  converged <- FALSE

  # Loop until convergence
  n_iterations <- 1
  while (!converged) {

    # Current search parameters
    current_block_sequence <- block_sequence(max.blocks, min.blocks, sigma)
    current_partition_sequence <- list(start_partition)
    current_merge_sequence <- -1 * diff(current_block_sequence)
    current_dl_sequence <- calculate_dl(graph, start_partition, degree_correction)/E

    # Stop after this iteration if search is exhaustive
    if(all(current_merge_sequence == 1)) converged <- TRUE

    # Loop trough all merges
    if (verbose) cat(paste("\nMerging from", max.blocks, "to", min.blocks, "block(s)...\n"))
    for(i in 1:length(current_merge_sequence)) {
      collapse_result <- collapse_step(graph, current_partition_sequence[[i]],
                                       degree_correction = degree_correction,
                                       n.merges = current_merge_sequence[i],
                                       n.moves = n.moves, n.sweeps = n.sweeps,
                                       eps = eps, beta = beta, verbose = verbose)

      # Global bookkeeping.
      collapse_iterations[[n_iterations+i+1]] <- list(
        "number_of_blocks" = current_block_sequence[i+1],
        "entropy" = collapse_result$new_entropy/E,
        "description_length" = collapse_result$description_length/E,
        "partition" = collapse_result$new_partition
      )

      # Current collapse loop bookkeeping.
      current_partition_sequence[[i+1]] <- collapse_result$new_partition
      current_dl_sequence[i+1] <- collapse_result$description_length/E

    }

    # Update number of iterations
    n_iterations <- n_iterations + length(current_merge_sequence)

    # Minimum description length of current agglom. merge iteration
    minimum_dl <- min(current_dl_sequence)
    minimum_dl_iteration <- which(current_dl_sequence == minimum_dl)
    best_number_of_blocks <- current_block_sequence[minimum_dl_iteration]
    best_partition <- current_partition_sequence[[minimum_dl_iteration]]

    # Update for next agglom. search.
    max.blocks <- current_block_sequence[max(1, minimum_dl_iteration - 1)]
    min.blocks <- current_block_sequence[min(length(current_block_sequence), minimum_dl_iteration + 1)]
    start_partition <- current_partition_sequence[[max(1, minimum_dl_iteration - 1)]]

    if (verbose) {
      dl_msg <- "Minimal description length for"
      dl_msg <- paste(dl_msg,  best_number_of_blocks)
      dl_msg <- paste(dl_msg, "blocks")
      dl_msg <- paste0(dl_msg, " (DL/E = ", round(minimum_dl, digits=2) ,").\n")
      cat(dl_msg)
    }

  }

  if (verbose) if (verbose) {
    dl_msg <- "\nConverged at"
    dl_msg <- paste(dl_msg,  best_number_of_blocks)
    dl_msg <- paste(dl_msg, "blocks")
    dl_msg <- paste0(dl_msg, " (DL/E = ", round(minimum_dl, digits=2) ,").\n\n")
    cat(dl_msg)
  }

  degree_parameters <- block_degree_sequence(graph, best_partition)

  list("block_transmission_probs" = NULL, "B_opt" = best_number_of_blocks,
       "minimum_description_length" = minimum_dl, "best_partition" = best_partition,
       "degree_parameters" = degree_parameters,
       "iterations" = collapse_iterations)
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
