#' Stochastic block model
#'
#' Estimate a ...
#'
#' @param graph An igraph graph.
#' @param partition An vector of integers giving the block partition of nodes.
#' @param degree_correction Whether to use degree correction.
#' @param n.merges Number of block merges in collapse step.
#' @param n.moves Number of merge trials per block.
#' @param n.sweeps Number of sweeps after block merge.
#' @param eps A parameter controlling ...
#' @param beta A parameter controlling ....
#' @param verbose Wether to print verbose output to console.
#' @import igraph
#' @importFrom utils setTxtProgressBar txtProgressBar

collapse_step <- function(graph, partition, degree_correction = FALSE,
                          n.merges = 1, n.moves = 10, n.sweeps = 0,
                          eps = 0.1, beta = 1, verbose = TRUE)
{
  degree_correction = as.logical(degree_correction)
  n.merges <- as.integer(n.merges)
  n.moves <- as.integer(n.moves)
  n.sweeps <- as.integer(n.sweeps)

  partition <- check_partition(partition)
  B.start <- max(partition)

  if(B.start - n.merges < 1){
    warning(paste("Cannot merge", B.start, "blocks", n.moves, "time(s). Skipping merge."))
    return(list("new_partition" = partition, "g1" = NA, "g2" = NA, "entropy_delta" = 0))
  }

  # Block graph
  block.graph <- contract(graph, partition)
  block.partition <- 1:B.start
  block.graph.edges <- block_edge_list(block.graph, block.partition, B.start)

  # Exhaustive search?
  exhaustive_search <- n.moves >= (B.start - 1)
  n.moves <- min(n.moves, B.start)

  # Prepare merge results dataframe
  merge.results <- data.frame(
    "g1" = rep(0, B.start * n.moves),
    "g2" = rep(0, B.start * n.moves),
    "dS" = rep(0, B.start * n.moves)
  )

  # Start partition properties
  old_entropy <- get_entropy(graph, partition, degree_correction)

  # Init Progressbar
  if(verbose) {
    pbmessage <- paste("  [", B.start, "->", B.start - n.merges, "blocks ]")
    pbwidth <- min(80, getOption("width")) - nchar(pbmessage) - 4
    pb = txtProgressBar(min = 0, max = B.start * n.moves, width = pbwidth,
                        initial = 0, style = 3)
    cat(pbmessage)
  }

  for (b in 1:B.start) {
    for (i in 1:n.moves) {

      if (exhaustive_search) {

        # Just go through each block and always accept.
        proposed_merge <- i

      } else {

        # Repeat until new block is proposed (hacky!)
        proposed_merge <- b
        while (b == proposed_merge) {
          # Get a merge proposal using the mcmc sweep move proposal function.
          proposed_merge <- propose_move(V(block.graph)[b], block.graph, 1:B.start, block.graph.edges, eps)
        }

        # If proposed merge was already proposed,
        # we don't need to waste time checking because decision will always
        # result in same state.
        if (any(merge.results[merge.results[,1] == b, 2] == proposed_merge)) next
        if (any(merge.results[merge.results[,1] == proposed_merge, 2] == b)) next

      }

      # Calculate entropy delta for merge and place into results
      proposed_partition <- replace(partition, which(partition == b), proposed_merge)
      proposed_entropy <- get_entropy(graph, proposed_partition, degree_correction)
      entropy_delta <- proposed_entropy - old_entropy

      # Place into results
      merge.idx <- (b-1)*n.moves + i
      merge.results[merge.idx,] <- c(b, proposed_merge, entropy_delta)

    }
    if (verbose) setTxtProgressBar(pb, b * i)
  }
  if (verbose) close(pb)

  # Clean up merge.results
  if(any(merge.results$g1 == merge.results$g2)){
    merge.results <- merge.results[-which(merge.results$g1 == merge.results$g2),]
  }
  merge_results_ordered <- merge.results[order(merge.results$dS, decreasing = FALSE),]

  n.merged <- 0
  new_partition <- partition

  while (n.merged < n.merges) {
    if (length(merge_results_ordered[1,]) == 0){
      warning(paste("Ran out of unique blocks to merge. Skipping after", n.merged, "merge(s)."))
      break
    }

    # Get best merge pair
    best_merge_pair <- merge_results_ordered[1,]

    # Update results list. Remove merges of same blocks.
    merge_results_ordered <- merge_results_ordered[-1,]

    # Update new partition
    new_partition <- replace(new_partition, new_partition == best_merge_pair$g1, best_merge_pair$g2)

    #message(best_merge_pair$g1, " -> ", best_merge_pair$g2, " (", length(unique(new_partition)), ")")

    # Update merge results
    merge_results_ordered$g1 <- replace(merge_results_ordered$g1, which(merge_results_ordered$g1 == best_merge_pair$g1), best_merge_pair$g2)
    merge_results_ordered$g2 <- replace(merge_results_ordered$g2, which(merge_results_ordered$g2 == best_merge_pair$g1), best_merge_pair$g2)
    merge_results_ordered <- merge_results_ordered[which(merge_results_ordered$g1 != merge_results_ordered$g2),]

    # Update n.merged iterator
    n.merged <- n.merged + 1
  }

  new_entropy_delta <- get_entropy(graph, new_partition, degree_correction) - old_entropy
  new_partition <- check_partition(new_partition)

  # Perform mcmc sweeps to settle partition (costly)
  if (n.sweeps > 0) {
    sweep_results <- mcmc_sweep(graph, new_partition, max(new_partition), degree_correction, n.sweeps, eps, beta)
    new_partition <- sweep_results$best_partition
    new_entropy_delta <- new_entropy_delta + sweep_results$best_entropy_delta
    new_partition <- check_partition(new_partition)
  }

  list("new_partition" = new_partition, "entropy_delta" = new_entropy_delta)
}
