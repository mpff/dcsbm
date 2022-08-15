


collapse_step <- function(graph, partition, n.merges = 1, n.moves = 10, n.sweeps = 0, eps = 0.1, beta = 1)
{
  n.merges <- as.integer(n.merges)
  n.moves <- as.integer(n.moves)
  n.sweeps <- as.integer(n.sweeps)

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
  exhaustive_search <- n.moves >= B.start
  n.moves <- min(n.moves, B.start)

  # Prepare merge results dataframe
  merge.results <- data.frame(
    "g1" = rep(0, B.start * n.moves),
    "g2" = rep(0, B.start * n.moves),
    "dS" = rep(0, B.start * n.moves)
  )

  # Start partition properties
  old_entropy <- get_entropy(graph, partition)

  # Init Progressbar
  pbmessage <- paste("  [", B.start, "->", B.start - n.merges, "blocks ]")
  pbwidth <- getOption("width") - nchar(pbmessage) - 4
  pb = txtProgressBar(min = 0, max = B.start * n.moves, width = pbwidth,
                      initial = 0, style = 3)
  cat(pbmessage)

  for (b in 1:B.start) {
    for (i in 1:n.moves) {

      if (exhaustive_search) {

        # Just go through each block and always accept.
        proposed_merge <- i

      } else {

        # Get a merge proposal using the mcmc sweep move proposal function.
        proposed_merge <- propose_move(V(block.graph)[b], block.graph, 1:B.start, block.graph.edges, eps)

        # If proposed merge is the current block, or was already proposed,
        # we don't need to waste time checking because decision will always
        # result in same state.
        if (b == proposed_merge) next
        if (any(merge.results[merge.results[,1] == b, 2] == proposed_merge)) next
        if (any(merge.results[merge.results[,1] == proposed_merge, 2] == b)) next

      }

      # Calculate entropy delta for merge and place into results
      proposed_partition <- replace(partition, which(partition == b), proposed_merge)
      proposed_entropy <- get_entropy(graph, proposed_partition)
      entropy_delta <- proposed_entropy - old_entropy

      # Place into results
      merge.idx <- (b-1)*n.moves + i
      merge.results[merge.idx,] <- c(b, proposed_merge, entropy_delta)

    }
    setTxtProgressBar(pb, b * i)
  }
  close(pb)

  # Clean up merge.results
  if(any(merge.results$dS == 0)) merge.results <- merge.results[-which(merge.results$dS == 0),]
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
    filter_g1 <- (merge_results_ordered$g1 == best_merge_pair$g1 | merge_results_ordered$g2 == best_merge_pair$g1)
    filter_g2 <- (merge_results_ordered$g1 == best_merge_pair$g2 | merge_results_ordered$g2 == best_merge_pair$g2)
    merge_results_ordered <- merge_results_ordered[-(filter_g1|filter_g2),]

    # Update new partition
    new_partition <- replace(new_partition, new_partition == best_merge_pair$g1, best_merge_pair$g2)

    n.merged <- n.merged + 1
  }

  # Apply moves
  new_partition <- check_partition(new_partition)
  new_entropy_delta <- get_entropy(graph, new_partition) - old_entropy

  # Perform mcmc sweeps to settle partition (costly)
  if (n.sweeps > 0) {
    sweep_results <- mcmc_sweep(graph, new_partition, B.start-1, n.sweeps, eps, beta)
    new_partition <- sweep_results$best_partition
    new_entropy_delta <- new_entropy_delta + sweep_results$best_entropy_delta
  }

  new_partition <- check_partition(new_partition)

  list("new_partition" = new_partition, "entropy_delta" = new_entropy_delta)
}
