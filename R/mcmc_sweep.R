#' Run multiple MCMC sweep overs nodes
#'
#' Runs a n.sweeps MCMC sweep across all nodes (for algorithm details see Piexoto,
#' 2018). Each node is given a chance to move blocks or stay in current block
#' and all nodes are processed in random order for each sweep.
#'
#' @param G An igraph graph.
#' @param p Vector of integer values giving the initial block membership
#' of each vertex
#' @param B Number of blocks.
#' @param dc Wether to use degree correction.
#' @param n.sweeps Number of sweeps to run.
#' @param eps (optional) A number giving the ...
#' @param beta (optional) A number giving the greediness of the moves.
#' @return A new partition given as a vector of integer values.
#' @export
#' @import igraph

mcmc_sweep <- function(G, p, B, dc = FALSE, n.sweeps = 1, eps = 0.1, beta = 1)
{
  # Initial checks
  stopifnot(is.igraph(G))
  stopifnot(all(degree(G) > 0))

  B <- as.integer(B)
  p <- as.integer(p)
  ds <- as.logical(dc)

  old_entropy <- get_entropy(G, p, dc)

  # Bookkeper variables
  entropy_delta <- rep(0, n.sweeps + 1)
  new_partition <- p
  best_partition <- p
  best_entropy_delta <- 0
  best_entropy <- old_entropy

  for (i in 1:n.sweeps) {
    sweep_results <- mcmc_single_sweep(G, new_partition, B, dc, eps, beta)
    entropy_delta[i+1] <- entropy_delta[i] + sweep_results$entropy_delta
    new_partition <- sweep_results$new_partition
    if (sweep_results$best_entropy < best_entropy) {
      best_partition <- sweep_results$best_partition
      best_entropy_delta <- sweep_results$best_entropy_delta
      best_entropy <- sweep_results$best_entropy
    }
  }

  list("best_partition" = best_partition, "best_entropy_delta" = best_entropy_delta,
       "best_entropy" = best_entropy, "new_partiton" = new_partition,
       "entropy_delta" = entropy_delta)
}



#' Run a single MCMC sweep over nodes
#'
#' Runs a single MCMC sweep across all nodes (for algorithm details see Piexoto,
#' 2018). Each node is given a chance to move blocks or stay in current block
#' and all nodes are processed in random order.
#'
#' @param G An igraph graph.
#' @param p Vector of integer values giving the block membership of each
#' vertex
#' @param B Number of blocks.
#' @param dc Wether to use degree correction.
#' @param eps (optional) A number giving the ...
#' @param beta (optional) A number giving the greediness of the moves.
#' @return A new partition given as a vector of integer values.
#' @export
#' @import igraph

mcmc_single_sweep <- function(G, p, B, dc = FALSE, eps = 0.1, beta = 1)
{
  # Initial checks
  stopifnot(is.igraph(G))
  stopifnot(all(degree(G) > 0))

  B <- as.integer(B)
  p <- as.integer(p)
  dc <- as.logical(dc)

  # Book keeper variables for this sweeps stats
  entropy_delta <- 0
  old_entropy <- get_entropy(G, p, dc)

  # Book keeper variable for best state
  original_entropy <- old_entropy
  best_entropy_delta <- 0
  best_entropy <- old_entropy
  best_partition <- p

  # Keep a list of block adjacent edges
  block.edges <- block_edge_list(G, p, B)

  for(curr_v in sample(V(G))){

    # Get a move proposal
    proposed_new_block = propose_move(curr_v, G, p, block.edges, eps)
    old_block = p[curr_v]

    # If proposed block is the current block, we don't need to waste
    # time checking because decision will always result in same state.
    if (old_block == proposed_new_block) next

    # Calculate acceptance probability based on posterior changes
    proposal_results = get_proposal_results(curr_v, proposed_new_block,
                                            old_entropy, G, p, block.edges,
                                            eps, beta, dc = dc)

    # Make movement decision
    move_accepted = proposal_results$prob_to_accept > stats::runif(1)

    if (move_accepted) {
      p <- swap_blocks(p, curr_v, proposed_new_block)
      entropy_delta <- entropy_delta + proposal_results$entropy_delta
      old_entropy <- proposal_results$new_entropy
      block.edges <- update_block_edge_list(block.edges, curr_v, G, proposed_new_block, old_block)

      if (old_entropy < best_entropy) {
        best_partition <- p
        best_entropy_delta <- entropy_delta
        best_entropy <- old_entropy
      }
    }
  }
  list("new_partition" = p, "entropy_delta" = entropy_delta, "new_entropy" = old_entropy,
       "best_partition" = best_partition, "best_entropy_delta" = best_entropy_delta, "best_entropy" = best_entropy)
}



#' A single move proposal during an MCMC sweep
#'
#' @param curr_v A node in graph.
#' @param G An igraph graph.
#' @param p Vector of integer values giving the block membership of each
#' vertex
#' @param block.edges a list of all incident edge ids per block
#' @param eps (optional) A number giving the ...
#' @return A new group membership for vertex.
#' @import igraph

propose_move <- function(curr_v, G, p, block.edges, eps = 1) {

  # Get block membership of a random neighbor
  vertex_neighborhood <- neighbors(G, curr_v)
  next_vertex <- resample(vertex_neighborhood, 1)
  neighbor_block <- p[next_vertex]

  # Decide if we are going to choose a random block for our node
  ergo_amnt = eps * length(block.edges)
  draw_from_neighbor = stats::runif(1) > ergo_amnt / (length(block.edges[[neighbor_block]]) + ergo_amnt)

  # Draw from potential candidates
  if (draw_from_neighbor) {
    # Sample a random edge from all edges connected to neighbor block
    random_edge_to_neighbor_block <- sample(block.edges[[neighbor_block]], 1)
    g1 <- p[head_of(G, random_edge_to_neighbor_block)]
    g2 <- p[tail_of(G, random_edge_to_neighbor_block)]
    move_proposal <- ifelse(neighbor_block != g1, g1, g2)
  } else {
    # Sample a random block from all blocks
    move_proposal <- sample(1:length(block.edges), 1)
  }
  move_proposal
}



#' Calculate entropy delta of the SBM before and after the proposed move and
#' the ratio of the probabilities of moving to the proposed block before the move and
#' moving back to the original block after the move.
#'
#' @param curr_v A node in graph.
#' @param proposed_new_block An integer giving the proposed block id.
#' @param old_entropy A number giving the pre move entropy
#' @param G An igraph graph.
#' @param old_partition Vector of integer values giving the block membership of each
#' vertex pre move.
#' @param block.edges a list of all incident edge ids per block
#' @param eps (optional) A number giving the ...
#' @param beta (optional) A number giving the greediness of the moves.
#' @param dc Wether to use degree correction.
#' @return A new group membership for vertex.
#' @import igraph

get_proposal_results <- function(curr_v, proposed_new_block,
                                 old_entropy, G, old_partition, block.edges,
                                 eps = 0.1, beta = 1, dc = FALSE){

  old_block <- old_partition[curr_v]

  # No need to go on if we're "swapping" to the same group
  if (proposed_new_block == old_block) return(list("entropy_delta" = 0,
                                                   "prob_to_accept" = 1,
                                                   "new_entropy" = old_entropy))

  # Calcualte old entropy (speed this up!)
  old_edge_matrix <- block_edge_counts(G, old_partition, n.blocks = length(block.edges))
  old_node_counts <- block_node_counts(old_partition, n.blocks = length(block.edges))

  # Get new partition
  new_partition <- replace(old_partition, curr_v, proposed_new_block)

  # Calculate new entropy (speed this up!)
  new_edge_matrix <- block_edge_counts(G, new_partition, n.blocks = length(block.edges))
  if(dc == FALSE) {
    new_node_counts <- block_node_counts(new_partition, n.blocks = length(block.edges))
    new_entropy <- entropy_trad(new_edge_matrix, new_node_counts, is.directed(G), is.simple(G))
  } else {
    # TODO: Fix degree!
    if (!is.directed(G)) {
      d <- list("total" = degree(G))
    } else {
      d <- list("in" = degree(G, mode = "in"), "out" = degree(G, mode = "out"))
    }
    new_entropy <- entropy_corrected(new_edge_matrix, d, directed = is.directed(G))
  }

  entropy_delta <- new_entropy - old_entropy

  # Calculate R_t (for all t)
  new_block_edges <- update_block_edge_list(block.edges, curr_v, G, proposed_new_block, old_block)
  new_edge_counts <- colSums(new_edge_matrix)
  new_Rt <- eps * length(block.edges) / (new_edge_counts + eps * length(block.edges))

  old_edge_counts <- colSums(old_edge_matrix)
  old_Rt <- eps * length(block.edges) / (old_edge_counts + eps * length(block.edges))

  # Calculate fraction of neighbors belonging to each block
  vertex_neighborhood <- neighbors(G, curr_v)
  vertex_neighborhood_block_memberships <- old_partition[vertex_neighborhood]
  vertex_neighborhood_blocks <- sort(unique(vertex_neighborhood_block_memberships))
  vertex_neighborhood_size <- length(vertex_neighborhood)

  prob_return_to_old <- 0
  prob_move_to_new <- 0

  for(t in vertex_neighborhood_blocks){
    pti <- sum(vertex_neighborhood_block_memberships == t)/vertex_neighborhood_size

    # Calculated after move
    if( new_edge_counts[t] != 0 ) {
      prob_return_to_old <- prob_return_to_old + pti *
        ( (1-new_Rt[t]) * new_edge_matrix[t,old_block] / new_edge_counts[t]  + new_Rt[t]/length(new_edge_counts))
    }

    # Calculated before move
    if( old_edge_counts[t] != 0 ) {
      prob_move_to_new <- prob_move_to_new + pti *
        ( (1-old_Rt[t]) * old_edge_matrix[t,proposed_new_block] / old_edge_counts[t]  + old_Rt[t]/length(old_edge_counts))
    }
  }

  prob_ratio <- prob_return_to_old/prob_move_to_new
  prob_to_accept <- min(exp(- beta * entropy_delta) * prob_ratio, 1)

  return(list("entropy_delta" = entropy_delta, "prob_to_accept" = prob_to_accept, "new_entropy" = new_entropy))
}




#' Create list of edges per block
#'
#' @param G An igraph graph.
#' @param p Vector of integer values giving the block membership of each
#' vertex
#' @param B Number of blocks.

block_edge_list <- function(G, p, B = NULL)
{
  if(is.null(B)) {
    blocks <- order(unique(p))
  } else {
    blocks <- 1:B
  }
  lapply(blocks, function(b){
    b.vertices <- V(G)[which(p == b)]
    b.edge.ids <- lapply(incident_edges(G, b.vertices, mode = "all"), as_ids)
    b.edge.ids <- unlist(b.edge.ids)
    sort(unique(b.edge.ids))
  })
}


#' Create list of edges per block
#'
#' @param block.edges ...
#' @param v vertex id
#' @param G An igraph graph.
#' @param new.b ...
#' @param old.b ...
#' @return A new block.edges list.

update_block_edge_list <- function(block.edges, v, G, new.b, old.b)
{
  v.edges <- as_ids(incident(G, v, mode = "all"))
  block.edges[[old.b]] <- setdiff(block.edges[[old.b]], v.edges)
  block.edges[[new.b]] <- union(block.edges[[new.b]], v.edges)
  block.edges
}


#' Create list of edges per block
#'
#' @param p Vector of integer values giving the block membership of each
#' vertex
#' @param curr_v vertex id
#' @param proposed_new_block block id

swap_blocks <- function(p, curr_v, proposed_new_block){
  p[curr_v] <- proposed_new_block
  p
}
