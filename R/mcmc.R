#' Run a single MCMC sweep over nodes
#'
#' Runs a single MCMC sweep across all nodes (for algorithm details see Piexoto,
#' 2018). Each node is given a chance to move blocks or stay in current block
#' and all nodes are processed in random order.
#'
#' @param graph An igraph graph.
#' @param partition Vector of integer values giving the block membership of each
#' vertex
#' @param eps (optional) A number giving the ...
#' @return A new partition given as a vector of integer values.
#' @export
#' @import igraph

mcmc_sweep <- function(graph, partition, B = NULL, eps = 1, beta = 1)
{
  # Initial checks
  stopifnot(is.igraph(graph))
  partition <- as.integer(partition)
  # Remove vertices with no connections, Rename partition
  zero.vertices <- which(degree(graph) == 0)
  # Renumber blocks in partition.
  partition <- check_partition(partition)
  # Final checks
  stopifnot(length(graph) == length(partition))
  if(length(zero.vertices) == length(partition)) return(rep.int(1, length(partition)))
  # Get parameters
  E <- block_edge_counts(graph, partition)
  n <- block_node_counts(partition)
  PE <- edge_probs(E)
  Pr <- random_probs(E, eps)
  if(!is.null(B)){
    B <- as.integer(B)
  } else {
    B <- length(Pr)
  }
  block.edges <- block_edge_list(G, partition)
  # Get sweep order
  v.random <- sample(V(graph))
  if(length(zero.vertices) > 0) v.random <- v.random[-which(v.random %in% zero.vertices)]
  S.start <- entropy_undirected_trad(E, n)
  for(vertex in v.random){
    #move.proposal <- mcmc_step(vertex, graph, partition, PE, Pr, B, eps)
    move.proposal <- mcmc_step3(vertex, graph, partition, block.edges, Pr, B, eps)
    if(move.proposal == partition[vertex]) next()
    S.new <-
    block.edges <- update_block_edge_list(block.edges, vertex, graph, move.proposal, partition[vertex])

    # Update partition
    partition[vertex] <- move.proposal
    partition <- check_partition(partition)
    #E <- block_edge_counts(graph, partition)
    #PE <- edge_probs(E)
    #Pr <- random_probs(E, eps)
    Pr <- eps * B / (sapply(block.edges, length) + eps * B)
  }
  #new.partition <- sapply(v.random, function(vertex) {
  #  mcmc_step(vertex, graph, partition, PE, Pr, B, eps)
  #})
  #new.partition <- new.partition[order(v.random)]
  if(length(zero.vertices) > 0) {
    #partition[-zero.vertices] <- new.partition
    partition[zero.vertices] <- 1  # sort zero nodes into group 1.
  } #else {
  #  partition <- new.partition
  #}
  partition
}


#' A single MCMC step over one node
#'
#' @param vertex A node in graph.
#' @param graph An igraph graph.
#' @param partition Vector of integer values giving the block membership of each
#' vertex
#' @param PE asdf
#' @param Pr asdf
#' @param B Number of groups to sample from
#' @param eps (optional) A number giving the ...
#' @return A new group membership for vertex.
#' @import igraph

mcmc_step <- function(vertex, graph, partition, PE, Pr, B, eps = 1) {
  vertex_neighborhood <- neighbors(graph, vertex)
  next_vertex <- resample(vertex_neighborhood, 1)
  s <- partition[next_vertex]
  if(stats::runif(1) > Pr[s]){
    next_vertex_neighborhood <- neighbors(graph, next_vertex)
    t <- partition[resample(next_vertex_neighborhood, 1)]
  } else {
    t <- sample(1:B, 1)
  }
  t
}

#' A single MCMC step over one node (alternative)
#'
#' @param vertex A node in graph.
#' @param graph An igraph graph.
#' @param partition Vector of integer values giving the block membership of each
#' vertex
#' @param PE asdf
#' @param Pr asdf
#' @param B Number of groups to sample from
#' @param eps (optional) A number giving the ...
#' @return A new group membership for vertex.
#' @import igraph

mcmc_step2 <- function(vertex, graph, partition, PE, Pr, B, eps = 1) {
  vertex_neighborhood <- neighbors(graph, vertex)
  next_vertex <- resample(vertex_neighborhood, 1)
  s <- partition[next_vertex]
  if(stats::runif(1) > Pr[s]){
    t <- sample(1:length(Pr), 1, prob = PE[,s])
  } else {
    t <- sample(1:B, 1)
  }
  t
}

#' A single MCMC step over one node (alternative 3)
#'
#' @param vertex A node in graph.
#' @param graph An igraph graph.
#' @param partition Vector of integer values giving the block membership of each
#' vertex
#' @param block.edges asdf
#' @param Pr asdf
#' @param B Number of groups to sample from
#' @param eps (optional) A number giving the ...
#' @return A new group membership for vertex.
#' @import igraph

mcmc_step3 <- function(vertex, graph, partition, block.edges, Pr, B, eps = 1) {
  vertex_neighborhood <- neighbors(graph, vertex)
  next_vertex <- resample(vertex_neighborhood, 1)
  s <- partition[next_vertex]
  if(stats::runif(1) > Pr[s]){
    random.edge.s <- sample(block.edges[[s]], 1)
    g1 <- partition[head_of(G, random.edge.s)]
    g2 <- partition[tail_of(G, random.edge.s)]
    t <- ifelse(s != g1, g1, g2)
  } else {
    t <- sample(1:B, 1)
  }
  t
}



#' The probability of sampling a fully random (or new) group in a single MCMC step
#'
#' @param E Block adjacency matrix
#' @param eps (optional) A number giving the ...
#' @return A vector of length B giving the above probability per group.
#' @import igraph

random_probs <- function(E, eps = 1) {
  e <- colSums(E)
  eps * (nrow(E)+1) / (e + eps * (nrow(E)+1))
}


#' The probability of ...
#'
#' @param E Block adjacency matrix
#' @return A B x B matrix with rows giving the edge probabilities per group.

edge_probs <- function(E){
  e <- colSums(E)
  E %*% diag(1/e, nrow = nrow(E))
}

