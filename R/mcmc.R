mcmc_sweep <- function(graph, partition, n.blocks = NULL, eps = 1)
{
  if(is.null(n.blocks)) n.blocks <- max(partition)
  E <- block_edge_counts(graph, partition, n.blocks)
  e <- colSums(E)
  pnew <- sapply(V(graph), function(vertex) {
    mcmc_step(vertex, graph, partition, n.blocks, e, eps)
  })
  pnew
}


mcmc_step <- function(vertex, graph, partition, B, e, eps = 1) {
  vertex_neighborhood <- neighbors(graph, vertex)
  if(length(vertex_neighborhood) == 0) return(partition[vertex])  # if no neighbors, don't change.
  next_vertex <- resample(vertex_neighborhood, 1)
  t <- partition[next_vertex]
  if(runif(1) > Rt(t, e, B, eps)){
    next_vertex_neighborhood <- neighbors(graph, next_vertex)
    s <- partition[resample(next_vertex_neighborhood, 1)]
  } else {
    s <- sample(1:B, 1)
  }
  s
}



Rt <- function(t, e, B, eps) {
  eps * B / (e[t] + eps * B)
}
