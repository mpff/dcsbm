#' Sample planted partition model (Piexoto 2014)
#'
#' Sampling from the planted partition model as in Piexoto 2014.
#'
#' This function samples graphs from a stochastic block model by building a
#' block adjacency matrix from the Parameters N, k, B and c.
#'
#' @param N Number of vertices in the graph.
#' @param c A number between 0 and 1 controlling the
#' @param k The average degree of a vertex.
#' @param B The number of blocks.
#' @param directed Logical scalar, whether to generate a directed graph.
#' @param loops Logical scalar, whether self-loops are allowed in the graph.
#' @return An igraph graph.
#' @keywords graphs, sample, planted partition
#' @examples
#' ## Three groups with only a few connection between groups
#' G <- sample_ppm(300, c=0.9, k=10, B=3)
#' p <- c(rep(1,100), rep(2, 200), rep(3, 300))
#' plot(G, vertex.label=NA, vertex.color=p)
#' @export
#' @import igraph

sample_ppm <- function (N, c, k, B, directed = FALSE, loops = FALSE)
{
  # Input checks
  N <- as.integer(N)
  c <- as.double(c)
  k <- as.integer(k)
  B <- as.integer(B)
  directed <- as.logical(directed)
  loops <- as.logical(loops)

  # Create block.sizes
  n <- floor(N / B)
  block.sizes <- rep(n, B)
  while(sum(block.sizes) < N) block.sizes[1] <- block.sizes[1] + 1  # pad group 1

  # Create pref.matrix
  E <- N * k  # avg. no of edges
  pm <- diag(2 * E * c / B, B)
  pm[upper.tri(pm) | lower.tri(pm)] <- 2 * E * (1-c) / (B * (B-1))

  # Sample from SBM according to pm
  res <- make_empty_graph(N, directed = directed)
  for(b in 1:B){
    b.first <- sum(block.sizes[0:(b-1)]) + 1
    b.last  <- sum(block.sizes[0:b])
    b.vertices <- V(res)[b.first:b.last]
    for(b.to in 1:B){
      if(!directed & b < b.to) next()
      b.to.first <- sum(block.sizes[0:(b.to-1)]) + 1
      b.to.last  <- sum(block.sizes[0:b.to])
      b.to.vertices <- V(res)[b.to.first:b.to.last]
      if(b != b.to){
        vertices.fro <- sample(b.vertices, pm[b,b.to], replace = TRUE)
        vertices.to <- sample(b.to.vertices, pm[b,b.to], replace = TRUE)
      } else {
        vertices.fro <- sample(b.vertices, pm[b,b.to]/2, replace = TRUE)
        vertices.to <- sample(b.to.vertices, pm[b,b.to]/2, replace = TRUE)
      }
      res <- add_edges(res, c(rbind(vertices.fro, vertices.to)))
      if(!loops & any(is.loop(res))){
        lp.ends <- ends(res, which(is.loop(res)))
        lp.ends[,2] <- sapply(lp.ends[,1], function(v){
          sample(b.vertices[-v], 1)
        })
        res <- delete_edges(res, which(is.loop(res)))
        lp.ends.fix <- c(rbind(lp.ends[,1], lp.ends[,2]))
        res <- add_edges(res, lp.ends.fix)
      }
    }
  }

  # Add parameters
  if (igraph_opt("add.params")) {
    res$name <- "Planted partition model"
    res$loops <- loops
  }

  # Return
  res
}


#' Sample planted partition model with degree variability (Piexoto 2020)
#'
#' This function samples graphs from a stochastic block model by building a
#' block adjacency matrix from the Parameters N, k, B and c and accounting
#' for a exponential degree sequence inside each block controlled by \code{k_coef}.
#'
#' @param N Number of vertices in the graph.
#' @param c A number between 0 and 1 controlling the
#' @param k The average degree of a vertex.
#' @param B The number of blocks.
#' @param k_coef A number describing the degree variability within each block.
#' A value of 0 means that all vertices have the same expected degree.#'
#' @param directed Logical scalar, whether to generate a directed graph.
#' @param loops Logical scalar, whether self-loops are allowed in the graph.
#' @return An igraph graph.
#' @keywords graphs, sample, planted partition
#' @examples
#' ## Three groups with only a few connection between groups
#' G <- sample_dcppm(300, c=0.9, k=10, B=3, k_coef=2)
#' p <- c(rep(1,100), rep(2, 200), rep(3, 300))
#' plot(G, vertex.label=NA, vertex.color=p)
#' @export
#' @import igraph

sample_dcppm <- function (N, c, k, B, k_coef = 0, directed = FALSE, loops = FALSE)
{
  # Input checks
  N <- as.integer(N)
  c <- as.double(c)
  k <- as.integer(k)
  B <- as.integer(B)
  k_coef <- as.double(k_coef)
  directed <- as.logical(directed)
  loops <- as.logical(loops)

  # Create block.sizes
  n <- floor(N / B)
  block.sizes <- rep(n, B)
  while(sum(block.sizes) < N) block.sizes[1] <- block.sizes[1] + 1  # pad group 1

  # Create edge count and pref. matrices
  E <- N * k  # avg. no of edges
  em <- diag(2 * E * c / B, B)
  em[upper.tri(em) | lower.tri(em)] <- 2 * E * (1-c) / (B * (B-1))
  pm <- 2*em/sum(em)

  # Sample a degree sequence
  deg <- rep(1, N)  # at least degree 1 for each vertex

  if (!directed) {
    bdeg <- colSums(em/2)
  } else {
    bdeg <- colSums(em)
  }

  for (b in 1:B) {
    # Get vertices inbn block b
    b.first <- sum(block.sizes[0:(b-1)]) + 1
    b.last  <- sum(block.sizes[0:b])

    # build variability distribution over vertices in block b
    prob <- (1:(b.last - b.first + 1))**k_coef
    prob <- prob/sum(prob)

    # Sample degrees from variability distribution
    while (sum(deg[b.first:b.last]) < bdeg[b]) {
      i <- sample(b.first:b.last, 1, prob = prob)
      deg[i] <- deg[i] + 1
    }
  }

  # Sample from SBM according to em, pm and deg.
  res <- make_empty_graph(N, directed = directed)

  for(b in 1:B){

    # Get vertices in block b
    b.first <- sum(block.sizes[0:(b-1)]) + 1
    b.last  <- sum(block.sizes[0:b])
    b.vertices <- V(res)[b.first:b.last]

    # Iterate over each vertex
    for(v.fro in b.vertices) {
      for (i in 1:deg[v.fro]) {
        # Sample adjacent block according to pm
        b.to <- sample(1:B, 1, prob = pm[b,])

        # Sample a vertex from this block according to deg
        b.to.first <- sum(block.sizes[0:(b.to-1)]) + 1
        b.to.last  <- sum(block.sizes[0:b.to])
        b.to.vertices <- V(res)[b.to.first:b.to.last]
        v.to <- sample(b.to.vertices, 1, prob = deg[b.to.first:b.to.last])

        ## TODO: Check if v.to = v and account for loops!
        if(!loops) {
          while(v.to == v.fro) v.to <- sample(b.to.vertices, 1, prob = deg[b.to.first:b.to.last])
        }

        # Add edge to graph
        res <- add_edges(res, c(v.fro, v.to))
      }
    }
  }

  # Add parameters
  if (igraph_opt("add.params")) {
    res$name <- "Planted partition model with degree variability"
    res$loops <- loops
  }

  # Return
  res
}
