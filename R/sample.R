#' Sample planted partition model
#'
#' Sampling from the planted partition model
#'
#' This function samples graphs from a stochastic block model by (doing the
#' equivalent of) \code{n.trial} Bernoulli trials for each potential edge with
#' probabilities given by \code{p} for in-group edges and \code{q} for between
#' group edges.
#'
#' @param n Number of vertices in the graph.
#' @param p Probability of creating an edge between vertices of the same group.
#' @param q Probability of creating an edge between vertices of different groups.
#' @param block.sizes Numeric vector giving the number of vertices in each group.
#' The sum of the vector must match the number of vertices.
#' @param n.trials Number of repeated Bernoulli trials per edge. If
#' \code{n.trials} is greater than 1, generates a weighted graph.
#' @param directed Logical scalar, whether to generate a directed graph.
#' @param loops Logical scalar, whether self-loops are allowed in the graph.
#' @return An igraph graph.
#' @keywords graphs, sample, planted partition
#' @examples
#' ## Three groups with only a few connection between groups
#' g1 <- sample_ppm(1000, p=0.3, q=0.01, block.sizes=c(100,600,300))
#' g1
#' ## Three groups with weighted connections.
#' g2 <- sample_ppm(1000, p=0.1, q=0.003, block.sizes=c(100,600,300), n.trials=3)
#' g2
#' @export
#' @import igraph

sample_ppm <- function (n, p, q, block.sizes,
                        n.trials = 1, directed = FALSE, loops = FALSE)
{
  n <- as.integer(n)
  p <- as.double(p)
  q <- as.double(q)
  block.sizes <- as.integer(block.sizes)
  n.trials <- as.integer(n.trials)
  directed <- as.logical(directed)
  loops <- as.logical(loops)
  pm <- diag(p, length(block.sizes))
  pm[upper.tri(pm) | lower.tri(pm)] <- q
  res <- sample_sbm(n, pm, block.sizes, directed, loops)
  if(n.trials > 1){
    # See https://stackoverflow.com/a/27766128 (Combine two graphs with weights)
    E(res)$weight <- 1
    for (i in seq(n.trials - 1)){
      resnew <- sample_sbm(n, pm, block.sizes, directed, loops)
      E(resnew)$weight <- 1
      res <- union(res, resnew)
      E(res)$weight_1[is.na(E(res)$weight_1)] <- 0
      E(res)$weight_2[is.na(E(res)$weight_2)] <- 0
      E(res)$weight <- E(res)$weight_1 + E(res)$weight_2
      res <- remove.edge.attribute(res, "weight_1")
      res <- remove.edge.attribute(res, "weight_2")
      if (igraph_opt("add.params")) {
        res$name <- res$name_1
        res <- remove.graph.attribute(res, "name_1")
        res <- remove.graph.attribute(res, "name_2")
        res$loops <- res$loops_1
        res <- remove.graph.attribute(res, "loops_1")
        res <- remove.graph.attribute(res, "loops_2")
      }
    }
  }
  if (igraph_opt("add.params")) {
    res$name <- "Planted partition model"
    if(n.trials > 1) res$name <- "Weighted planted partition model"
    res$loops <- loops
  }
  res
}


#' Sample planted partition model (Piexoto 2014)
#'
#' Alternative sampling from the planted partition model (see Piexoto 2014)
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
#' G <- sample_ppm2(10000, c=0.99, k=10, B=3)
#' p <- c(rep(1,3334), rep(2, 3333), rep(3, 3333))
#' plot(G, vertex.label=NA, vertex.color=p)
#' @export
#' @import igraph

sample_ppm2 <- function (N, c, k, B, directed = FALSE, loops = FALSE)
{
  # Input checks
  N <- as.integer(N)
  c <- as.double(c)
  k <- as.integer(k)
  B <- as.integer(B)
  directed <- as.logical(directed)
  stopifnot(!directed)
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
    res$c <- c
    res$k <- k
    res$B <- B
  }

  # Return
  res
}
