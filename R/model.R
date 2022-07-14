#' Degree corrected stochastic block model
#'
#' Estimate a ...
#'
#' @param graph An igraph graph.
#' @param n.blocks Number of blocks.
#' @param degree_correction Type of degree correction to use. "oneway" for one
#' parameter per vertex, "twoways" for two parameters per vertex (input/output),
#' "none" for no degree correction.
#' @param directed Logical scalar, Whether to use the directed or undirected SBM.
#' Ignored for undirected graphs.
#' @param control List of parameters for the inference algorithm.
#' @return Todo
#' @keywords graphs, inference, stochastic block model, degree correction
#' @examples
#' ## Three groups with weighted connections.
#' g1 <- sample_ppm(1000, p=0.1, q=0.003, block.sizes=c(100,600,300))
#' model <- sbm(g1, 3)
#' @export
#' @import igraph

sbm <- function (graph, n.blocks,
                   degree_correction = c("none", "oneway", "twoways"),
                   directed = FALSE, control = list())
{
  stopifnot(is.igraph(graph))
  n.blocks <- as.integer(n.blocks)
  degree_correction <- match.arg(degree_correction)
  directed <- as.logical(directed)
  control <- NULL
  NULL
}
