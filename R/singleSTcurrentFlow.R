#' a single source-sink current flow example for an undirected graph from an online tutorial
#' 
#' The tutorial behind this implementation of current flow has an example for current flow for a single source-target. 
#' It is included primarily for tests.
#' 
#' @docType data
#' @usage data(singleSTcurrentFlow)
#'
#' @format An undirected igraph object with 8 nodes and 10 edges
#' \describe{
#'   \item{name}{A unique node identifier}
#'   \item{isSource}{Boolean - details if the node is the (single) source or unit current}
#'   \item{isSink}{Boolean - details if the node is the (single) sink or unit current}
#'   \item{absoluteFlow}{A node attribute detailing how much current flows through the node}
#'   \item{potential}{A node attribute detailign the electrical potential at each node for the given source/sink setup}
#'   \item{absoluteCurrent}{An edge attribute detailing how much current flows along the edge}
#' }
#' @source \url{http://www.sci.unich.it/~francesc/teaching/network/flowcentrality.html}
"singleSTcurrentFlow"