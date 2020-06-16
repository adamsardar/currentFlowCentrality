#' Calculate the current flow distance and betweeness, as popularised by Newman et al.
#' 
#' This is a nice little wrapper function that delegates out computation to c++ implementations
#'
#' @param network The network for wich to calculate current flow centrality for. If directed then it will be made undirected with a warning.
#' @import igraph
#' @export
#' @return list A list of two fatures: betweeness and distance 
#' @references Newman, M. E. J. (2005). A measure of betweenness centrality based on random walks. Social Networks, 27(1)
calculateCurrentFlow <- function(network){
  
  if(! is.connected(network)){stop("Input network MUST be connected! Try taking the largest connected component?")}
  
  adjMat <- network %>% get.adjacency(sparse=TRUE)
  Glapl <- network %>% graph.laplacian(sparse=TRUE)

  Gpinv <- dcPinv(Glapl)
  
  betweenessVector <- calculateAverageCurrentFlowBetweeness(Gpinv, adjMat)
  names(betweenessVector) <- colnames(Glapl) 
  
  distanceVector <- calculateCurrentFlowDistance(Gpinv)
  names(distanceVector) <- colnames(Glapl) 
  
  returnList <- list(
                      betweeness = betweenessVector,
                      distance = distanceVector
  )
  
  return(returnList)
}