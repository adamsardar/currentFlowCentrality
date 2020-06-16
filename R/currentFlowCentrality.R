#' currentFlowCentrality: A package for a betweeness cetrality based on current flow
#'
#' Most centralities are based on the idea of counting shortest paths from A to B and then iterating over all AB pairs.
#' Newman et al suggests the idea of using the principle of electrical current flow, hence counting all paths.
#'
#' @docType package
#' @name currentFlowCentrality
#' @author Adam Sardar
#' @references Newman, M. E. J. (2005). A measure of betweenness centrality based on random walks. Social Networks, 27(1)
#' @import methods
#' @importFrom Rcpp evalCpp
#' @importFrom magrittr "%<>%"
#' @importFrom magrittr add
#' @useDynLib currentFlowCentrality
NULL