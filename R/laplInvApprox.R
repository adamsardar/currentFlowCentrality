is.whole.number <- function(x){ x == as.integer(x) }
is.positive.semidefinite <- function(x){ all(x >= 0) }
is.laplacian <- function(x){ all( rowSums(x) == 0) & 
    all(colSums(x) == 0) & 
    isSymmetric(x)}

globalVariables(c("."))
#' Approximate the psuedo-inverse of a laplacian matrix from a handful of eigenvectors and values
#'
#' Based on work by Bozzo & Franceschet who showed that one could closely approximate the psuedo-inverse of a graph laplacian using a small (k << n)
#' number of eigenvalues and vectors. Consult their paper for details
#'
#' @family approximatePsuedoInverse
#'
#' @param laplMat A laplacian matrix for which to calculate the approximate psudeo-inverse. Since we only handle undirected graphs, this must be symmetric
#' @param kApprox Number of eigenvectors/values to use in construction of the psuedo-inverse
#' @param exact Logical flag to turn off use of rARPACK in generation of the desired few eigenvectors/values. This is very useful for testing
#'
#' @return A list of two entries: ginv - a matrix of the approximate graph laplacian approximated using the first k eigenvectors and twoNormError,
#' the spectral norm error in the matrix estimate
#' @import Matrix
#' @importFrom ensurer ensure
#' @importFrom magrittr multiply_by_matrix
#' @importFrom magrittr subtract
#' @importFrom RSpectra eigs_sym
#' @export
#'
#' @examples
#' data(dolphinPodLusseau)
#' 
#' dolphnLapl <- igraph::graph.laplacian(dolphinPodLusseau)
#' 
#' str(kthCutoffGinvApproximation(dolphnLapl))
#' 
#' @references Bozzo, E., & Franceschet, M. (2012). Effective and efficient approximations of the generalized inverse of the graph Laplacian matrix with an application to current-flow betweenness centrality. Physics and Society
kthCutoffGinvApproximation <- function(laplMat,kApprox = 8, exact = FALSE){
  
  
  laplMat %<>% ensure(is.laplacian,
                      err_desc = "Expecting laplMat to be a symetric matrix and for all rows and columns to sum to 0")
  
  kApprox %<>% ensure(is.numeric,
                      is.whole.number,
                      is.positive.semidefinite,
                      . < nrow(laplMat),
                      err_desc = "kApprox must be a semi-definite whole number less than the number of rows/columns in laplMat")
  
  exact %<>% ensure(is.logical,
                    err_desc = "exact shoudl be a logical flag")
  
  
  if(exact){
    
    fewEigenvecs <- eigen(laplMat, symmetric = TRUE)
    fewEigenvecs$values <- fewEigenvecs$values[( nrow(laplMat) - kApprox ):nrow(laplMat)]
    fewEigenvecs$vectors <- fewEigenvecs$vectors[,( nrow(laplMat) - kApprox):nrow(laplMat)]
  }else{
    fewEigenvecs <- eigs_sym(laplMat,kApprox+1,which="SM",opts = list(maxitr=2000))    
  }
  
  ### K'th cutoff - sum_j (1/lambda_j) * V(:,j)V(:.j)^t
  approx.graphLapl.ginv <- fewEigenvecs$vectors[,kApprox:2] %>%
    multiply_by_matrix(diag(x=1/fewEigenvecs$values[kApprox:2] , nrow = length(kApprox:2))) %>%   
    tcrossprod(y=fewEigenvecs$vectors[,kApprox:2])
  
  if(! is.null(row.names(laplMat))){row.names(approx.graphLapl.ginv) <- row.names(laplMat)  }
  if(! is.null(colnames(laplMat))){colnames(approx.graphLapl.ginv) <- colnames(laplMat)  }
  
  res <- list(ginv = approx.graphLapl.ginv,
              twoNormError = fewEigenvecs$values[kApprox]/fewEigenvecs$values[1])
  
  class(res) <- "kthCutoffGinvApproximation"
  
  return(res)
}

globalVariables(c("."))
#' Approximate the psuedo-inverse of a laplacian matrix from a handful of eigenvectors and values and a cunning guess for the remaining values
#'
#' Based on work by Bozzo & Franceschet who showed that one could closely approximate the psuedo-inverse of a graph laplacian using a small (k << n)
#' number of eigenvalues and vectors alongside an educated guess at the harmonic mean of the remainign values (called the stretch approximation). Consult
#' their paper for details
#' 
#' @family approximatePsuedoInverse
#' 
#' @inheritParams kthCutoffGinvApproximation
#'
#' @return A list of two entries: ginv - a matrix of the approximate graph laplacian approximated using the first k eigenvectors and twoNormError,
#' the spectral norm error in the matrix estimate
#' @import Matrix
#' @importFrom ensurer ensure
#' @importFrom magrittr multiply_by_matrix
#' @importFrom magrittr subtract
#' @importFrom RSpectra eigs_sym
#' @export
#'
#' @examples
#' data(dolphinPodLusseau)
#' 
#' dolphnLapl <- igraph::graph.laplacian(dolphinPodLusseau)
#' 
#' str(stretchGinvApproximation(dolphnLapl))
#' 
#' @references Bozzo, E., & Franceschet, M. (2012). Effective and efficient approximations of the generalized inverse of the graph Laplacian matrix with an application to current-flow betweenness centrality. Physics and Society
stretchGinvApproximation <- function(laplMat,kApprox = 8, exact = FALSE){
  
  laplMat %<>% ensure(is.laplacian,
                      err_desc = "Expecting laplMat to be a symetric matrix and for all rows and columns to sum to 0")
  
  kApprox %<>% ensure(is.numeric,
                      is.whole.number,
                      is.positive.semidefinite,
                      . < nrow(laplMat),
                      err_desc = "kApprox must be a semi-definite whole number less than the number of rows/columns in laplMat")
  
  exact %<>% ensure(is.logical,
                    err_desc = "exact shoudl be a logical flag")
  
  if(exact){
    
    fewEigenvecs <- eigen(laplMat)
    largestEV <- fewEigenvecs
    largestEV$values <- largestEV$values[1]
    fewEigenvecs$values <- fewEigenvecs$values[( nrow(laplMat) - kApprox):nrow(laplMat)]
    fewEigenvecs$vectors <- fewEigenvecs$vectors[,( nrow(laplMat) - kApprox):nrow(laplMat)]
  }else{
    
    fewEigenvecs <- eigs_sym(laplMat,kApprox+1,which="SM",opts = list(maxitr=2000))
    largestEV <- eigs_sym(laplMat,1,which="LM",opts = list(retvec = FALSE))
  }
  
  #1/sigma in original paper
  stretchApprox <- 0.5*(1/fewEigenvecs$values[1] + 1/largestEV$values[1])  
  #stretchApprox <- 2*fewEigenvecs$values[kApprox]
  
  strectchApprox.graphLapl.ginv <- fewEigenvecs$vectors[,kApprox:2] %>%
    multiply_by_matrix(diag(1/fewEigenvecs$values[kApprox:2]- stretchApprox, nrow = length(kApprox:2)))  %>%                   tcrossprod(y=fewEigenvecs$vectors[,kApprox:2]) %>% 
    add(diag(nrow(laplMat))*stretchApprox) %>%
    subtract(  tcrossprod(fewEigenvecs$vectors[,kApprox+1])*stretchApprox)
  
  
  if(! is.null(row.names(laplMat))){row.names(strectchApprox.graphLapl.ginv) <- row.names(laplMat)  }
  if(! is.null(colnames(laplMat))){colnames(strectchApprox.graphLapl.ginv) <- colnames(laplMat)  }
  
  res <- list(ginv = strectchApprox.graphLapl.ginv,
              twoNormError = (1/fewEigenvecs$values[1] - 1/largestEV$values[1])*fewEigenvecs$values[kApprox])
  
  class(res) <- "stretchGinvApproximation"
  
  return(res)
}

#' Calculates the 2-norm or spectral norm of a matrix
#'
#' The spectral norm of a matrix is the largest singular value
#'
#' @param mat Input matrix
#' 
#' @importFrom RSpectra svds
#' 
#' @return norm The spectral norm of the matrix
spNorm <- function(mat){
  
  svds(mat,k=1)$d[1]
}
