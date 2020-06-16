#' Calculate the Moore-Penrose psuedo-inverse of a matrix using the divide & conquer algorithm in LAPACK
#'
#' The ubiquitous LAPACK library provides several implementations for the singular-value decomposition (SVD), which underlies
#' psuedo-inverse calculation. The divide and conquer algorithm is particularly quick.
#'
#' @param matrixIn Matrix for which to calculate the psuedo-inverse
#'
#' @return ginv The matrix representation of the psuedo-inverse
#' @export
#'
#' @references \url{http://gallery.rcpp.org/articles/divide-and-concquer-svd/}
#' @references \url{https://en.wikipedia.org/wiki/Moore-Penrose_pseudoinverse}
setGeneric("dcPinv",function(matrixIn) {
  standardGeneric("dcPinv")
})

#' @describeIn dcPinv Psuedo-inverse of base matrix class objects
setMethod("dcPinv",
          c(matrixIn = "matrix"),
          function(matrixIn){
              
            ginv <- dcPinvDense(matrixIn)
            if(! is.null(row.names(matrixIn))){row.names(ginv) <- row.names(matrixIn)  }
            if(! is.null(colnames(matrixIn))){colnames(ginv) <- colnames(matrixIn)  }
            
             return(ginv)
          }
)

#' @describeIn dcPinv Psuedo-inverse of Matrix class objects
setMethod("dcPinv",
          c(matrixIn = "Matrix"),
          function(matrixIn){
            
            ginv <- dcPinvDense(as.matrix(matrixIn))
            if(! is.null(row.names(matrixIn))){row.names(ginv) <- row.names(matrixIn)  }
            if(! is.null(colnames(matrixIn))){colnames(ginv) <- colnames(matrixIn)  }
            
            return(ginv)
          }
)