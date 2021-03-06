% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/laplInvApprox.R
\name{kthCutoffGinvApproximation}
\alias{kthCutoffGinvApproximation}
\title{Approximate the psuedo-inverse of a laplacian matrix from a handful of eigenvectors and values}
\usage{
kthCutoffGinvApproximation(laplMat, kApprox = 8, exact = FALSE)
}
\arguments{
\item{laplMat}{A laplacian matrix for which to calculate the approximate psudeo-inverse. Since we only handle undirected graphs, this must be symmetric}

\item{kApprox}{Number of eigenvectors/values to use in construction of the psuedo-inverse}

\item{exact}{Logical flag to turn off use of rARPACK in generation of the desired few eigenvectors/values. This is very useful for testing}
}
\value{
A list of two entries: ginv - a matrix of the approximate graph laplacian approximated using the first k eigenvectors and twoNormError,
the spectral norm error in the matrix estimate
}
\description{
Based on work by Bozzo & Franceschet who showed that one could closely approximate the psuedo-inverse of a graph laplacian using a small (k << n)
number of eigenvalues and vectors. Consult their paper for details
}
\examples{
data(dolphinPodLusseau)

dolphnLapl <- igraph::graph.laplacian(dolphinPodLusseau)

str(kthCutoffGinvApproximation(dolphnLapl))

}
\references{
Bozzo, E., & Franceschet, M. (2012). Effective and efficient approximations of the generalized inverse of the graph Laplacian matrix with an application to current-flow betweenness centrality. Physics and Society
}
\seealso{
Other approximatePsuedoInverse: \code{\link{stretchGinvApproximation}}
}
\concept{approximatePsuedoInverse}
