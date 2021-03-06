% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calculateCurrentFlow.R
\name{calculateCurrentFlow}
\alias{calculateCurrentFlow}
\title{Calculate the current flow distance and betweeness, as popularised by Newman et al.}
\usage{
calculateCurrentFlow(network)
}
\arguments{
\item{network}{The network for wich to calculate current flow centrality for. If directed then it will be made undirected with a warning.}
}
\value{
list A list of two fatures: betweeness and distance
}
\description{
This is a nice little wrapper function that delegates out computation to c++ implementations
}
\references{
Newman, M. E. J. (2005). A measure of betweenness centrality based on random walks. Social Networks, 27(1)
}
