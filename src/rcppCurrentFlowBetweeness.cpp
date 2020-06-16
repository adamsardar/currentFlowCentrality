#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;


//' @title generate the Moore-Penrose psueoinverse using the armadillo package
//' @description The Armadillo package provides a low level interface to native LAPACK/BLAS routines
//' and we plan to use the divide-and-conquer implementation of SVD to generate the psuedo-inverse
//' @param X A matrix to be psuedo-inverted
//' @return An armadillo (dense) matrix containg the psuedo-inverse values 
//' @references \url{http://en.wikipedia.org/wiki/Moore-Penrose_pseudoinverse}
//' @references \url{http://gallery.rcpp.org/articles/divide-and-concquer-svd/}
// [[Rcpp::export]]
arma::mat dcPinvDense(arma::mat X) {
  
  double  tolerance = 1e-11;
  
  arma::mat U, V;
  
  arma::vec S,invS;
  arma::svd(U, S, V, X, "dc");
  
  invS.zeros(X.n_rows); //Preallocate
  
  arma::uvec indiciesAboveTol = find(S >= tolerance);
  
  invS.elem(indiciesAboveTol) = 1/S.elem(indiciesAboveTol);
  
  return V * diagmat(invS) * U.t();
}

//' @title generate the Moore-Penrose psueoinverse using the armadillo package
//' @description The Armadillo package provides a low level interface to native LAPACK/BLAS routines
//' It provides several implementations of SVD and this function uses the standard algorithm
//' @param X A matrix to be psuedo-inverted
//' @return An armadillo (dense) matrix containg the psuedo-inverse values 
//' @references \url{http://en.wikipedia.org/wiki/Moore-Penrose_pseudoinverse}
//' @references \url{http://gallery.rcpp.org/articles/divide-and-concquer-svd/}
// [[Rcpp::export]]
arma::mat standardPinv(arma::mat X) {
  
  double  tolerance = 1e-8;
  
  arma::mat U, V;
  arma::vec S, invS;
  arma::svd(U, S, V, X, "standard");
  
  invS.zeros(X.n_rows); //Preallocate
  
  arma::uvec indiciesAboveTol = find(S >= tolerance);
  
  invS.elem(indiciesAboveTol) = 1/S.elem(indiciesAboveTol);
  
  return V * diagmat(invS) * U.t();
}

//' @title generate the average resistive distance per node
//' @description We take advantage of the notion of resistor network to define current-flow distance.
//' Given nodes \eqn{i} and \eqn{j}, the resistance distance \eqn{R_{i,j}}{R_ij} between \eqn{i} and \eqn{j}
//' is the effective (or equivalent) resistance between \eqn{i} and \eqn{j} when a unit of current is injected
//' from source \eqn{s} and removed from target \eqn{t}.
//' @param graphLaplacianInv The psuedoinverse of the graph laplacian. The graph MUST be fully connected - I make no checks for you ...
//' @return currentFlowDistances A column matrix of resistive distances from the network
//' @references \url{http://www.sci.unich.it/~francesc/teaching/network/flowcentrality.html}
//' @export
// [[Rcpp::export]]
NumericVector calculateCurrentFlowDistance(NumericMatrix graphLaplacianInv){
  
  arma::mat  Gpinv = as<arma::mat>(wrap(graphLaplacianInv));
  
  if(!Gpinv.is_square()){throw std::invalid_argument("Laplacian matricies (and hence their psuedo-inverses) are square.");}
  
  arma::uvec rowSumsAboveEta = find(abs(sum(Gpinv,0)) >= 1E-11);
  arma::uvec colSumsAboveEta = find(abs(sum(Gpinv,1)) >= 1E-11);
  
  if( rowSumsAboveEta.n_elem > 0){throw std::invalid_argument("Rowsums of a graph laplacian (and it's psuedo-inverse - they share the same null space) must be zero!");}
  if( colSumsAboveEta.n_elem > 0){throw std::invalid_argument("Colsums of a graph laplacian (and it's psuedo-inverse - they share the same null space) must be zero!");}
  
  int numberOfNodes = Gpinv.n_rows;
  double commonTerm = trace(Gpinv)/((double) numberOfNodes);
  
  arma::vec currentFlowDistances = arma::vec(numberOfNodes);
  currentFlowDistances.fill(commonTerm);
  
  currentFlowDistances += Gpinv.diag();
  
  NumericVector returnCurrentFlowDistances = NumericVector(currentFlowDistances.begin(),currentFlowDistances.end());
  returnCurrentFlowDistances.names() = rownames(graphLaplacianInv);
  
  return returnCurrentFlowDistances;
}



//' @title for a single source and single sink node of unit current, calculate the per-node flow
//' @description We take advantage of the notion of resistor network to define current-flow. It is defined as
//' \deqn{F_{i}^{(st)}=\frac{1}{2} \sum_{j} A_{ij} \left | V_{i}^{(st)} - V_{j}^{(st)}  \right |}
//' Given nodes \eqn{i} and \eqn{j} and a unit of current injected from source \eqn{s} and removed from target \eqn{t}.
//' @param Gpinv The psuedoinverse of the (undirected) graph laplacian. The graph MUST be fully connected - I make no checks for you ...
//' @param adjacency The adjacency matrix of a (undirected) graph
//' @param sourceIndex The integer index of the source of unit flow
//' @param sinkIndex The integer index of the sink of unit flow
//' @return currentFlowDistances A column matrix of resistive distances from the network
//' @references \url{http://www.sci.unich.it/~francesc/teaching/network/flowcentrality.html}
// [[Rcpp::export]]
arma::vec singleSourceTargetFlow(arma::mat Gpinv, arma::sp_mat adjacency, int sourceIndex, int sinkIndex){

  if(sourceIndex+1 > Gpinv.n_cols){throw std::invalid_argument("sourceIndex must be less than or equal to the number of rows");}
  
  if(sinkIndex+1 > Gpinv.n_cols){throw std::invalid_argument("sinkIndex must be less than or equal to the number of rows/columns");}
  
  int nNodes = Gpinv.n_rows;
  arma::vec flowForSingleSourceSink = arma::vec(nNodes);

  for(int i; i < nNodes; ++i){

    arma::rowvec sourceRow = Gpinv.row(sourceIndex);  
    arma::rowvec sinkRow = Gpinv.row(sinkIndex);
  
    double potentialAti =  sourceRow[i] - sinkRow[i];
    
    arma::sp_mat::const_row_iterator row_it     = adjacency.begin_row(i);
    arma::sp_mat::const_row_iterator row_it_end = adjacency.end_row(i);
    
    //Non-zero columns of the adjacency matrix correspond to nodes that are adjacent to the node i
    double singleNodeFlow = 0;
    for(; row_it != row_it_end; ++row_it){
      
      if(*row_it > 0){
        
        double potentialsAtAdjacentNode = sourceRow(row_it.col()) - sinkRow(row_it.col());
        
        singleNodeFlow += std::abs(potentialAti - potentialsAtAdjacentNode);
      }
    }
    
    flowForSingleSourceSink[i] = 0.5*singleNodeFlow;
  }

  flowForSingleSourceSink[sourceIndex] *= 2;
  flowForSingleSourceSink[sinkIndex] *= 2;
  
  return flowForSingleSourceSink;  
}

//' @title calculate the average current flow betweeness
//' @description We take advantage of the notion of resistor network to define current-flow. It is defined as
//' \deqn{F_{i}^{(st)}=\frac{1}{2} \sum_{j} A_{ij} \left | V_{i}^{(st)} - V_{j}^{(st)}  \right |}
//' Given nodes \eqn{i} and \eqn{j} and a unit of current injected from source \eqn{s} and removed from target \eqn{t}.
//' The current flow betweeness is defined as: \deqn{b_{i} = \frac{\sum_{s<t} F_{i}^{(st)} }{\frac{1}{2}n(n-1))}}. Hence
//' averageing over all \eqn{s} - \eqn{t} pairs.
//' @inheritParams singleSourceTargetFlow
//' @export
//' @param graphLaplacianInv The psuedoinverse of the graph laplacian. The graph MUST be fully connected - I make no checks for you ...
//' @return currentFlowBetweeness A vector of current flow betweeness, averaged over all source-sink pairs
//' @references \url{http://www.sci.unich.it/~francesc/teaching/network/flowcentrality.html}
// [[Rcpp::export]]
NumericVector calculateAverageCurrentFlowBetweeness(NumericMatrix graphLaplacianInv, arma::sp_mat adjacency){
  
  arma::mat Gpinv = as<arma::mat>(wrap(graphLaplacianInv));
  
  if(!Gpinv.is_square()){throw std::invalid_argument("Laplacian matricies (and hence their psuedo-inverses) are square.");}
  if(!adjacency.is_square()){throw std::invalid_argument("Adjacency matricies are by definition square.");}
  
  arma::uvec rowSumsAboveEta = find(abs(sum(Gpinv,0)) >= 1E-11);
  arma::uvec colSumsAboveEta = find(abs(sum(Gpinv,1)) >= 1E-11);
  
  if( rowSumsAboveEta.n_elem > 0){throw std::invalid_argument("Rowsums of a graph laplacian (and it's psuedo-inverse - they share the same null space) must be zero!");}
  if( colSumsAboveEta.n_elem > 0){throw std::invalid_argument("Colsums of a graph laplacian (and it's psuedo-inverse - they share the same null space) must be zero!");}
  
  if( Gpinv.n_rows != adjacency.n_rows || Gpinv.n_cols != adjacency.n_cols ){throw std::invalid_argument("Gpinv and adjacency matricies must be the same size!");}
  
  int nNodes = Gpinv.n_rows;
  arma::vec currentFlowBetweeness = arma::vec(nNodes).zeros();
  
  //Iterate over all unique source-sink pairs and keep a cumulative sum
  for(int source = 0; source < nNodes; ++source){
    for(int sink = source+1; sink < nNodes; ++sink){
      
      currentFlowBetweeness += singleSourceTargetFlow(Gpinv,adjacency,source,sink);
    }
  }
  
  //divide by a normalising factor
  currentFlowBetweeness /= 0.5*((double) nNodes*(nNodes - 1));
  
  NumericVector returnCurrentFlowBetweeness = NumericVector(currentFlowBetweeness.begin(),currentFlowBetweeness.end());
  returnCurrentFlowBetweeness.names() = rownames(graphLaplacianInv) ;
 
  return returnCurrentFlowBetweeness;
}


//' @title generate the pairwise resistive distance
//' @description We take advantage of the notion of resistor network to define current-flow distance.
//' Given nodes \eqn{i} and \eqn{j}, the resistance distance \eqn{R_{i,j}}{R_ij} between \eqn{i} and \eqn{j}
//' is the effective (or equivalent) resistance between \eqn{i} and \eqn{j} when a unit of current is injected
//' from source \eqn{s} and removed from target \eqn{t}.
//' @inheritParams calculateCurrentFlowDistance
//' @return currentFlowDistances A matrix of pairwise resistive distances from the network
//' @references \url{http://www.sci.unich.it/~francesc/teaching/network/flowcentrality.html}
//' @export
// [[Rcpp::export]]
NumericMatrix calculateResistiveDistances(NumericMatrix graphLaplacianInv){
  
  arma::mat  Gpinv = as<arma::mat>(wrap(graphLaplacianInv));
  
  int nNodes = Gpinv.n_rows;
  
  if(!Gpinv.is_square()){throw std::invalid_argument("Laplacian matricies (and hence their psuedo-inverses) are square.");}
  
  if(!Gpinv.is_symmetric() ){throw std::invalid_argument("At current this method is only defined for undirected networks (hence symmetric laplacian psuedoinverse matricies)");}
  
  arma::mat resistiveDistances;
  
  resistiveDistances.zeros(nNodes,nNodes);
  
  arma::vec GpinvDiag = Gpinv.diag();
  
  for(int i =0; i < nNodes; ++i){
    for(int j =0; j < nNodes; ++j){
      
      resistiveDistances(i,j) = GpinvDiag[i] + GpinvDiag[j] - 2*Gpinv(i,j); 
    }
  }
  
  NumericMatrix returnResistiveDistances = as<NumericMatrix>(wrap(resistiveDistances));
  rownames(returnResistiveDistances) = rownames(graphLaplacianInv);
  colnames(returnResistiveDistances) = colnames(graphLaplacianInv);
  
  return returnResistiveDistances;
}
