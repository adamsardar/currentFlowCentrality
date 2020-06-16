context("Testing resistive distance methods")
library(igraph)


set.seed(4+8+15+16+23+24)
testGraph <- simplify(as.undirected(barabasi.game(100)))
testLapl <- graph.laplacian(testGraph,sparse=FALSE)
testLaplPinv <- MASS::ginv(testLapl)

for(i in 1:nrow(testLaplPinv)){
  for(j in 1:ncol(testLaplPinv)){
    
    if( abs(testLaplPinv[i,j]-testLaplPinv[j,i]) > 1E-14) message("i: ",i , "j: ", j)
     
  }
}

test_that("Testing for correctness of solution using matrix properties",{
  
  resDists <- calculateResistiveDistances(testLaplPinv)
  
  for(i in 1:nrow(resDists)){
    
    #Nodes should be distance 0 from themselves
    expect_equal(resDists[i,i],0)
    
    #Don't test ALL nodes - just sample 10 at uniform random
    for(j in sample(1:nrow(resDists),10)){
      
      #For an undirected network, the matrix should be symmetric
      expect_equal(resDists[i,j],resDists[j,i])
      
      for(k in sample(1:nrow(resDists),10)){
        
        diff <- resDists[i,j] - resDists[i,k] + resDists[k,j]
        
        #The triangle inequality should hold for a distnace
        #Allow for tolerance
        expect_true(diff >= -1E-10)
      }
      
    }
  }
  
})


test_that("Test error handling when illegal arguments (such as asymetric laplacians) are provided",{

  testGraphDir <- barabasi.game(30)
  testLaplDir <- graph.laplacian(testGraphDir,sparse=FALSE)
  
  expect_error(calculateResistiveDistances(testLaplDir),regexp = "symmetric")
})


