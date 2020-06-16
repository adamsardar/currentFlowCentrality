context("The resistive distance is a popular way of representing average disnace of a ndoe from a network but using all paths")

data(karate, package = "igraphdata")
karateGraphPinv <- dcPinv(igraph::graph.laplacian(karate))

test_that("Test error handling of errors",{

  expect_error(calculateCurrentFlowDistance(matrix(0,2,3)),info = "Graph laplacians (and their inverses) must be square.")

  expect_error(calculateCurrentFlowDistance(matrix(rnorm(5*5),5,5)),info= "A random matrix is not a graph laplacian and an error should be thrown")
  
  rowSumsZeroMat <- matrix(0,11,11)
  colSumsZeroMat <- matrix(0,11,11)
  for(i in 1:11){
    rowSumsZeroMat[i,] <- sample(seq(-5,5),11)
    colSumsZeroMat[,i] <- sample(seq(-5,5),11)
  }
  
  expect_error(calculateCurrentFlowDistance(rowSumsZeroMat),info = "A matrix without zero colsums should throw an error")
  expect_error(calculateCurrentFlowDistance(colSumsZeroMat),info = "A matrix without zero rowsums should throw an error")
  
  currentFlowDistances <- calculateCurrentFlowDistance(karateGraphPinv)
  expect_equal(names(currentFlowDistances),V(karate)$name)
  expect_output(str(currentFlowDistances),"num")
  
  karate %<>% delete_vertex_attr("name")
  karateGraphPinv <- dcPinv(igraph::graph.laplacian(karate))
  currentFlowDistances <- calculateCurrentFlowDistance(karateGraphPinv)
  expect_true(is.null(names(currentFlowDistances)))
  expect_output(str(currentFlowDistances),"num")
})