context("Looking at a single source/sink scenario")
library(igraph)
library(magrittr)

data("singleSTcurrentFlow")

set.seed(4+8+15+16+23+42)

test_that("Compare calculated results to a reference",{
  
  adjMat <- get.adjacency(singleSTcurrentFlow)
  adjMat %<>% .[V(singleSTcurrentFlow)$name,V(singleSTcurrentFlow)$name]
  Glapl <- graph.laplacian(singleSTcurrentFlow)
  Glapl %<>% .[V(singleSTcurrentFlow)$name,V(singleSTcurrentFlow)$name]
  
  Gpinv <- MASS::ginv( as.matrix(Glapl) )
  colnames(Gpinv) <- colnames(Glapl) 
  row.names(Gpinv) <- row.names(Glapl) 
  
  sourceID <- which(V(singleSTcurrentFlow)$isSource) 
  sinkID <- which(V(singleSTcurrentFlow)$isSink)
  
  flowVals <- singleSourceTargetFlow( Gpinv, Matrix(adjMat),sourceID-1,sinkID-1) %>% as.vector()
  names(flowVals) <- colnames(adjMat)
  
  referenceValue = V(singleSTcurrentFlow)$absoluteFlow
  names(referenceValue) <- V(singleSTcurrentFlow)$name
  #covr seems to affect compiled code ALOT
  if(!"covr" %in% loadedNamespaces()){
  expect_equal(round(flowVals,digits=2),referenceValue)
  }
    
  directedGraph <- barabasi.game(50)
  directedGraphAdj <- directedGraph %>% get.adjacency(sparse=FALSE)
  directedGraphLapl <- directedGraph %>% graph.laplacian(sparse=FALSE)
  directedGraphLaplPinv <- directedGraphLapl %>% MASS::ginv()
  
  expect_error(singleSourceTargetFlow(directedGraphLaplPinv,directedGraphAdj,20,40),regexp = "directed",
               info = "Function should spit and error if we give a directed graph to it")
  
  undirectedGraph <- barabasi.game(50) %>% as.undirected()
  undirectedGraphAdj <- undirectedGraph %>% get.adjacency(sparse=FALSE)
  undirectedGraphLapl <- undirectedGraph %>% graph.laplacian(sparse=FALSE)
  undirectedGraphLaplPinv <- undirectedGraphLapl %>% MASS::ginv()
  
  expect_output(singleSourceTargetFlow(undirectedGraphLaplPinv,undirectedGraphAdj,20,40) %>% str,"num",
                info = "Just making sure that the function works with an undirected graph")
  
  expect_equal(singleSourceTargetFlow(undirectedGraphLaplPinv,undirectedGraphAdj,20,20) %>% as.vector,rep(0,50),
               info = "For the same source and sink node, there should be no current flowing anywhere")
  
  expect_error(singleSourceTargetFlow(undirectedGraphLaplPinv,undirectedGraphAdj,20000,40),regexp = "less",
               info = "Check that we create an error if we ask for a source value that's far out")
  
  expect_error(singleSourceTargetFlow(undirectedGraphLaplPinv,undirectedGraphAdj,40,20000),regexp = "less",
               info = "Check that we create an error if we ask for a sink value that's far out")
})


test_that("Check what happens on a disconnected network",{

  set.seed(4+8+15+16+23+42)
  
  disjointGraph <- graph.disjoint.union(barabasi.game(20),barabasi.game(20)) %>% as.undirected
  
  V(disjointGraph)$name <- 1:vcount(disjointGraph)
  
  undirectedGraphAdj <- disjointGraph %>% get.adjacency(sparse=FALSE)
  undirectedGraphLapl <- disjointGraph %>% graph.laplacian(sparse=FALSE)
  
  undirectedGraphLaplPinv <- undirectedGraphLapl %>% MASS::ginv()
  
  sourceNode = 19
  sinkNode = 18
  disjointGraphSingeFlow <- singleSourceTargetFlow(undirectedGraphLaplPinv,undirectedGraphAdj,sourceNode-1,sinkNode-1) %>% as.vector
  names(disjointGraphSingeFlow) <- V(disjointGraph)$name
  
  seperateNodeIds <- 21:40 %>% as.character
  emptyFlowRef <- rep(0,20)
  names(emptyFlowRef) <- seperateNodeIds
  
  expect_equal(disjointGraphSingeFlow[seperateNodeIds],emptyFlowRef, info = "For the given pair of node ids, no current should flow through the disjoint part of the network")

  nodesOnPath <- rep(1,5)
  names(nodesOnPath) <- c("19","1","2","3","18")
  #covr seems to affect compiled code ALOT
  if(!"covr" %in% loadedNamespaces()){
   expect_equal(disjointGraphSingeFlow[c("19","1","2","3","18")],nodesOnPath,info="For a simple small world network with few branches, ther current flow along a path should be 1")
  }
})
  

test_that("Some simple tests that make use of properties of disconnected graphs",{
  
  set.seed(4+8+15+16+23+42)
  
  orphanNodeGraph <- graph.disjoint.union(barabasi.game(20),barabasi.game(1)) %>% as.undirected
  V(orphanNodeGraph)$name <- 1:21 %>% as.character
  
  orphanNodeGraphAdj <- orphanNodeGraph %>% get.adjacency(sparse=FALSE)
  orphanNodeGraphLapl <- orphanNodeGraph %>% graph.laplacian(sparse=FALSE)
  
  orphanNodeGraphLaplPinv <- orphanNodeGraphLapl %>% MASS::ginv()
  
  sourceNode = 17
  sinkNode = 2
  orphanNodeGraphSingeFlow <- singleSourceTargetFlow(orphanNodeGraphLaplPinv,orphanNodeGraphAdj,sourceNode-1,sinkNode-1) %>% as.vector
  names(orphanNodeGraphSingeFlow) <- V(orphanNodeGraph)$name
  
  nodesOnPath <- rep(1,4)
  names(nodesOnPath) <- c("17","8","1","2")
  if(!"covr" %in% loadedNamespaces()){
  expect_equal(orphanNodeGraphSingeFlow[c("17","8","1","2")],nodesOnPath,info="For a simple small world network with few branches, ther current flow along a path should be 1")
  }
    
  sourceNode = 21
  sinkNode = 18
  orphanNodeGraphSingeFlow <- singleSourceTargetFlow(orphanNodeGraphLaplPinv,orphanNodeGraphAdj,sourceNode-1,sinkNode-1) %>% as.vector
  names(orphanNodeGraphSingeFlow) <- V(orphanNodeGraph)$name
})