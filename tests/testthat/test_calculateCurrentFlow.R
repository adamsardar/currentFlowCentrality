context("Looking at a single source/sink scenario")
library(igraph)

test_that("Testing the structure and error handling of the main user interface to the current flow methods",{
  
  set.seed(4+8+15+16+23+42)
  disjointGraph <- graph.disjoint.union(barabasi.game(20),barabasi.game(20)) %>% as.undirected
  V(disjointGraph)$name <- as.character( 1:vcount(disjointGraph))

  expect_error(calculateCurrentFlow(disjointGraph),rexp = "connected",info = "Only connected graphs are accepted at this time.")
  
  connectedGraph <- graph.disjoint.union(barabasi.game(40)) %>% as.undirected
  V(connectedGraph)$name <-as.character( 1:vcount(connectedGraph))
  
  connectedGraphFlowBetweeness <-  calculateCurrentFlow(connectedGraph)
  expect_output(str(connectedGraphFlowBetweeness),regexp = "List of 2")
  expect_output(str(connectedGraphFlowBetweeness),regexp = "distance")
  expect_output(str(connectedGraphFlowBetweeness),regexp = "betweeness")
  expect_equal(length(connectedGraphFlowBetweeness$betweeness),length(connectedGraphFlowBetweeness$distance))
  
  expect_equal(names(connectedGraphFlowBetweeness$distance),V(connectedGraph)$name)
  expect_equal(names(connectedGraphFlowBetweeness$betweeness),V(connectedGraph)$name)
})