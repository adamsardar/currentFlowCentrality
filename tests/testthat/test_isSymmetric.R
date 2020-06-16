context("Some quick checks of isSymmetricMatrix")
library(igraph)

test_that("Just make sure that the isSymmetricMatrix method works",{
  
  set.seed(4+8+15+16+23+24)
  directedGraph <- barabasi.game(50)
  directedGraphAdj <- directedGraph %>% get.adjacency(sparse=FALSE)
  
  expect_false(isSymmetricMatrix(directedGraphAdj))
  expect_true(isSymmetricMatrix(directedGraph %>% as.undirected %>% get.adjacency(sparse=FALSE)))
  
  expect_error(isSymmetricMatrix(matrix(rnorm(3*10),3,10)),regexp = "square")
})
