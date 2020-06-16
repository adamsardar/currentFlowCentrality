context("The psuedo-inverse matrix (aka the Moore-Penrose inverse) is at the heart of these current-flow centralities")
library(igraph)
library(Matrix)

set.seed(4+8+15+16+23+24)

test_that("Testing for correstness of solution",{
  
  M <- matrix(rnorm(50*50), 50, 50)
  expect_output(str(dcPinv(M)),regexp="num")
  
  expect_equal(MASS::ginv(M),dcPinv(M))
  expect_equal(MASS::ginv(M),standardPinv(M))
})

test_that("Checking sparse dcPinv",{
  
  data(dolphinPodLusseau)
  
  sampleSparseMatrix <- dolphinPodLusseau %>% graph.laplacian(sparse=TRUE)
  laplGinvSp <- dcPinv(sampleSparseMatrix)
  expect_output(str(laplGinvSp),regexp="num")
  expect_false(is.null(row.names(laplGinvSp)))
  expect_false(is.null(colnames(laplGinvSp)))
  
  sampleDenseMatrix <- dolphinPodLusseau %>% graph.laplacian(sparse=FALSE)
  laplGinvDen <- dcPinv(sampleDenseMatrix)
  expect_output(str(laplGinvDen),regexp="num")  
  expect_false(is.null(row.names(laplGinvDen)))
  expect_false(is.null(colnames(laplGinvDen)))
  
  standardMat <- matrix(sampleDenseMatrix,nrow=nrow(sampleDenseMatrix))
  laplGinvStan <- dcPinv(sampleDenseMatrix)
  expect_output(str(laplGinvDen),regexp="num")
  expect_false(is.null(row.names(laplGinvStan)))
  expect_false(is.null(colnames(laplGinvStan)))
  
  laplStan <- dcPinv(laplGinvDen)
  
  expect_lt(norm(laplStan - sampleDenseMatrix),0.01)
  
  expect_equal(laplGinvSp,laplGinvDen)
  expect_equal(laplGinvStan,laplGinvDen)
})