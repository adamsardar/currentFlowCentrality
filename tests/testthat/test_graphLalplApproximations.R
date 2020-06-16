context("Ensuring that our graph laplacian approximation functions are reliable")
library(igraph)
library(data.table)

test_that("Examining Psuedo-inverse Approximations On The karate graph",{
  
  data(karate,package = "igraphdata")
  karate.igraph <- karate  
  
  graphLapl <- karate.igraph %>% graph.laplacian(sparse = TRUE)
  
  exact.graphLapl.ginv <- graphLapl %>% dcPinv
  
  expect_error(kthCutoffGinvApproximation("Not A Matrix"), regexp = "matrix")
  expect_error(kthCutoffGinvApproximation(graphLapl,k = 10^6), regexp = "less than")
  
  approxKCutoff <- kthCutoffGinvApproximation(graphLapl)
  expect_true(is.matrix(approxKCutoff$ginv))
  expect_equal(row.names(approxKCutoff$ginv),V(karate.igraph)$name)
  expect_equal(colnames(approxKCutoff$ginv),V(karate.igraph)$name)
  expect_true(is.double(approxKCutoff$twoNormError))
  expect_equal(class(approxKCutoff),"kthCutoffGinvApproximation")

  expect_error(stretchGinvApproximation("Not A Matrix"), regexp = "matrix")
  expect_error(stretchGinvApproximation(graphLapl,k = 10^6), regexp = "less than")
  
  expect_error(stretchGinvApproximation("Not A Matrix"))  
  approxStretch <- stretchGinvApproximation(graphLapl)
  expect_true(is.matrix(approxStretch$ginv))
  expect_equal(row.names(approxStretch$ginv),V(karate.igraph)$name)
  expect_equal(colnames(approxStretch$ginv),V(karate.igraph)$name)
  expect_true(is.double(approxStretch$twoNormError))
  expect_equal(class(approxStretch),"stretchGinvApproximation")
  
  approxPerf <- data.table(k = 2:8)
  
  iCount <- 1L
  for(kApp in 2:8){
    set(approxPerf,iCount,"kthAprrox",list(list(kthCutoffGinvApproximation(graphLapl,kApp))))
    iCount <- iCount + 1L
  }
  
  approxPerf[,predictedKError := kthAprrox[[1]]$twoNormError,by=k]
  approxPerf[,actualKError := spNorm(exact.graphLapl.ginv - kthAprrox[[1]]$ginv)/spNorm(exact.graphLapl.ginv),by=k]
  
  approxPerf[,stretchAprrox := list(list(stretchGinvApproximation(graphLapl,k,exact=TRUE))),by = k]
  approxPerf[,predictedStretchError := stretchAprrox[[1]]$twoNormError,by=k]
  approxPerf[,actualStretchError := spNorm(exact.graphLapl.ginv - stretchAprrox[[1]]$ginv)/spNorm(exact.graphLapl.ginv),by=k]
  
  expect_true(approxPerf[,all(round(predictedKError,digits=3) == round(actualKError,digits=3))], info = "Spectral norm error should match perfectly with prediction")
  expect_true(approxPerf[,all(predictedStretchError >= actualStretchError)], info = "Stretch approx. error is simply and upper bound")

  expect_true(approxPerf[,all(predictedKError >= predictedStretchError)], info = "Stretch approximation is supposed to be better than the k'th cutoff")  
  expect_true(approxPerf[,all(actualKError >= actualStretchError)], info = "Stretch approximation is supposed to be better than the k'th cutoff")
})


test_that("Look at the dolphin social network as Bozzo et al do in their paper",{

  data(dolphinPodLusseau)

  exactDolphinLapl <- dolphinPodLusseau %>% graph.laplacian(sparse=FALSE) %>% dcPinv
  expect_output(str(exactDolphinLapl),regexp="num")
  expect_false(is.null(row.names(exactDolphinLapl)))
  expect_false(is.null(colnames(exactDolphinLapl)))
  
  exactCurrentFlowCentralities <- calculateAverageCurrentFlowBetweeness(exactDolphinLapl,get.adjacency(dolphinPodLusseau,sparse=FALSE)) %>% as.vector
  names(exactCurrentFlowCentralities) <- V(dolphinPodLusseau)$name
  
  exactTopTenDolphins <- sort(exactCurrentFlowCentralities,decreasing = TRUE) %>% head(n=10) %>% names
  expect_equal(exactTopTenDolphins[1:6],c("Beescratch","SN100","Jet","SN9","Web","DN63"))
  
  ####################
  
  approxDolphinLapl <- dolphinPodLusseau %>% graph.laplacian(sparse=TRUE) %>% stretchGinvApproximation(kApprox = 3)
  expect_output(str(approxDolphinLapl),regexp="num")
  expect_false(is.null(row.names(approxDolphinLapl$ginv)))
  expect_false(is.null(colnames(approxDolphinLapl$ginv)))
  
  approxCurrentFlowCentralities <- calculateAverageCurrentFlowBetweeness(approxDolphinLapl$ginv,get.adjacency(dolphinPodLusseau,sparse=FALSE)) %>% as.vector
  names(approxCurrentFlowCentralities) <- V(dolphinPodLusseau)$name
  
  
  
  approxTopTenDolphins <- sort(approxCurrentFlowCentralities,decreasing = TRUE) %>% head(n=10) %>% names
  
  expect_true(length(which(approxTopTenDolphins %in% exactTopTenDolphins))/length(approxTopTenDolphins) >= 0.7)
  ####################
  
  approxDolphinLapl <- dolphinPodLusseau %>% graph.laplacian(sparse=TRUE) %>% stretchGinvApproximation(kApprox = 6)
  expect_output(str(approxDolphinLapl),regexp="num")
  expect_false(is.null(row.names(approxDolphinLapl$ginv)))
  expect_false(is.null(colnames(approxDolphinLapl$ginv)))
  
  approxCurrentFlowCentralities <- calculateAverageCurrentFlowBetweeness(approxDolphinLapl$ginv,get.adjacency(dolphinPodLusseau,sparse=FALSE)) %>% as.vector
  names(approxCurrentFlowCentralities) <- V(dolphinPodLusseau)$name
  
  approxTopTenDolphins <- sort(approxCurrentFlowCentralities,decreasing = TRUE) %>% head(n=10) %>% names
  expect_true(length(which(approxTopTenDolphins %in% exactTopTenDolphins))/length(approxTopTenDolphins) >= 0.7)
  ####################
  
  approxDolphinLapl <- dolphinPodLusseau %>% graph.laplacian(sparse=TRUE) %>% stretchGinvApproximation(kApprox = 12)
  approxCurrentFlowCentralities <- calculateAverageCurrentFlowBetweeness(approxDolphinLapl$ginv,get.adjacency(dolphinPodLusseau,sparse=FALSE)) %>% as.vector
  names(approxCurrentFlowCentralities) <- V(dolphinPodLusseau)$name
  
  approxTopTenDolphins <- sort(approxCurrentFlowCentralities,decreasing = TRUE) %>% head(n=10) %>% names
  
  expect_true(length(which(approxTopTenDolphins %in% exactTopTenDolphins))/length(approxTopTenDolphins) >= 0.7)
  ####################
  
  approxDolphinLapl <- dolphinPodLusseau %>% graph.laplacian(sparse=TRUE) %>% kthCutoffGinvApproximation(kApprox = 12)
  expect_output(str(approxDolphinLapl),regexp="num")
  expect_false(is.null(row.names(approxDolphinLapl$ginv)))
  expect_false(is.null(colnames(approxDolphinLapl$ginv)))
  
  approxCurrentFlowCentralities <- calculateAverageCurrentFlowBetweeness(approxDolphinLapl$ginv,get.adjacency(dolphinPodLusseau,sparse=FALSE)) %>% as.vector
  names(approxCurrentFlowCentralities) <- V(dolphinPodLusseau)$name
  
  approxTopTenDolphins <- sort(approxCurrentFlowCentralities,decreasing = TRUE) %>% head(n=10) %>% names
  
  expect_true(length(which(approxTopTenDolphins %in% exactTopTenDolphins))/length(approxTopTenDolphins) >= 0.7)
})