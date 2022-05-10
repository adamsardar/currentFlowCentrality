<!-- badges: start -->
[![R-CMD-check](https://github.com/adamsardar/currentFlowCentrality/workflows/R-CMD-check/badge.svg)](https://github.com/adamsardar/currentFlowCentrality/actions)
[![Codecov test coverage](https://codecov.io/gh/adamsardar/currentFlowCentrality/branch/master/graph/badge.svg)](https://codecov.io/gh/adamsardar/currentFlowCentrality?branch=master)
[![License: MPL 2.0](https://img.shields.io/badge/License-MPL%202.0-brightgreen.svg)](https://opensource.org/licenses/MPL-2.0)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-blue.svg)](https://www.tidyverse.org/lifecycle/#experimental)
  <!-- badges: end -->
  
A little package for calculating the current-flow centrality of nodes in a network. This extends the concept
of betweeness and closeness based in shortest paths to those based on all paths, with the intuitive interpretation
of electrical current in an electrical network.

### Installation

```
> remotes::install_github("adamsardar/currentFlowCentrality")
```

### Reference (for the algorithm):

Newman, M. E. J. (2005). A measure of betweenness centrality based on random walks. Social Networks, 27(1), 39-54. doi:10.1016/j.socnet.2004.11.009


