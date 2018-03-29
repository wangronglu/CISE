# CISE
CISE is an R package that fits the multiple-gragh factorization (M-GRAF) model ( [Wang et al., 2017](https://arxiv.org/abs/1707.06360))
for separating common and individual low rank structure for multiple graphs with a common set of nodes.

## Installation
```
install.packages("CISE")

# Or the development version from GitHub:
# install.packages("devtools")
devtools::install_github("wangronglu/CISE")
```

## Usage
The package contains a dataset of binary brain networks of 212 subjects. After loading the data, we can apply 3 variants of M-GRAF model 
to decompose the data. Model checking is needed to evaluate which model provides the best fit.
```
library(CISE)
# load data
# A is an array of adjacency matrices of multiple graphs
data(A)

? MGRAF2
# specify low rank K=5
res = MGRAF2(A = A, K=5, tol=0.01, maxit=5)

# check the common structure - baseline propensities of edges between pairwise nodes
res$Z
# check low rank component for each graph
res$Q
# check joint log-likelihood across iterations
res$LL
```
