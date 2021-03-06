% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MGRAF1.R
\name{MGRAF1}
\alias{MGRAF1}
\title{First variant of M-GRAF model}
\usage{
MGRAF1(A, K, tol, maxit)
}
\arguments{
\item{A}{Binary array with size VxVxn storing the VxV symmetric adjacency
matrices of n graphs.}

\item{K}{An integer that specifies the latent dimension of the graphs}

\item{tol}{A numeric scalar that specifies the convergence threshold of CISE
algorithm. CISE iteration continues until the absolute percent change in
joint log-likelihood is smaller than this value. Default is tol = 0.01.}

\item{maxit}{An integer that specifies the maximum number of iterations.
Default is maxit = 5.}
}
\value{
A list is returned containing the ingredients below from M-GRAF1
  model corresponding to the largest log-likelihood over iterations.
  \item{Z}{A numeric vector containing the lower triangular entries in the
  estimated matrix Z.} \item{Lambda}{Kxn matrix where each column stores the
  diagonal entries in \eqn{\Lambda_i}.} \item{Q}{VxKxn array containing the
  estimated VxK orthonormal matrix \eqn{Q_i}, i=1,...,n.} \item{D_LT}{Lxn
  matrix where each column stores the lower triangular entries in \eqn{D_i =
  Q_i * \Lambda_i * Q_i^{\top}}; L=V(V-1)/2.} \item{LL_max}{Maximum
  log-likelihood across iterations.} \item{LL}{Joint log-likelihood at each
  iteration.}
}
\description{
\code{MGRAF1} returns the estimated common structure Z and subject-specific
low rank components \eqn{Q_i} and \eqn{\Lambda_i} for multiple undirected
graphs.
}
\details{
The subject-specific deviation \eqn{D_i} is decomposed into \deqn{D_i = Q_i
* \Lambda_i * Q_i^{\top},} where each \eqn{Q_i} is a VxK orthonormal matrix and
each \eqn{\Lambda_i} is a KxK diagonal matrix.
}
\examples{
data(A)
res = MGRAF1(A=A, K=5, tol=0.01, maxit=5)

}
