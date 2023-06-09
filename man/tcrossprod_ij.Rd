% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sparse.R
\name{tcrossprod_ij}
\alias{tcrossprod_ij}
\title{Compute cross product between two sets of rows of a matrix.}
\usage{
tcrossprod_ij(X, id1, id2)
}
\arguments{
\item{X}{matrix}

\item{id1}{vector of integers specifying the list of rows of \code{X}
(first set)}

\item{id2}{vector of integers specifying the list of rows of \code{X},
(second set), same length as \code{id1}.}
}
\value{
matrix containing the cross product of \code{X[id1, ]} and
\code{X[id2, ]}.
}
\description{
\code{tcrossprod_ij} computes cross product between two sets of rows of a
matrix.
}
\examples{

set.seed(42)
X <- matrix(rnorm(5 * 3), 5, 3)

id1 <- c(1, 3)
id2 <- c(5, 4)

(s1 <- matric::tcrossprod_ij(X, id1, id2))

(s2 <- tcrossprod(X)[id1, id2])

all.equal(s1, s2)
}
