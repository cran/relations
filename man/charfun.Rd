\name{charfun}
\alias{relation_charfun}
\title{Relation Characteristic Functions}
\description{
  Determine the characteristic function of a relation.
}
\usage{
relation_charfun(x, components = FALSE)
}
\arguments{
  \item{x}{an object inheriting from class \code{\link{relation}}.}
  \item{components}{a logical indicating whether the characteristic
    function created should take vectors
    (each vector corresponding to one domain) as argument,
    or a data frame (with the elements in the rows). In the former case,
    all vectors are recycled to fit the longest vector in case of binary
    relations.}
}
\seealso{
  \code{\link{relation}()}
}
\examples{
## Relation 'a divides b':
divides <- function(a, b) b \%\% a == 0
R <- relation(domain = list(1 : 10, 1 : 10), charfun = divides)
R
## 'Recover' characteristic function:
"\%|\%" <- relation_charfun(R)

## Use it.
2L \%|\% 6L
2:4 \%|\% 6L
2L \%|\% c(2:3, 6L)

## This also works:
"\%|\%"(2L, 6L)
## (and more generally does for arities > 2).
}
\keyword{math}
