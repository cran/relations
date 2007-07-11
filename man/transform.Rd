\name{transform}
\alias{transform_incidences}
\title{Transform incidences}
\description{
  Carry out transformations between incidence matrices from
  endorelations and other codings.
}
\usage{
transform_incidences(x, from = c("PO","SO","01","-1+1"),
                        to = c("PO","SO","01","-1+1"))
}
\arguments{
  \item{x}{An incidence matrix from an endorelation.}
  \item{from, to}{The coding scheme (see details).}
}
\details{
In the following, we consider an incidence matrix \eqn{X} with cells
\eqn{x_{jk}} of a relation \eqn{R} with tuples \eqn{(a_j, b_k)}.

For the \code{"PO"} (\dQuote{Preference Order}) coding,
\eqn{X} is a 0/1 matrix, and 
\eqn{a_j R b_k} iff \eqn{x_{jk} = 1}. It follows in particular
that if both \eqn{x_{jk}} and \eqn{x_{kj}} are 0, the corresponding pair
\eqn{(a_j, b_k)} is not contained in R, i.e., \eqn{a_j} and \eqn{b_k}
are unrelated.

For the \code{"SO"} (\dQuote{"Strict Order"}) coding,
\eqn{X} is a 0/1 matrix with possible
\code{NA} values. As for \code{"PO"}, \eqn{a_j R b_k} iff
\eqn{x_{jk} = 1}, but at most one of \eqn{x_{jk}} and \eqn{x_{kj}} can
be 1. If both are missing (\code{NA}), \eqn{a_j} and \eqn{b_k}
are unrelated.

For the \code{"01"} coding, \eqn{X} is a matrix with values 0, 1, or
0.5. The coding is similar to \code{"SO"}, except that \code{NA} is
represented by 0.5.

For the \code{"-1+1"} coding, \eqn{X} is a matrix with values -1, 0, or 1.
The coding is similar to \code{"SO"}, except that \code{NA} is
represented by 0, and \eqn{x_{jk} = -1} if \emph{not} \eqn{a_j R b_k}.

}
\seealso{
  \code{\link{relation_incidence}}.
}
\examples{
x <- relation(domain = 1:4,
              graph = set(pair(1,2), pair(4,2), pair(1,3), pair(1,4),
                          pair(3,2), pair(2,1)))
inc <- relation_incidence(x)
print(inc)

transform_incidences(inc, to = "SO")
transform_incidences(inc, to = "01")
transform_incidences(inc, to = "-1+1")

## transformations should be loss-free:
inc2 <- transform_incidences(inc, from = "PO", to = "-1+1")
inc2 <- transform_incidences(inc2, from = "-1+1", to = "SO")
inc2 <- transform_incidences(inc2, from = "SO", to = "01")
inc2 <- transform_incidences(inc2, from = "01", to = "PO")
stopifnot(identical(inc, inc2))
}
\keyword{math}
