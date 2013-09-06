\name{closure}
\alias{transitive_closure}
\alias{reflexive_closure}
\alias{closure.relation}
\title{Transitive and Reflexive Closure}
\description{
  Computes transitive and reflexive closure of an endorelation.
}
\usage{
transitive_closure(x)
reflexive_closure(x)
\method{closure}{relation}(x, operation = c("transitive", "reflexive"), ...)
}
\arguments{
  \item{x}{an \R object inheriting from class \code{\link{relation}},
    representing an endorelation.}
  \item{operation}{character string indicating the kind of closure.}
  \item{\dots}{currently not used.}
}
\details{
  Let \eqn{R} be an endorelation on \eqn{X} and \eqn{n} be the number of
  elements in \eqn{X}.
  
  The \emph{transitive closure} of \eqn{R} is the smallest transitive
  relation on \eqn{X} that contains \eqn{R}.  The code implements
  Warshall's Algorithm which is of complexity \eqn{O(n^3)}.

  The \emph{reflexive closure} of \eqn{R} is computed by setting the
  diagonal of the incidence matrix to 1.

}
\seealso{
  \code{\link{relation}},
  \code{\link{reflexive_reduction}},
  \code{\link{transitive_reduction}},
  \code{\link[sets]{closure}}.
}
\examples{
R <- as.relation(1 : 5)
relation_incidence(R)

## transitive closure/reduction
RR <- transitive_reduction(R)
relation_incidence(RR)
R == transitive_closure(RR)

## same
require("sets")				# closure() and reduction()
R == closure(reduction(R))

## reflexive closure/reduction

RR <- reflexive_reduction(R)
relation_incidence(RR)
R == reflexive_closure(RR)
## same:
R == closure(reduction(R, "reflexive"), "reflexive")
}
\references{
  S. Warshall (1962),
  A theorem on Boolean matrices.
  \emph{Journal of the ACM}, \bold{9}/1, 11--12.
}
\keyword{math}
