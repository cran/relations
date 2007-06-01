\name{transfuns}
\alias{transitive_closure}
\alias{transitive_reduction}
\title{Transitive Closure and Reduction}
\description{
  Computes transitive closure and reduction of an endorelation.
}
\usage{
transitive_reduction(x)
transitive_closure(x)
}
\arguments{
  \item{x}{an \R object inheriting from class \code{\link{relation}},
    representing an endorelation.}
}
\details{
  Let \eqn{R} be an endorelation on \eqn{X} and \eqn{n} be the number of
  elements in \eqn{X}.
  
  The \emph{transitive closure} of \eqn{R} is the smallest transitive
  relation on \eqn{X} that contains \eqn{R}.  The code implements
  Warshall's Algorithm which is of complexity \eqn{O(n^3)}.

  The \emph{transitive reduction} of \eqn{R} is the smallest relation
  \eqn{R'} on \eqn{X} so that the the transitive closure of \eqn{R'} is
  the same than the transitive closure of \eqn{R}.  The function is
  implemented using a depth-first-search approach, also with complexity
  \eqn{O(n^3)}.
}
\seealso{
  \code{\link{relation}}, \code{\link{reflexive_reduction}}.
}
\examples{
R <- as.relation(1 : 5)
relation_incidence(R)
RR <- transitive_reduction(R)
relation_incidence(RR)
R == transitive_closure(RR)
}
\references{
  Warshall, Stephen (1962).
  A theorem on Boolean matrices.
  \emph{Journal of the ACM}, \bold{9(1)}, 11--12.
}
\keyword{math}
