\name{reduction}
\alias{transitive_reduction}
\alias{reflexive_reduction}
\alias{reduction.relation}
\title{Transitive and Reflexive Reduction}
\description{
  Computes transitive and reflexive reduction of an endorelation.
}
\usage{
transitive_reduction(x)
reflexive_reduction(x)
\method{reduction}{relation}(x, operation = c("transitive", "reflexive"), ...)
}
\arguments{
  \item{x}{an \R object inheriting from class \code{\link{relation}},
    representing an endorelation.}
  \item{operation}{character string indicating the kind of reduction.}
  \item{\dots}{currently not used.}
}
\details{
  Let \eqn{R} be an endorelation on \eqn{X} and \eqn{n} be the number of
  elements in \eqn{X}.
  
  The \emph{transitive reduction} of \eqn{R} is the smallest relation
  \eqn{R'} on \eqn{X} so that the transitive closure of \eqn{R'} is the
  same than the transitive closure of \eqn{R}.  The function is
  implemented using a depth-first-search approach with complexity
  \eqn{O(n^3)}.  Currently, it can only be used for crisp relations.

  The \emph{reflexive reduction} of \eqn{R} is computed by setting the
  diagonal of the incidence matrix to 0.

}
\seealso{
  \code{\link{relation}},
  \code{\link{reflexive_reduction}},
  \code{\link{transitive_reduction}},
  \code{\link[sets:closure]{reduction}}.
}
\examples{
R <- as.relation(1 : 5)
relation_incidence(R)

## transitive closure/reduction
RR <- transitive_reduction(R)
relation_incidence(RR)
R == transitive_closure(RR)

## same
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
