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
  same than the transitive closure of \eqn{R}.

  The transitive reduction of an \emph{acyclic} relation can be obtained
  by subtracting from \eqn{R} the composition of \eqn{R} with its transitive closure.

  The transitive reduction of a \emph{cyclic} relation is the transitive
  reduction of the condensation, combined with the component
  representation of \eqn{R}. (Note that the transitive reduction of a
  cyclic relation is cyclic.)
  
  The \emph{reflexive reduction} of \eqn{R} is computed by setting the
  diagonal of the incidence matrix to 0.

}
\seealso{
  \code{\link{relation}()},\cr
  \code{\link{reflexive_reduction}()},
  \code{\link{transitive_reduction}()},
  \code{\link[sets]{reduction}()},\cr
  \code{\link{relation_condensation}()},
  \code{\link{relation_component_representation}()}.
}
\examples{
R <- as.relation(1 : 5)
relation_incidence(R)

## transitive closure/reduction
RR <- transitive_reduction(R)
relation_incidence(RR)
R == transitive_closure(RR)

## same
require("sets")				# closure() and reduction() etc.
R == closure(reduction(R))

## reflexive closure/reduction

RR <- reflexive_reduction(R)
relation_incidence(RR)
R == reflexive_closure(RR)

## same:
R == closure(reduction(R, "reflexive"), "reflexive")

## transitive reduction of a cyclic relation:
## (example from La Poutre and van Leeuwen)

require("sets")				# set(), pair() etc.
if(require("Rgraphviz")) {
  G <- set(pair(1L, 2L), pair(2L, 1L), pair(1L, 3L), pair(3L, 1L),
           pair(3L, 7L), pair(2L, 5L), pair(2L, 6L), pair(6L, 5L),
           pair(5L, 7L), pair(4L, 6L), pair(5L, 4L), pair(4L, 7L))
  R <- endorelation(graph = G)
  plot(relation_ensemble(R, R), type = c("raw", "simplified"), main =
           c("original graph", "transitive reduction"))
  }
}
\references{
  \bibshow{La_Poutré+van_Leeuwen:1988, Warshall:1962}
}
\keyword{math}
