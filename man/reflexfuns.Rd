\name{reflexfuns}
\alias{reflexive_closure}
\alias{reflexive_reduction}
\title{Reflexive Closure and Reduction}
\description{
  Computes reflexive closure and reduction of an endorelation.
}
\usage{
reflexive_reduction(x)
reflexive_closure(x)
}
\arguments{
  \item{x}{an \R object inheriting from class \code{\link{relation}},
    representing an endorelation.}
}
\details{
  Let \eqn{R} be an endorelation on \eqn{X}.
  
  The \emph{reflexive closure} (\emph{reflexive reduction})
  of \eqn{R} is computed by setting the
  diagonal of the incidence matrix to 0 (1).

}
\seealso{
  \code{\link{relation}}, \code{\link{transitive_reduction}}.
}
\examples{
R <- as.relation(1 : 5)
relation_incidence(R)
RR <- reflexive_reduction(R)
relation_incidence(RR)
R == reflexive_closure(RR)
}
\keyword{math}
