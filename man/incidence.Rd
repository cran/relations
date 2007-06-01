\name{incidence}
\alias{relation_incidence}
\title{Relation Incidences}
\description{
  Determine the incidences of a relation.
}
\usage{
relation_incidence(x)
}
\arguments{
  \item{x}{an object inheriting from class \code{\link{relation}}.}
}
\value{
  For a \eqn{k}-ary relation, a \eqn{k}-dimensional 0/1 array inherting
  from class \code{relation_incidence} whose elements indicate whether
  the corresponding \eqn{k}-tuples are contained in the relation or not.
}
\seealso{
  \code{\link{relation}};
  \code{\link{relation_incidence<-}} for modifying the incidences.
}
\examples{
R <- as.relation(c(A = 1, B = 2, C = 3))
relation_incidence(R)
}
\keyword{math}
  
