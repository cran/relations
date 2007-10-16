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
  For a \eqn{k}-ary relation, a \eqn{k}-dimensional numeric array with
  values in the unit interval inherting from class
  \code{relation_incidence} whose elements give the memberships of the
  corresponding \eqn{k}-tuples are contained in the relation (for a
  crisp relation, a binary (0/1) array with elements indicating whether
  the corresponding tuples are contained in the relation or not).
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
  
