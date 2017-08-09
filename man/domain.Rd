\name{domain}
\alias{relation_arity}
\alias{relation_domain}
\alias{relation_domain_names}
\alias{relation_size}
\title{Relation Domain, Arity, and Size}
\description{
  Determine the domain, domain names, arity, or size of a relation or a
  relation ensemble.
}
\usage{
relation_arity(x)
relation_domain(x)
relation_domain_names(x)
relation_size(x)
}
\arguments{
  \item{x}{an \R object inheriting from class \code{\link{relation}} or
    \code{\link{relation_ensemble}}.}
}
\value{
  For determining the domain, an object of class \code{relation_domain},
  inheriting from \code{\link[sets]{tuple}}.
}
\seealso{
  \code{\link[sets]{tuple}()};
  \code{\link{relation}()};
  \code{\link{relation_domain<-}()} and
  \code{\link{relation_domain_names<-}()}
  for modifying the domain and domain names of a relation, respectively.
}
\examples{
## A simple relation:
R <- as.relation(c(A = 1, B = 2, C = 3))
relation_incidence(R)
relation_arity(R)
relation_domain(R)
relation_domain_names(R)
relation_size(R)
}
\keyword{math}
