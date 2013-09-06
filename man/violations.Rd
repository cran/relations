\name{violations}
\encoding{UTF-8}
\alias{relation_violations}
\title{Violations of Relation Properties}
\description{Computes a measure of remoteness of a relation from a
  specified property.}
\usage{
relation_violations(x,
                    property =
                    c("complete", "match",
                      "reflexive", "irreflexive", "coreflexive",
                      "symmetric", "antisymmetric", "asymmetric",
                      "transitive", "negatively_transitive",
                      "Ferrers", "semitransitive",
                      "trichotomous",
                      "Euclidean"),
                    tuples = FALSE,
                    na.rm = FALSE)
}
\arguments{
  \item{x}{an endorelation.}
  \item{property}{a character string specifying one of the properties
    for which the number of violations can be computed.}
  \item{tuples}{a logical indicating whether to return the amount of
    violations (default), or the tuples for which the property is
    violated.}
  \item{na.rm}{a logical indicating whether to remove tuples for which
    the property information is not available (due to missing
    memberships).}
}
\value{
  If \code{tuples} is false (default), the amount of violations for the
  specified property: for crisp relations, the minimum number of object
  tuples involved in the definition of the property (e.g., singletons
  for reflexivity, pairs for antisymmetry, and triples for transitivity)
  that must be modified/added/removed to make the relation satisfy the
  property.

  If \code{tuples} is true, a set of tuples of objects for which the
  respective property is violated. 
}
\seealso{
  \link{predicates} for the definitions of the properties.
}
\examples{
## partial order:
R <- as.relation(1:3)
relation_incidence(R)
## R clearly is transitive, but not symmetric:
relation_violations(R, "transitive")
relation_violations(R, "symmetric")
## Pairs for which symmetry is violated:
relation_violations(R, "symmetric", TRUE)

## create a simple relation:
require("sets")				# set(), pair() etc.
R <- relation(domain = letters[1:2],
              graph = set(pair("a","b"), pair("b","a")))
relation_incidence(R)
## R is clearly symmetric, but not antisymmetric:
relation_violations(R, "symmetric")
relation_violations(R, "antisymmetric")
}
\keyword{math}
