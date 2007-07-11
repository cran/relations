\name{violations}
\encoding{UTF-8}
\alias{relation_violations}
\title{Violations of Relation Properties}
\description{Computes a measure of remoteness of a relation from a
  specified property.}
\usage{
relation_violations(x, family = c("T", "transitive",
                                  "S", "symmetric",
                                  "A", "antisymmetric",
                                  "C", "complete",
                                  "R", "reflexive"))

}
\arguments{
  \item{x}{an endorelation.}
  \item{family}{a character string specifying one of the
    relation families for which the number of violations can be
    computed.}
}
\value{
  The number of violations for the specified family, i.e. the
  minimum number of pairs (for transitive relations: triples)
  that must be modified/added/removed to make
  the relation comply to the families' restrictions.
}
\examples{
## partial order:
R <- as.relation(1:3)
relation_incidence(R)
## R clearly is transitive, but not symmetric:
relation_violations(R, "transitive")
relation_violations(R, "symmetric")

## create a simple relation:
R <- relation(domain = letters[1:2],
              graph = set(pair("a","b"), pair("b","a")))
relation_incidence(R)
## R is clearly symmetric, but not antisymmetric:
relation_violations(R, "symmetric")
relation_violations(R, "antisymmetric")
}
\keyword{math}
