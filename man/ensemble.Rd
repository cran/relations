\name{ensemble}
\alias{relation_ensemble}
\alias{as.relation_ensemble}
\alias{is.relation_ensemble}
\title{Relation Ensembles}
\description{Creation and manipulation of relation ensembles.}
\usage{
relation_ensemble(..., list = NULL)
as.relation_ensemble(x)
is.relation_ensemble(x)
}
\arguments{
  \item{\dots}{\R objects representing relations, or coercible to such.}
  \item{list}{a list of \R objects as in \code{\dots}.}
  \item{x}{for coercion with \code{as.relation_ensemble()}, an \R object
    as in \code{\dots}; for testing with \code{is.relation_ensemble()},
    an arbitrary \R object.}
}
\details{
  \code{relation_ensemble()} creates non-empty \dQuote{relation
    ensembles}, i.e., collections of relations \eqn{R_i = (D, G_i)} with
  the same domain \eqn{D} and possibly different graphs \eqn{G_i}.

  Such ensembles are implemented as suitably classed lists of relation
  objects, making it possible to use \code{lapply()} for computations on
  the individual relations in the ensemble.  Available methods for
  relation ensembles include those for subscripting, \code{c()},
  \code{t()}, \code{rep()}, and \code{print()}.
}
\examples{
data("Cetacea")
## Consider each variable an equivalence relation on the objects.
## Note that 2 variables (LACHRYMAL_AND_JUGAL_BONES and HEAD_BONES) have
## missing values, and hence are excluded.
ind <- sapply(Cetacea, function(s) all(!is.na(s)))
relations <- as.relation_ensemble(Cetacea[, ind])
## This gives a relation ensemble of length 14 (number of complete
## variables in the data set).
print(relations)
## Are there any duplicated relations?
any(duplicated(relations))
## Replicate and combine ...
thrice <- c(rep(relations, 2), relations)
## Extract unique elements again:
all.equal(unique(thrice), relations)
## Note that unique() does not preserve attributes, and hence names.
## In case we want otherwise:
all.equal(thrice[!duplicated(thrice)], relations)
## Subscripting:
relation_dissimilarity(relations[1 : 2], relations["CLASS"])
## Which relation is "closest" to the classification?
d <- relation_dissimilarity(relations)
sort(as.matrix(d)[, "CLASS"])[-1]
}
\keyword{math}
