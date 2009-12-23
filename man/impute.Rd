\name{impute}
\alias{relation_impute}
\title{Impute relations}
\description{
  Impute missing incidences in relations by averaging all possible
  completions within a specified family.
}
\usage{
relation_impute(x, method = NULL, control = list(), \dots)
}
\arguments{
  \item{x}{an endorelation or an ensemble of endorelations.}
  \item{method}{character string specifying the method to be used (see
    \bold{Details}).  If \code{NULL}, it is guessed from the
    relation with missing \emph{objects} removed.}
  \item{control}{named list of control arguments. Currently, only
    \code{n} is accepted by the \code{any/\var{F}} methods, indicating
    the number of solutions to be returned. Default is 1; \code{"all"}
    or \code{NA} will generate all possible completions. Note that
    \code{n} is currently ignored if \code{x} is a relation ensemble.}
  \item{\dots}{named list of control arguments, overriding the ones
    in \code{control}.}

}
\value{
  If \code{x} is an ensemble or more than one solution is requested
  using the \code{n} control argument: an ensemble
  of endorelations.  An endorelation otherwise.
}
\details{
  Endorelations with missing elements (i.e., whose incidence is
  \code{NA}) are imputed using one of the methods described as follows.

  \describe{
    \item{\code{"omit"}}{Missing incidences are replaced by zeros, i.e.,
      the corresponding elements are removed from the graph.}
    \item{\code{"any/\var{F}"}}{The incidences are replaced by
      arbitrary values suitable for family \var{F}, with possible values:
      \describe{
	\item{\code{G}}{General (unrestricted) relations.}
	\item{\code{L}}{Linear orders.}
	\item{\code{W}}{Weak orders.}
	\item{\code{O}}{Partial orders.}
      }
      \code{L}, \code{W}, and \code{O} can optionally be complemented by
      \code{/first} or \code{/last} to further restrict imputed elements
      to be placed on top or bottom of the given ordering.
    }
    \item{\code{"average/\var{F}"}}{Computes the relation with average
      incidences, based on all possible completions as indicated for the
      \code{any/\var{F}} methods. Note that these completions are not
      explicitely generated to compute the averages, and that the
      resulting relation will typically be fuzzy.}
    }
}
\examples{
## create a relation with a missing object
R <- ranking(1:2, 1:3)
print(R)
R <- as.relation(R)

## find all suitable completions within L
plot(relation_impute(R, method = "any/L", n = "all"))

## find 3 suitable partial orders
plot(relation_impute(R, method = "any/O", n = 3))

## compute average completion
R1 <- relation_impute(R, method = "average/O")
relation_incidence(R1)

## check correctness of averaging
R2 <- mean(relation_impute(R, "any/O", n = "all"))
stopifnot(all.equal(R1, R2))
}
\keyword{math}
