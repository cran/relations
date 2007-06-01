\name{dissimilarity}
\alias{relation_dissimilarity}
\title{Dissimilarity Between Relations}
\description{Compute the dissimilarity between (ensembles of)
  relations.}
\usage{
relation_dissimilarity(x, y = NULL, method = "symdiff", \dots)
}
\arguments{
  \item{x}{an ensemble of relations, or something which can be coerced
    to such (see \code{\link{relation_ensemble}}).} 
  \item{y}{\code{NULL} (default), or as for \code{x}.}
  \item{method}{a character string specifying one of the built-in
    methods for computing dissimilarity, or a function to be taken as
    a user-defined method.  If a character string, its lower-cased
    version is matched against the lower-cased names of the available
    built-in methods using \code{\link{pmatch}}.  See \bold{Details} for
    available built-in methods.}
  \item{\dots}{further arguments to be passed to methods.}
}
\value{
  If \code{y} is \code{NULL}, an object of class \code{\link{dist}}
  containing the dissimilarities between all pairs of elements of
  \code{x}.  Otherwise, a matrix with the dissimilarities between the
  elements of \code{x} and the elements of \code{y}.
}
\details{
  Currently, \code{"symdiff"} is the only built-in method available.
  This computes the cardinality of the symmetric difference of two
  relations, i.e., the number of tuples contained in exactly one of two
  relations.  For complete endorelations, this coincides with the
  \emph{Kemeny-Snell} metric.
}
\references{
J. G. Kemeny and J. L. Snell (1962),
\emph{Mathematical Models in the Social Sciences},
chapter \dQuote{Preference Rankings: An Axiomatic Approach}.
MIT Press: Cambridge.
}
\keyword{math}
