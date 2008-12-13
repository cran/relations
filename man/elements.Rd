\name{elements}
\alias{relation_elements}
\title{Elements of Relation Domains}
\description{
  Obtain elements of endorelation domains which have certain properties.
}
\usage{
relation_elements(x, which, ...)
}
\arguments{
  \item{x}{an endorelation.}
  \item{which}{a character string specifying the property to be tested
    for.  Currently, one of \code{"minimal"}, \code{"first"},
    \code{"last"}, or \code{"maximal"}, or a unique abbreviation
    thereof.}
  \item{...}{additional arguments to be employed in the property tests.}
}
\details{
  Let \eqn{R} be an endorelation with domain \eqn{(X, X)} and consider
  elements \eqn{x} and \eqn{y} of \eqn{X}.  We say that \eqn{x} is
  \describe{
    \item{minimal:}{there is no \eqn{y \ne x}{y != x} with \eqn{y R x}.}
    \item{a first element:}{\eqn{x R y} for all \eqn{y \ne x}{y != x}.}
    \item{a last element:}{\eqn{y R x} for all \eqn{y \ne x}{y != x}.}
    \item{maximal:}{there is no \eqn{y \ne x}{y != x} with \eqn{x R y}.}
  }
  When computing the tests for the above properties, an additional
  \code{na.rm} argument can be given to control the handling of missing
  incidences.  By default, these are treated as false, to the effect
  that they invalidate \dQuote{for all} tests (corresponding to
  \code{na.rm = FALSE}) and pass the \dQuote{there is no} tests
  (corresponding to \code{na.rm = TRUE}).
}
\value{
  A set with the elements having the specified property.
}
\keyword{math}
