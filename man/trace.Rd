\name{trace}
\alias{relation_trace}
\title{Traces of Endorelations}
\description{
  Compute the left or right trace of an endorelation.
}
\usage{
relation_trace(x, which)
}
\arguments{
  \item{x}{an endorelation.}
  \item{which}{one of \code{"left"} or \code{"right"}, or a unique
    abbreviation thereof.}
}
\details{
  Let \eqn{R} be a crisp endorelation.  The left and right trace of
  \eqn{R} contain all pairs \eqn{x, y} for which \eqn{z R x} implies
  \eqn{z R y} for all \eqn{z} (left trace) or \eqn{y R z} implies
  \eqn{x R z} for all \eqn{z} (right trace), respectively.  These are
  the largest (in the natural ordering of relations with the same
  domain) relations such that \eqn{R * S \le R} or \eqn{S * R \le R},
  respectively (where \eqn{*} denotes composition).  In the fuzzy case,
  the memberships of the traces can be defined as the infima over the
  corresponding fuzzy membership implications.  See Chapter 2.3 in Fodor
  and Roubens (1994) for more information.
}
\references{
  J. Fodor and M. Roubens (1994),
  \emph{Fuzzy Preference Modelling and Multicriteria Decision Support}.
  Kluwer Academic Publishers: Dordrecht.
}
\keyword{math}
