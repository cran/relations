\name{cover}
\alias{relation_cover}
\title{Covering Relations}
\description{
  Compute the covering relation of an endorelation.
}
\usage{
relation_cover(x)
}
\arguments{
  \item{x}{an endorelation.}
}
\details{
  Let \eqn{R} be an endorelation with domain \eqn{(X, X)} and \eqn{P} be
  the asymmetric part of \eqn{R} for which \eqn{x P y} iff \eqn{x R y}
  and not \eqn{y R x}.  (If \eqn{R} is a \eqn{\le} order relation,
  \eqn{P} is the associated strict order.)  We say that \eqn{x} is
  covered by \eqn{y} if \eqn{x P y} and there is no \eqn{z} such that
  \eqn{x P z} and \eqn{z P y}.  One also says that \eqn{y} covers
  \eqn{x}, or that it is a successor of \eqn{x}.

  The covering relation of \eqn{R} consists of all pairs \eqn{(x, y)}
  for which \eqn{x} is covered by \eqn{y}.
}
\keyword{math}
