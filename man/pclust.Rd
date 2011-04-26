\name{pclust}
\encoding{UTF-8}
\alias{relation_pclust}
\title{Prototype-Based Partitions of Relations}
\description{
  Compute prototype-based partitions of a relation ensemble by minimizing
  \eqn{\sum w_b u_{bj}^m d(x_b, p_j)^e}, the sum of the case-weighted and
  membership-weighted \eqn{e}-th powers of the dissimilarities between
  the elements \eqn{x_b} of the ensemble and the prototypes \eqn{p_j},
  for suitable dissimilarities \eqn{d} and exponents \eqn{e}.
}
\usage{
relation_pclust(x, k, method, m = 1, weights = 1,
                control = list())
}
\arguments{
  \item{x}{an ensemble of relations (see
    \code{\link{relation_ensemble}}),
    or something which can be coerced to this.}
  \item{k}{an integer giving the number of classes to be used in the
    partition.}
  \item{method}{the consensus method to be employed, see
    \code{\link{relation_consensus}}.}
  \item{m}{a number not less than 1 controlling the softness of the
    partition (as the \dQuote{fuzzification parameter} of the fuzzy
    \eqn{c}-means algorithm).  The default value of 1 corresponds to
    hard partitions obtained from a generalized \eqn{k}-means problem;
    values greater than one give partitions of increasing softness
    obtained from a generalized fuzzy \eqn{c}-means problem.}
  \item{weights}{a numeric vector of non-negative case weights.
    Recycled to the number of elements in the ensemble given by \code{x}
    if necessary.}
  \item{control}{a list of control parameters.  See \bold{Details}.}
}
\value{
  An object of class \code{\link[clue:partition]{cl_partition}}.
}
\details{
  For \eqn{m = 1}, a generalization of the Lloyd-Forgy variant of the
  \eqn{k}-means algorithm is used, which iterates between reclassifying
  objects to their closest prototypes, and computing new prototypes as
  consensus relations (generalized \dQuote{central relations}, Régnier
  (1965)) for the classes.  This procedure was proposed in Gaul and
  Schader (1988) as the \dQuote{Clusterwise Aggregation of Relations}
  (CAR).

  For \eqn{m > 1}, a generalization of the fuzzy \eqn{c}-means recipe
  is used, which alternates between computing optimal memberships for
  fixed prototypes, and computing new prototypes as the consensus
  relations for the classes.

  This procedure is repeated until convergence occurs, or the maximal
  number of iterations is reached.

  Consensus relations are computed using
  \code{\link{relation_consensus}}.
  
  Available control parameters are as follows.
  \describe{
    \item{\code{maxiter}}{an integer giving the maximal number of
      iterations to be performed.  Defaults to 100.}
    \item{\code{reltol}}{the relative convergence tolerance.  Defaults
      to \code{sqrt(.Machine$double.eps)}.}
    \item{\code{control}}{control parameters to be used in
      \code{\link{relation_consensus}}.}
  }

  The dissimilarities \eqn{d} and exponent \eqn{e} are implied by the
  consensus method employed, and inferred via a registration mechanism
  currently only made available to built-in consensus methods.  For the
  time being, all optimization-based consensus methods use the symmetric
  difference dissimilarity (see \code{\link{relation_dissimilarity}})
  for \eqn{d} and \eqn{e = 1}.

  The fixed point approach employed is a heuristic which cannot be
  guaranteed to find the global minimum.  Standard practice would
  recommend to use the best solution found in \dQuote{sufficiently many}
  replications of the base algorithm.
}
\references{
  S. Régnier (1965).
  Sur quelques aspects mathématiques des problèmes de classification
  automatique.
  \emph{ICC Bulletin}, \bold{4}, 175--191.
  
  W. Gaul and M. Schader (1988).
  Clusterwise aggregation of relations.
  \emph{Applied Stochastic Models and Data Analysis}, \bold{4}, 273--282.
}
\keyword{math}
\keyword{cluster}
