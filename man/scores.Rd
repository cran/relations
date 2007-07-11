\name{scores}
\encoding{UTF-8}
\alias{relation_scores}
\title{Relation Scores}
\description{
  Compute scores for the tuples of an endorelation.
}
\usage{
relation_scores(x,
                method = c("Barthelemy/Monjardet", "Wei", "Borda",
                           "Kendall", "differential", "rankB"),
                normalize = FALSE)
}
\arguments{
  \item{x}{an object inheriting from class \code{\link{relation}},
    representing an endorelation.}
  \item{method}{character string indicating the method (see details).}
  \item{normalize}{logical indicating whether the score vector should be
    normalized to sum up to 1.}
}
\details{
  In the following, consider the pairs of relation \eqn{R}, represented
  by \code{x}.  Available built-in methods are as follows:
  
  \describe{
    \item{\code{"Barthelemy/Monjardet"}}{for each \eqn{a} in the domain
      of \eqn{R}, \eqn{(M(a) + N(a) - 1) / 2}, where \eqn{M(a)} and
      \eqn{N(a)} are the numbers of objects \eqn{b} such that
      \eqn{b R a}, and \eqn{b R a} and not \eqn{a R b}, respectively.
      If \eqn{R} is a \eqn{\le} preference relation, we get the number
      of dominated objects plus half the number of the equivalent
      objects minus 1, i.e., the usual average ranks minus 1.
      See Barthélemy and Monjardet (1981).
    }
    \item{\code{"Wei"}}{the eigenvector corresponding to the greatest
      eigenvalue of the incidence matrix of the complement of \eqn{R}.
      See Wei (1952).
    }
    \item{\code{"Kendall"}, \code{"Borda"}}{the column sums of the
      incidence matrix, i.e., for each object \eqn{a}, the number of
      objects \eqn{b} such that \eqn{b R a}.
      See Borda (1770) and Kendall (1955).
    }
    \item{\code{"differential"}}{the column sums of the incidence
      matrix minus the row sums.  For each object \eqn{a}, this is the 
      difference of the numbers of objects \eqn{b} such that \eqn{b R a}
      and \eqn{a R b}, respectively.  In the case of a \eqn{\le}
      preference relation, this is the difference between the
      \emph{indegree} and the \emph{outdegree} of each object/node in
      the digraph corresponding to the relation incidence.
      See Regenwetter and Rykhlevskaia (2004).
    }
    \item{\code{"rankB"}}{a linear transformation of the differential
      \eqn{D(R)}, so that the scores of the \eqn{n} objects
      start with 1: \eqn{(n + 1 + D(R)) / 2}.}
  }
}
\value{
  A vector of scores, with names taken from the relation domain labels.
}
\references{
  J.-P. Barthélemy and B. Monjardet (1981),
  The median procedure in cluster analysis and social choice theory.
  \emph{Mathematical Social Sciences}, \bold{1}:235--267.

  J. C. Borda (1781),
  Mémoire sur les élections au scrutin.
  Histoire de l'Académie Royale des Sciences.

  M. Kendall (1955),
  Further contributions to the theory of paired comparisons.
  \emph{Biometrics}, \bold{11}:43--62.

  M. Regenwetter and E. Rykhlevskaia (2004),
  On the (numerical) ranking associated with any finite binary relation.
  \emph{Journal of Mathematical Psychology}, \bold{48}:239--246.
  
  T. H. Wei (1952).
  \emph{The algebraic foundation of ranking theory}.
  Unpublished thesis, Cambridge University.
}
\examples{
## Example taken from Cook and Cress (1992, p.74)
I <- matrix(c(0, 0, 1, 1, 1,
              1, 0, 0, 0, 1,
              0, 1, 0, 0, 1,
              0, 1, 1, 0, 0,
              0, 0, 0, 1, 0),
            ncol = 5,
            byrow = TRUE)
R <- relation(domain = letters[1:5], incidence = I)

## Note that this is a "preference matrix", so take complement:
R <- !R

## Compare Kendall and Wei scores
cbind(
      Kendall = relation_scores(R, method = "Kendall", normalize = TRUE),
      Wei = relation_scores(R, method = "Wei", normalize = TRUE)
     )

}
\keyword{math}
