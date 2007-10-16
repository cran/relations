\name{scores}
\encoding{UTF-8}
\alias{relation_scores}
\title{Relation Scores}
\description{
  Compute scores for the tuples of an endorelation.
}
\usage{
relation_scores(x,
                method = c("ranks", "Barthelemy/Monjardet", "Borda",
                           "Kendall", "Wei", "differential"),
                normalize = FALSE, \dots)
}
\arguments{
  \item{x}{an object inheriting from class \code{\link{relation}},
    representing an endorelation.}
  \item{method}{character string indicating the method (see details).}
  \item{normalize}{logical indicating whether the score vector should be
    normalized to sum up to 1.}
  \item{\dots}{further arguments to be passed to methods.}
}
\details{
  In the following, consider an endorelation \eqn{R} on \eqn{n} objects.
  Let the \emph{in-degree} \eqn{I(a)} and \emph{out-degree} \eqn{O(a)}
  of an object \eqn{a} be defined as the numbers of objects \eqn{b} such
  that \eqn{b R a} and, respectively, \eqn{a R b}, and let \eqn{D(a) =
    I(a) - O(a)} be the \emph{differential} of \eqn{a} (see Regenwetter
  and Rykhlevskaia (2004)).  Note that \eqn{I} and \eqn{O} are given by
  the column sums and row sums of the incidence matrix of \eqn{R}.  If
  \eqn{R} is a preference relation with a \eqn{\le} intepretation,
  \eqn{D(a)} is the difference between the numbers of objects dominated
  by \eqn{a} (i.e., \eqn{< a}) and dominating \eqn{a} (i.e., \eqn{> a}),
  as ``ties'' cancel out.

  Available built-in score methods are as follows:
  
  \describe{
    \item{\code{"ranks"}}{generalized ranks.  A linear transformation of
      the differential \eqn{D} to the range from 1 to \eqn{n}.  An
      additional argument \code{decreasing} can be used to specify the
      order of the ranks.  By default, or if \code{decreasing} is true,
      objects are ranked according to decreasing differential (``from
      the largest to the smallest in the \eqn{\le} preference context)
      using \eqn{(n + 1 - D(a)) / 2}.  Otherwise, if \code{decreasing}
      is false, objects are ranked via \eqn{(n + 1 + D(a)) / 2} (``from
      the smallest to the largest).  See Regenwetter and Rykhlevskaia
      (2004) for more details on generalized ranks.
    }
    \item{\code{"Barthelemy/Monjardet"}}{\eqn{(M(a) + N(a) - 1) / 2},
      where \eqn{M(a)} and \eqn{N(a)} are the numbers of objects \eqn{b}
      such that \eqn{b R a}, and \eqn{b R a} and not \eqn{a R b},
      respectively.  If \eqn{R} is a \eqn{\le} preference relation, we
      get the number of dominated objects plus half the number of the
      equivalent objects minus 1 (the ``usual'' average ranks minus one
      if the relation is complete).  See Barthélemy and Monjardet
      (1981).
    }
    \item{\code{"Borda"}, \code{"Kendall"}}{the out-degrees.
      See Borda (1770) and Kendall (1955).
    }
    \item{\code{"Wei"}}{the eigenvector corresponding to the greatest
      eigenvalue of the incidence matrix of the complement of \eqn{R}.
      See Wei (1952).
    }
    \item{\code{"differential"}}{the differentials, equivalent to the
      negative \emph{net flow} of Bouyssou (1992).}
  }

  Definitions of scores for \dQuote{preference relations} \eqn{R} are
  somewhat ambiguous because \eqn{R} can encode \eqn{\le} or \eqn{\ge}
  (or strict variants thereof) relationships (and all such variants are
  used in the literature).  Package \pkg{relations} generally assumes a
  \eqn{\le} encoding, and that scores in the strict sense should
  increase with preference (the most preferred get the highest scores)
  whereas ranks decrease with preference (the most preferred get the
  lowest ranks).
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

  D. Bouyssou (1992),
  Ranking methods based on valued preference relations: A
  characterization of the net flow network.
  \emph{European Journal of Operational Research}, \bold{60}:61--67.
	
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
