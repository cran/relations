\name{scores}
\encoding{UTF-8}
\alias{relation_scores}
\alias{relation_scores.relation}
\alias{relation_scores.relation_ensemble}
\title{Relation Scores}
\description{
  Compute scores for the tuples of (ensembles of) endorelations.
}
\usage{
\method{relation_scores}{relation}(x,
                method = c("ranks", "Barthelemy/Monjardet",
                           "Borda", "Kendall", "Wei",
                           "differential", "Copeland"),
                normalize = FALSE, \dots)
\method{relation_scores}{relation_ensemble}(x,
                method = c("Borda", "Kendall", "differential",
                           "Copeland"),
                normalize = FALSE,
                weights = 1, \dots)
}
\arguments{
  \item{x}{an object inheriting from class \code{\link{relation}},
    representing an endorelation.}
  \item{method}{character string indicating the method (see
    \bold{Details}).}
  \item{normalize}{logical indicating whether the score vector should be
    normalized to sum up to 1.}
  \item{weights}{Numeric vector of weights used in incidence
    aggregation, recycled as needed.}
  \item{\dots}{further arguments to be passed to methods.}
}
\details{
  In the following, consider an endorelation \eqn{R} on \eqn{n} objects.
  Let the \emph{in-degree} \eqn{I(x)} and \emph{out-degree} \eqn{O(x)}
  of an object \eqn{x} be defined as the numbers of objects \eqn{y} such
  that \eqn{y R x} and, respectively, \eqn{x R y}, and let \eqn{D(x) =
    I(x) - O(x)} be the \emph{differential} of \eqn{x} (see Regenwetter
  and Rykhlevskaia (2004)).  Note that \eqn{I} and \eqn{O} are given by
  the column sums and row sums of the incidence matrix of \eqn{R}.  If
  \eqn{R} is a preference relation with a \eqn{\le} interpretation,
  \eqn{D(x)} is the difference between the numbers of objects dominated
  by \eqn{x} (i.e., \eqn{< x}) and dominating \eqn{x} (i.e., \eqn{> x}),
  as \dQuote{ties} cancel out.

  \code{relation_score()} is generic with methods for relations and
  relation ensembles. Available built-in score methods for the
  relation method are as follows:
  
  \describe{
    \item{\code{"ranks"}}{generalized ranks.  A linear transformation of
      the differential \eqn{D} to the range from 1 to \eqn{n}.  An
      additional argument \code{decreasing} can be used to specify the
      order of the ranks.  By default, or if \code{decreasing} is true,
      objects are ranked according to decreasing differential
      (\dQuote{from the largest to the smallest} in the \eqn{\le}
      preference context) using \eqn{(n + 1 - D(x)) / 2}.
      Otherwise, if \code{decreasing} is false, objects are ranked via
      \eqn{(n + 1 + D(x)) / 2} (\dQuote{from the smallest to the largest}).
      See Regenwetter and Rykhlevskaia (2004) for more details on
      generalized ranks.
    }
    \item{\code{"Barthelemy/Monjardet"}}{\eqn{(M(x) + N(x) - 1) / 2},
      where \eqn{M(x)} and \eqn{N(x)} are the numbers of objects \eqn{y}
      such that \eqn{y R x}, and \eqn{y R x} and not \eqn{x R y},
      respectively.  If \eqn{R} is a \eqn{\le} preference relation, we
      get the number of dominated objects plus half the number of the
      equivalent objects minus 1 (the \dQuote{usual} average ranks minus
      one if the relation is complete).  See Barthélemy and Monjardet
      (1981).
    }
    \item{\code{"Borda"}, \code{"Kendall"}}{the out-degrees.
      See Borda (1770) and Kendall (1955).
    }
    \item{\code{"Wei"}}{the eigenvector corresponding to the greatest
      eigenvalue of the incidence matrix of the complement of \eqn{R}.
      See Wei (1952).
    }
    \item{\code{"differential"}, \code{"Copeland"}}{the differentials,
    equivalent to the negative \emph{net flow} of Bouyssou (1992), and
    also to the Copeland scores.}
  }

  For relation ensembles, currently only
  \code{"differential"}/\code{"Copeland"} and
  \code{"Borda"}/\code{"Kendall"} are implemented. They are computed on
  the aggregated incidences of the ensembles' relations.
  
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
  \emph{Mathematical Social Sciences}, \bold{1}, 235--267.
  \doi{10.1016/0165-4896(81)90041-X}.

  J. C. Borda (1781),
  Mémoire sur les élections au scrutin.
  Histoire de l'Académie Royale des Sciences.

  D. Bouyssou (1992),
  Ranking methods based on valued preference relations: A
  characterization of the net flow network.
  \emph{European Journal of Operational Research}, \bold{60}, 61--67.
  \doi{10.1016/0377-2217(92)90333-5}.
	
  M. Kendall (1955),
  Further contributions to the theory of paired comparisons.
  \emph{Biometrics}, \bold{11}, 43--62.
  \doi{10.2307/3001479}.

  M. Regenwetter and E. Rykhlevskaia (2004),
  On the (numerical) ranking associated with any finite binary relation.
  \emph{Journal of Mathematical Psychology}, \bold{48}, 239--246.
  \doi{10.1016/j.jmp.2004.03.003}.
  
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

## Example taken from Cook and Cress (1992, p.136)
## Note that the results indicated for the Copeland scores have
## (erroneously?) been computed from the *unweighted* votes.
## Also, they report the votes as strict preferences, so we
## create the dual relations.

D <- letters[1:5]
X <- as.relation(ordered(D, levels = c("b", "c", "a", "d", "e")))
Y <- as.relation(ordered(D, levels = c("d", "a", "e", "c", "b")))
Z <- as.relation(ordered(D, levels = c("e", "c", "b", "a", "d")))
E <- relation_ensemble(X, Y, Z)
relation_scores(E, "Copeland")
relation_scores(E, "Borda", weights = c(4, 3, 2))
}
\keyword{math}
