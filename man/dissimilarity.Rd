\name{dissimilarity}
\alias{relation_dissimilarity}
\title{Dissimilarity Between Relations}
\description{Compute the dissimilarity between (ensembles of)
  relations.}
\usage{
relation_dissimilarity(x, y = NULL, method = "symdiff", \dots)
}
\arguments{
  \item{x}{an ensemble of relations (see
    \code{\link{relation_ensemble}()}),
    or something which can be coerced to such.} 
  \item{y}{\code{NULL} (default), or as for \code{x}.}
  \item{method}{a character string specifying one of the built-in
    methods for computing dissimilarity, or a function to be taken as
    a user-defined method.  If a character string, its lower-cased
    version is matched against the lower-cased names of the available
    built-in methods using \code{\link{pmatch}()}.
    See \bold{Details} for available built-in methods.}
  \item{\dots}{further arguments to be passed to methods.}
}
\value{
  If \code{y} is \code{NULL}, an object of class \code{\link{dist}}
  containing the dissimilarities between all pairs of elements of
  \code{x}.  Otherwise, a matrix with the dissimilarities between the
  elements of \code{x} and the elements of \code{y}.
}
\details{
  Available built-in methods are as follows.

  \describe{
    \item{\code{"symdiff"}}{symmetric difference distance.
      This computes the cardinality of the symmetric difference of two
      relations, i.e., the number of tuples contained in exactly one of
      two relations.  For preference relations, this coincides with the
      \emph{Kemeny-Snell} metric (Kemeny and Snell, 1962).  For linear
      orders, it gives Kendall's \eqn{\tau} metric (Diaconis, 1988).
      
      Can also be referred to as \code{"SD"}.

      Only applicable to crisp relations.
    }

    \item{\code{"manhattan"}}{the Manhattan distance between the
      incidences.
    }

    \item{\code{"euclidean"}}{the Euclidean distance between the
      incidences.
    }

    \item{\code{"CS"}}{Cook-Seiford distance, a generalization of the
      distance function of Cook and Seiford (1978).  Let the generalized
      ranks of an object \eqn{a} in the (first) domain of an
      endorelation \eqn{R} be defined as the number of objects \eqn{b}
      dominating \eqn{a} (i.e., for which \eqn{a R b} and not \eqn{b R
      a}), plus half the number of objects \eqn{b} equivalent to \eqn{a}
      (i.e., for which \eqn{a R b} and \eqn{b R a}).  For preference
      relations, this gives the usual Kendall ranks arranged according
      to decreasing preference (and averaged for ties).  Then the
      generalized Cook-Seiford distance is defined as the \eqn{l_1}
      distance between the generalized ranks.  For linear orders, this
      gives Spearman's footrule metric (Diaconis, 1988).

      Only applicable to crisp endorelations.
    }

    \item{\code{"CKS"}}{Cook-Kress-Seiford distance, a generalization of
      the distance function of Cook, Kress and Seiford (1986).  For each
      pair of objects \eqn{a} and \eqn{b} in an endorelation \eqn{R}, we
      can have \eqn{a R b} and not \eqn{b R a} or vice versa (cases of
      \dQuote{strict preference}), \eqn{a R b} and \eqn{b R a} (the case
      of \dQuote{indifference}), or neither \eqn{a R b} nor \eqn{b R a}
      (the case of \dQuote{incomparability}).  (Only the last two are
      possible if \eqn{a = b}.)  The distance by Cook, Kress and Seiford
      puts indifference as the metric centroid between both preference
      cases and incomparability (i.e., indifference is at distance one
      from the other three, and each of the other three is at distance
      two from the others).  The generalized Cook-Kress-Seiford distance
      is the paired comparison distance (i.e., a metric) based on these
      distances between the four paired comparison cases.  (Formula 3 in
      the reference must be slightly modified for the generalization
      from partial rankings to arbitrary endorelations.)

      Only applicable to crisp endorelations.
    }

    \item{\code{"score"}}{score-based distance.  This computes
      \eqn{\Delta(s(x), s(y))} for suitable score and distance functions
      \eqn{s} and \eqn{\Delta}, respectively.  These can be specified by
      additional arguments \code{score} and \code{Delta}.  If
      \code{score} is a character string, it is taken as the method for
      \code{\link{relation_scores}}.  Otherwise, if given it must be a
      function giving the score function itself.  If \code{Delta} is a
      number \eqn{p \ge 1}, the usual \eqn{l_p} distance is used.
      Otherwise, it must be a function giving the distance function.
      The defaults correspond to using the default relation scores and
      \eqn{p = 1}, which for linear orders gives Spearman's footrule
      distance.

      Only applicable to endorelations.
    }

    \item{\code{"Jaccard"}}{Jaccard distance: 1 minus the ratio of the
      cardinalities of the intersection and the union of the relations.
    }

    \item{\code{"PC"}}{(generalized) paired comparison distance.  This
      generalizes the symdiff and CKS distances to use a general set of
      discrepancies \eqn{\delta_{kl}} between the possible paired
      comparison results with \eqn{a,b}/\eqn{b,a} incidences 0/0, 1/0,
      0/1, and 1/1 numbered from 1 to 4 (in a preference context with a
      \eqn{\le} encoding, these correspond to incompatibility, strict
      \eqn{<} and \eqn{>} preference, and indifference), with
      \eqn{\delta_{kl}} the discrepancy between possible results \eqn{k}
      and \eqn{l}.  The distance is then obtained as the sum of the
      discrepancies from the paired comparisons of distinct objects,
      plus half the sum of discrepancies from the comparisons of
      identical objects (for which the only possible results are
      incomparability and indifference). 
      The distance is a metric provided that the \eqn{\delta_{kl}}
      satisfy the metric conditions (non-negativity and zero iff
      \eqn{k = l}, symmetry and sub-additivity).

      The discrepancies can be specified via the additional argument
      \code{delta}, either as a numeric vector of length 6 with the
      non-redundant values \eqn{\delta_{21}, \delta_{31}, \delta_{41},
      \delta_{32}, \delta_{42}, \delta_{43}}, or as a character string
      partially matching one of the following built-in discrepancies
      with corresponding parameter vector \eqn{\delta}:

      \describe{
	\item{\code{"symdiff"}}{symmetric difference distance, with
	  discrepancy between distinct results two between either
	  opposite strict preferences or indifference and
	  incomparability, and one otherwise:
	  \eqn{\delta = (1, 1, 2, 2, 1, 1)}
	  (default).
	  
	  Can also be referred to as \code{"SD"}.	  
	}
	\item{\code{"CKS"}}{Cook-Kress-Seiford distance, see above:
	  \eqn{\delta = (2, 2, 1, 2, 1, 1)}.}
	\item{\code{"EM"}}{the distance obtained from the generalization
	  of the Kemeny-Snell distance for complete rankings to partial
	  rankings introduced in Emond and Mason (2000).  This uses a
	  discrepancy of two for opposite strict preferences, and one
	  for all other distinct results: 
	  \eqn{\delta = (1, 1, 1, 2, 1, 1)}.}
	\item{\code{"JMB"}}{the distance with parameters as suggested by
	  Jabeur, Martel and Ben Khélifa (2004):
	  \eqn{\delta = (4/3, 4/3, 4/3, 5/3, 1, 1)}.}
	\item{\code{"discrete"}}{the discrete metric on the set of paired
	  comparison results:
	  \eqn{\delta = (1, 1, 1, 1, 1, 1)}.}
      }

      Only applicable to crisp endorelations.
    }
  }
    
  Methods \code{"symdiff"}, \code{"manhattan"}, \code{"euclidean"} and
  \code{"Jaccard"} take an additional logical argument \code{na.rm}: if
  true (default: false), tuples with missing memberships are excluded in
  the dissimilarity computations.
}
\references{
  W. D. Cook, M. Kress and L. M. Seiford (1986),
  Information and preference in partial orders: a bimatrix
  representation.
  \emph{Psychometrika} \bold{51}/2, 197--207.
  \doi{10.1007/BF02293980}.

  W. D. Cook and L. M. Seiford (1978),
  Priority ranking and consensus formation.
  \emph{Management Science}, \bold{24}/16, 1721--1732.
  \doi{10.1287/mnsc.24.16.1721}.

  P. Diaconis (1988),
  \emph{Group Representations in Probability and Statistics}.
  Institute of Mathematical Statistics: Hayward, CA.
  
  E. J. Emond and D. W. Mason (2000),
  \emph{A new technique for high level decision support}.
  Technical Report ORD Project Report PR2000/13, Operational
  Research Division, Department of National Defence, Canada.
  \url{https://cradpdf.drdc-rddc.gc.ca/PDFS/zbc84/p513851.pdf}.

  K. Jabeur, J.-M. Martel and S. Ben Khélifa (2004).
  A distance-based collective preorder integrating the relative
  importance of the groups members.
  \emph{Group Decision and Negotiation}, \bold{13}, 327--349.
  \doi{10.1023/B:GRUP.0000042894.00775.75}.

  J. G. Kemeny and J. L. Snell (1962),
  \emph{Mathematical Models in the Social Sciences},
  chapter \dQuote{Preference Rankings: An Axiomatic Approach}.
  MIT Press: Cambridge.
}
\keyword{math}
