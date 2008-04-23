\name{consensus}
\encoding{UTF-8}
\alias{relation_consensus}
\title{Consensus Relations}
\description{Compute the consensus relation of a relation ensemble.}
\usage{
relation_consensus(x, method = NULL, weights = 1, control = list(), ...)
}
\arguments{
  \item{x}{an ensemble of relations, or something which can be coerced
    to such (see \code{\link{relation_ensemble}}).}
  \item{method}{a character string specifying one of the built-in
    methods for computing consensus relations, or a function to be
    taken as a user-defined method, or \code{NULL} (default value).  If
    a character string, its lower-cased version is matched against the
    lower-cased names of the available built-in methods using
    \code{\link{pmatch}}.  See \bold{Details} for available built-in
    methods and defaults.}
  \item{weights}{a numeric vector with non-negative case weights.
    Recycled to the number of elements in the ensemble given by \code{x}
    if necessary.}
  \item{control}{a list of control parameters.  See \bold{Details}.}
  \item{\dots}{a list of control parameters (overruling those specified
    in \code{control}).}
}
\value{
  The consensus relation(s).
}
\details{
  Consensus relations \dQuote{synthesize} the information in the
  elements of a relation ensemble into a single relation, often by
  minimizing a criterion function measuring how dissimilar consensus
  candidates are from the (elements of) the ensemble (the so-called
  \dQuote{optimization approach}), typically of the form
  \eqn{L(R) = \sum w_b d(R_b, R) ^ p}, where \eqn{d} is a suitable
  dissimilarity measure (see \code{\link{relation_dissimilarity}}),
  \eqn{w_b} is the case weight given to element \eqn{R_b} of the
  ensemble, and \eqn{p \ge 1}.  Such consensus relations are called
  \dQuote{central relations} in Régnier (1965).  For \eqn{p = 1}, we
  obtain (generalized) medians; \eqn{p = 2} gives (generalized) means
  (least squares consensus relations).

  Available built-in methods are as follows.  Apart from Condorcet's and
  the unrestricted Manhattan and Euclidean consensus methods, these are
  applicable to ensembles of endorelations only. 

  \describe{
    \item{\code{"Borda"}}{the consensus method proposed by Borda (1781).
      For each relation \eqn{R_b} and object \eqn{x}, one determines the
      Borda/Kendall scores, i.e., the number of objects \eqn{y} such
      that \eqn{y R_b x}.  These are then aggregated across relations
      by weighted averaging.  Finally, objects are ordered according to
      their aggregated scores.}
    \item{\code{"Copeland"}}{the consensus method proposed by Copeland
      (1951).  For each relation \eqn{R_b} and object \eqn{x}, one
      determines the Copeland scores, i.e., the number of objects
      \eqn{y} such that \eqn{y R_b x}, minus the number of objects
      \eqn{y} such that \eqn{x R_b y}.  Like the Borda method, these are
      then aggregated across relations by weighted averaging.  Finally,
      objects are ordered according to their aggregated scores.}
    \item{\code{"Condorcet"}}{the consensus method proposed by Condorcet
      (1785).  For a given ensemble of crisp relations, this minimizes
      the criterion function \eqn{L} with \eqn{d} as symmetric
      difference distance and \eqn{p = 1} over all possible crisp 
      relations.  In the case of endorelations, consensus is obtained by
      weighting voting, such that \eqn{x R y} if the weighted number of
      times that \eqn{x R_b y} is no less than the weighted number of
      times that this is not the case.  Even when aggregating linear
      orders, this can lead to intransitive consensus solutions
      (\dQuote{effet Condorcet}).}
    \item{\code{"CS"}}{the consensus method of Cook and Seiford (1978)
      which determines a linear order minimizing the criterion function
      \eqn{L} with \eqn{d} as generalized Cook-Seiford (ranking)
      distance and \eqn{p = 1} via solving a linear sum assignment
      problem.  One can obtain a relation ensemble with \emph{all}
      consensus relations by setting the control parameter \code{all} to
      \code{TRUE}.}
    \item{\code{"SD/\var{F}"}}{an exact solver for determining the
      consensus relation of an ensemble of crisp endorelations by
      minimizing the criterion function \eqn{L} with \eqn{d} as
      symmetric difference distance (\dQuote{SD}) and \eqn{p = 1} over a
      suitable class (\dQuote{Family}) of crisp endorelations as
      indicated by \var{F}, with values:
      \describe{
	\item{\code{E}}{equivalence relations: reflexive, symmetric, and
	  transitive.}
	\item{\code{L}}{linear orders: complete, antisymmetric, and
	  transitive.}
	\item{\code{O}}{partial orders: reflexive, antisymmetric and
	  transitive.}
	\item{\code{P}}{complete preorders (preference relations,
	  \dQuote{orderings}): complete, reflexive and transitive.}
	\item{\code{T}}{tournaments: complete, irreflexive and
	  antisymmetric (i.e., complete and asymmetric).}
	\item{\code{C}}{complete relations.}
	\item{\code{A}}{antisymmetric relations.}
	\item{\code{S}}{symmetric relations.}
	\item{\code{M}}{matches: complete and reflexive.}
	\item{\code{R}}{tRansitive relations.}
      }
      Consensus relations are determined by reformulating the consensus
      problem as an integer program (for the relation incidences), which
      is solved via package \pkg{lpSolve}.  See Hornik and Meyer (2007)
      for details.

      For fitting equivalences and preferences (cases \code{E} and
      \code{P}) it is possible to specify a maximal number of classes
      \eqn{k} using the control parameter \code{k}.

      Additional constraints on the incidences of the consensus solution
      can be given via the control parameter \code{constraints}, in the
      form of a 3-column matrix whose rows give row and column indices
      \eqn{i} and \eqn{j} and the corresponding incidence \eqn{I_{ij}}.
      (I.e., incidences can be constrained to be zero or one on an
      object by object basis.)

      One can obtain a relation ensemble with \emph{all} consensus
      relations by setting the control parameter \code{all} to
      \code{TRUE}.
      (See the examples.)
    }
    \item{\code{"manhattan"}}{the (unrestricted) median of the
      ensemble, minimizing \eqn{L} with \eqn{d} as Manhattan (symmetric
      difference) distance and \eqn{p = 1} over all (possibly fuzzy)
      relations.}
    \item{\code{"euclidean"}}{the (unrestricted) mean of the ensemble,
      minimizing \eqn{L} with \eqn{d} as Euclidean distance and
      \eqn{p = 2} over all (possibly fuzzy) relations.}
    \item{\code{"euclidean/\var{F}"}}{an exact solver for determining
      the restricted least squares Euclidean consensus relation of an
      ensemble of endorelations by minimizing the criterion function
      \eqn{L} with \eqn{d} as Euclidean difference distance and
      \eqn{p = 2} over a suitable family of crisp endorelations as
      indicated by \var{F}, with available families control parameters
      as for methods \code{"SD/\var{F}"}.
    }
  }
}
\references{
  J. C. Borda (1781),
  Mémoire sur les élections au scrutin.
  Histoire de l'Académie Royale des Sciences.
  
  W. D. Cook and M. Kress (1992),
  \emph{Ordinal information and preference structures: decision models and
    applications}.
  Prentice-Hall: New York.
  ISBN: 0-13-630120-7.

  W. D. Cook and L. M. Seiford (1978),
  Priority ranking and consensus formation.
  \emph{Management Science}, \bold{24}/16, 1721--1732.
  
  M. J. A. de Condorcet (1785),
  Essai sur l'application de l'analyse à la probabilité des décisions
  rendues à la pluralité des voix. 
  Paris.

  A. H. Copeland (1951),
  A Reasonable Social Welfare Function.
  \emph{mimeo}, University of Michigan.
  
  K. Hornik and D. Meyer (2007),
  Deriving consensus rankings from benchmarking experiments.
  In R. Decker and H.-J. Lenz,
  \emph{Advances in Data Analysis}.
  Studies in Classification, Data Analysis, and Knowledge Organization.
  Springer-Verlag: Heidelberg, 163--170.
  
  F. Marcotorchino and P. Michaud (1982).
  Agrégation de similarités en classification automatique.
  \emph{Revue de Statistique Appliquée}, \bold{30}(2):21--44.
  \url{http://www.numdam.org/item?id=RSA_1982__30_2_21_0}.
  
  S. Régnier (1965),
  Sur quelques aspects mathématiques des problèmes de classification
  automatique.
  \emph{ICC Bulletin}, \bold{4}:175--191.
}
\examples{
## Consensus equivalence.
## (I.e., in fact, consensus partition.)
## Classification of 30 felines, see Marcotorchino and Michaud (1982).
data("Felines")
## Consider each variable an equivalence relation on the objects.
relations <- as.relation_ensemble(Felines)
## This gives a relation ensemble of length 14 (number of variables in
## the data set).
## Now fit an equivalence relation to this:
E <- relation_consensus(relations, "SD/E")
## And look at the equivalence classes:
ids <- relation_class_ids(E)
## Or, more nicely:
split(rownames(Felines), ids)
## Which is the same as in the paper ...

## Consensus linear order.
## Example from Cook and Kress, pages 48ff.
## Relation from paired comparisons.
pm <- matrix(c(0, 1, 0, 1, 1,
               0, 0, 0, 1, 1,
               1, 1, 0, 0, 0,
               0, 0, 1, 0, 0,
               0, 0, 1, 1, 0),
             nrow = 5,
             byrow = TRUE,
             dimnames = list(letters[1:5], letters[1:5]))
## Note that this is a Cook and Kress "preference" matrix.
R <- as.relation(1 - pm)
relation_incidence(R)
relation_is_tournament(R)
## Closest linear order:
L <- relation_consensus(R, "SD/L")
relation_incidence(L)
## Visualize provided that Rgraphviz is available.
if(require("Rgraphviz")) plot(L)
## But note that this linear order is not unique.
L <- relation_consensus(R, "SD/L", control = list(all = TRUE))
print(L)
if(require("Rgraphviz")) plot(L)
## (Oh no.)
## Closest preference relation with at most 3 indifference classes:
P3 <- relation_consensus(R, "SD/P", control = list(k = 3))
relation_incidence(P3)
}
\keyword{math}
