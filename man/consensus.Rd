\name{consensus}
\encoding{UTF-8}
\alias{relation_consensus}
\title{Consensus Relations}
\description{
  Compute consensus relations of a relation ensemble.
}
\usage{
relation_consensus(x, method = NULL, weights = 1,
                   control = list(), ...)
}
\arguments{
  \item{x}{an ensemble of relations (see
    \code{\link{relation_ensemble}()}),
    or something which can be coerced to such.}
  \item{method}{a character string specifying one of the built-in
    methods for computing consensus relations, or a function to be
    taken as a user-defined method, or \code{NULL} (default value).  If
    a character string, its lower-cased version is matched against the
    lower-cased names of the available built-in methods using
    \code{\link{pmatch}()}.  See \bold{Details} for available built-in
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
  \eqn{\Phi(R) = \sum w_b d(R_b, R) ^ e}, where \eqn{d} is a suitable
  dissimilarity measure (see \code{\link{relation_dissimilarity}()}),
  \eqn{w_b} is the case weight given to element \eqn{R_b} of the
  ensemble, and \eqn{e \ge 1}.  Such consensus relations are called
  \dQuote{central relations} in Régnier (1965).  For \eqn{e = 1}, we
  obtain (generalized) medians; \eqn{e = 2} gives (generalized) means
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
      their aggregated scores. Note that this may result in a weak order
      (i.e., with objects being tied).

      One can enforce a linear order by setting the control parameter
      \code{L} to \code{TRUE}, and obtain a relation ensemble with up to
      \emph{n} or all such solutions by additionally setting the control
      parameter \code{n} to some positive integer or \code{"all"},
      respectively.
    }
    \item{\code{"Copeland"}}{the consensus method proposed by Copeland
      (1951).  For each relation \eqn{R_b} and object \eqn{x}, one
      determines the Copeland scores, i.e., the number of objects
      \eqn{y} such that \eqn{y R_b x}, minus the number of objects
      \eqn{y} such that \eqn{x R_b y}.  Like the Borda method, these are
      then aggregated across relations by weighted averaging.  Finally,
      objects are ordered according to their aggregated scores.
      Note that this may result in a weak order
      (i.e., with objects being tied).

      One can enforce a linear order by setting the control parameter
      \code{L} to \code{TRUE}, and obtain a relation ensemble with up to
      \emph{n} or all such solutions by additionally setting the control
      parameter \code{n} to some positive integer or \code{"all"},
      respectively.
    }
    \item{\code{"Condorcet"}}{the consensus method proposed by Condorcet
      (1785).  For a given ensemble of crisp relations, this minimizes
      the criterion function \eqn{\Phi} with \eqn{d} as symmetric
      difference distance and \eqn{e = 1} over all possible crisp 
      relations.  In the case of endorelations, consensus is obtained by
      weighting voting, such that \eqn{x R y} if the weighted number of
      times that \eqn{x R_b y} is no less than the weighted number of
      times that this is not the case.  Even when aggregating linear
      orders, this can lead to intransitive consensus solutions
      (\dQuote{effet Condorcet}).

      One can obtain a relation ensemble with up to \emph{n} or all such
      solutions consensus relations by setting the control parameter
      \code{n} to some positive integer or \code{"all"}, respectively.
    } 
    \item{\code{"CS"}}{the consensus method of Cook and Seiford (1978)
      which determines a linear order minimizing the criterion function
      \eqn{\Phi} with \eqn{d} as generalized Cook-Seiford (ranking)
      distance and \eqn{e = 1} via solving a linear sum assignment
      problem.

      One can obtain a relation ensemble with up to \emph{n} or all such
      consensus relations by setting the control parameter \code{n} to
      some positive integer or \code{"all"}, respectively.
    }
    \item{\code{"symdiff/\var{F}"}}{an exact solver for determining the
      consensus relation of an ensemble of crisp endorelations by
      minimizing the criterion function \eqn{\Phi} with \eqn{d} as
      symmetric difference (\dQuote{symdiff}) distance and \eqn{e = 1}
      over a suitable class (\dQuote{Family}) of crisp endorelations as
      indicated by \var{F}, with values:
      \describe{
	\item{\code{G}}{general (crisp) endorelations.}
	\item{\code{A}}{antisymmetric relations.}
	\item{\code{C}}{complete relations.}
	\item{\code{E}}{equivalence relations: reflexive, symmetric, and
	  transitive.}
	\item{\code{L}}{linear orders: complete, reflexive,
	  antisymmetric, and transitive.}
	\item{\code{M}}{matches: complete and reflexive.}
	\item{\code{O}}{partial orders: reflexive, antisymmetric and
	  transitive.}
	\item{\code{S}}{symmetric relations.}	
	\item{\code{T}}{tournaments: complete, irreflexive and
	  antisymmetric (i.e., complete and asymmetric).}
	\item{\code{W}}{weak orders (complete preorders, preferences,
	  \dQuote{orderings}): complete, reflexive and transitive.} 
	\item{\code{preorder}}{preorders: reflexive and transitive.}
	\item{\code{transitive}}{transitive relations.}
      }

      Can also be referred to as \code{"SD/\var{F}"}.
      
      Consensus relations are determined by reformulating the consensus
      problem as a binary program (for the relation incidences), see
      Hornik and Meyer (2007) for details.  The solver employed can be
      specified via the control argument \code{solver}, with currently
      possible values \code{"glpk"}, \code{"lpsolve"}, \code{"symphony"}
      or \code{"cplex"} or a unique abbreviation thereof, specifying to
      use the solvers from packages \pkg{Rglpk} (default),
      \pkg{lpSolve}, \pkg{Rsymphony}, or \pkg{Rcplex}, respectively.
      Unless control option \code{sparse} is false, a sparse formulation
      of the binary program is used, which is typically more efficient.

      For fitting equivalences and weak orders (cases \code{E} and
      \code{W}) it is possible to specify the number of classes \eqn{k}
      using the control parameter \code{k}.  For fitting weak orders,
      one can also specify the number of elements in the classes via
      control parameter \code{l}.

      Additional constraints on the incidences of the consensus solution
      can be given via the control parameter \code{constraints}, in the
      form of a 3-column matrix whose rows give row and column indices
      \eqn{i} and \eqn{j} and the corresponding incidence \eqn{I_{ij}}.
      (I.e., incidences can be constrained to be zero or one on an
      object by object basis.)

      One can obtain a relation ensemble with up to \emph{n} or all such
      consensus relations by setting the control parameter \code{n} to
      some positive integer or \code{"all"}, respectively.
      (See the examples.)
    }
    \item{\code{"manhattan"}}{the (unrestricted) median of the
      ensemble, minimizing \eqn{\Phi} with \eqn{d} as Manhattan (symmetric
      difference) distance and \eqn{e = 1} over all (possibly fuzzy)
      relations.}
    \item{\code{"euclidean"}}{the (unrestricted) mean of the ensemble,
      minimizing \eqn{\Phi} with \eqn{d} as Euclidean distance and
      \eqn{e = 2} over all (possibly fuzzy) relations.}
    \item{\code{"euclidean/\var{F}"}}{an exact solver for determining
      the restricted least squares Euclidean consensus relation of an
      ensemble of endorelations by minimizing the criterion function
      \eqn{\Phi} with \eqn{d} as Euclidean difference distance and
      \eqn{e = 2} over a suitable family of crisp endorelations as
      indicated by \var{F}, with available families and control
      parameters as for methods \code{"symdiff/\var{F}"}.
    }
    \item{\code{"majority"}}{a generalized majority method for which the
      consensus relation contains of all tuples occurring with a
      relative frequency of more than \eqn{100 p} percent (of 100
      percent if \eqn{p = 1}).  The fraction \eqn{p} can be specified
      via the control parameter \code{p}.  By default, \eqn{p = 1/2} is
      used.
    }
    \item{\code{"CKS/\var{F}"}}{an exact solver for determining the
      consensus relation of an ensemble of crisp endorelations by
      minimizing the criterion function \eqn{\Phi} with \eqn{d} as
      Cook-Kress-Seiford (\dQuote{CKS}) distance and \eqn{e = 1} over a
      suitable class (\dQuote{Family}) of crisp endorelations as
      indicated by \var{F}, with available families and control
      parameters as for methods \code{"symdiff/\var{F}"}.

      For fitting equivalences and weak orders (cases \code{E} and
      \code{W}) it is possible to specify the number of classes \eqn{k}
      using the control parameter \code{k}.
      
      One can obtain a relation ensemble with up to \emph{n} or all such
      consensus relations by setting the control parameter \code{n} to
      some positive integer or \code{"all"}, respectively.
    }
    \item{\code{"PC/\var{F}"}}{an exact solver for determining the
      consensus relation of an ensemble of crisp endorelations by
      minimizing the criterion function \eqn{\Phi} with \eqn{d} as
      (generalized) paired comparison (\dQuote{PC}) distance and
      \eqn{e = 1} over a suitable class (\dQuote{Family}) of crisp
      endorelations as indicated by \var{F}, with available families and
      control parameters as for methods \code{"symdiff/\var{F}"}, and
      control option \code{delta} for specifying the paired comparison
      discrepancies.

      For fitting equivalences and weak orders (cases \code{E} and
      \code{W}) it is possible to specify the number of classes \eqn{k}
      using the control parameter \code{k}.
      
      One can obtain a relation ensemble with up to \emph{n} or all such
      consensus relations by setting the control parameter \code{n} to
      some positive integer or \code{"all"}, respectively.
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
  \doi{10.1287/mnsc.24.16.1721}.
  
  M. J. A. de Condorcet (1785),
  Essai sur l'application de l'analyse à la probabilité des décisions
  rendues à la pluralité des voix. 
  Paris.

  A. H. Copeland (1951),
  A Reasonable Social Welfare Function.
  \emph{mimeo}, University of Michigan.
  
  E. J. Emond and D. W. Mason (2000),
  \emph{A new technique for high level decision support}.
  Technical Report ORD Project Report PR2000/13, Operational
  Research Division, Department of National Defence, Canada.
  %% \url{https://cradpdf.drdc-rddc.gc.ca/PDFS/zbc84/p513851.pdf}.

  K. Hornik and D. Meyer (2007),
  Deriving consensus rankings from benchmarking experiments.
  In R. Decker and H.-J. Lenz,
  \emph{Advances in Data Analysis}.
  Studies in Classification, Data Analysis, and Knowledge Organization.
  Springer-Verlag: Heidelberg, 163--170.
  
  F. Marcotorchino and P. Michaud (1982).
  Agrégation de similarités en classification automatique.
  \emph{Revue de Statistique Appliquée}, \bold{30}/2, 21--44.
  \url{https://eudml.org/doc/106132}.
  
  S. Régnier (1965),
  Sur quelques aspects mathématiques des problèmes de classification
  automatique.
  \emph{ICC Bulletin}, \bold{4}, 175--191.
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
E <- relation_consensus(relations, "symdiff/E")
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
## Note that this is a Cook and Kress "preference matrix" where entry
## (i,j) is one iff object i is preferred to object j (i > j).
## Set up the corresponding '<' relation:
R <- as.relation(t(pm))
relation_incidence(R)
relation_is(R, "tournament")
## Closest linear order:
L <- relation_consensus(R, "symdiff/L")
relation_incidence(L)
## Visualize provided that Rgraphviz is available.
if(require("Rgraphviz")) plot(L)
## But note that this linear order is not unique.
L <- relation_consensus(R, "symdiff/L", control = list(n = "all"))
print(L)
if(require("Rgraphviz")) plot(L)
## (Oh no: c is once first and once last.)
## Closest weak order relation with at most 3 indifference classes:
W3 <- relation_consensus(R, "symdiff/W", control = list(k = 3))
relation_incidence(W3)

## Consensus weak orders.
## Example from Emond and Mason, pages 28f.
## The reference provides 21 partial rankings of 15 objects,
## in 3 groups of 7 rankings (corresponding to three different
## ranking criteria) with respective weights 4, 5, and 7.
wei <- rep.int(c(4, 5, 7), rep(7, 3))
## The rankings are written by listing the object labels from the
## best to the worst, with a leading minus indicating a tie with
## the previous object:
EM_inputs <-
    c("6 1 -7 -9 10 3 8 11 5 -12 2 -4 -13",
      "6 10 9 3 4 -8 7 1 -5 -11 2 12 13 14 15",
      "6 10 3 7 8 11 5 14 15 12 1 -4 -13 2 -9",
      "6 9 -11 10 3 14 12 7 4 5 2 1 8 13 15",
      "10 6 7 1 11 -13 4 2 3 9 12 14 -15 8 5",
      "6 9 8 -10 11 4 1 5 7 15 2 12 14 13 3",
      "1 -6 -10 7 -12 9 3 4 -11 -14 -15 2 -13 8",
      "4 -10 1 -7 6 -9 -13 5 -14 3 12 8 11 -15 2",
      "4 -9 5 1 14 11 8 3 6 2 -13 10 12 7 15",
      "4 2 -5 8 15 7 11 -14 1 -12 -13 10 9 6",
      "2 -11 -12 -14 -15 6 -13 3 -4 9 8 -10 1 -5 -7",
      "4 14 10 2 5 3 1 13 12 7 15 8 11 6 9",
      "4 2 5 1 15 7 13 14 3 -12 8 11 6 9 10",
      "12 1 3 -4 2 11 -13 -15 9 14 6 8 7 -10 5",
      "5 4 9 2 -7 14 8 -11 3 1 15 12 6 10 13",
      "11 9 -14 15 12 3 4 13 8 6 7 10 5",
      "12 11 2 1 3 9 8 10 13 -14 6 4 -15 5 7",
      "4 -5 10 -12 3 8 -11 6 -7 -9 13 14 15",
      "12 5 -13 14 3 8 15 4 9 -10 11 6 7",
      "4 -5 -8 11 6 14 7 1 -2 -15 10 3 13 9 -12",
      "10 8 5 -11 6 -14 9 4 -13 -15 3 -12 2 1")
## Using the Emond-Mason paired comparison dissimilarity, there
## are three consensus rankings when using the above weights:
EM_solutions <-
    c("4 10 5-11 1 -2-14 3-12 9 8 6 7 13-15",
      "4 10 5-11 1 -2 9 14 3-12 8 6 7 13-15",
      "4 10 5-11 2-14 1 3-12 9 8 6 7 13-15")
## We can reproduce this as follows.
## We first provide a reader for the rankings, and a maker for
## creating the (possibly partial) ranking with the appropriate
## domain:
reader <- function(s) {
    strsplit(unlist(strsplit(gsub(" *-", "-", s),
                             " +")),
             "-",
             fixed = TRUE)
}
maker <- function(s) {
    ranking(lapply(reader(s), as.numeric),
            domain = as.numeric(1 : 15))
}
EM_inputs <- lapply(EM_inputs, maker)
EM_solutions <- lapply(EM_solutions, maker)
## Package 'relations' uses NA for non-diagonal incidences
## featuring unranked objects.
## Following the reference, we impute these by zeroes:
ens <- relation_impute(relation_ensemble(list = EM_inputs), "omit")
## We can now obtain all consensus weak orders (corresponding to
## complete rankings) as follows:
con <- relation_consensus(ens, "PC/W", wei, delta = "EM", all = TRUE)
## To verify that these agree with the solutions given in the
## reference:
sets::set_outer(con, relation_ensemble(list = EM_solutions), `==`)
}
\keyword{math}
