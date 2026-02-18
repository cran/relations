\name{predicates}
\alias{predicates}
\alias{relation_is}
\alias{relation_is_Euclidean}
\alias{relation_is_Ferrers}
\alias{relation_is_acyclic}
\alias{relation_is_antisymmetric}
\alias{relation_is_asymmetric}
\alias{relation_is_bijective}
\alias{relation_is_binary}
\alias{relation_is_complete}
\alias{relation_is_coreflexive}
\alias{relation_is_crisp}
\alias{relation_is_cyclic}
\alias{relation_is_endorelation}
\alias{relation_is_equivalence}
\alias{relation_is_functional}
\alias{relation_is_homogeneous}
\alias{relation_is_injective}
\alias{relation_is_interval_order}
\alias{relation_is_irreflexive}
\alias{relation_is_left_total}
\alias{relation_is_linear_order}
\alias{relation_is_match}
\alias{relation_is_negatively_transitive}
\alias{relation_is_partial_order}
\alias{relation_is_preference}
\alias{relation_is_preorder}
\alias{relation_is_quasiorder}
\alias{relation_is_quasitransitive}
\alias{relation_is_quaternary}
\alias{relation_is_reflexive}
\alias{relation_is_right_total}
\alias{relation_is_semiorder}
\alias{relation_is_semitransitive}
\alias{relation_is_strict_linear_order}
\alias{relation_is_strict_partial_order}
\alias{relation_is_strongly_complete}
\alias{relation_is_surjective}
\alias{relation_is_symmetric}
\alias{relation_is_ternary}
\alias{relation_is_tournament}
\alias{relation_is_transitive}
\alias{relation_is_trichotomous}
\alias{relation_is_weak_order}
\alias{relation_has_missings}
\title{Relation Predicates}
\description{
  Predicate functions for testing for binary relations and
  endorelations, and special kinds thereof.
}
\usage{
relation_is(x, predicate, ...)
relation_is_Euclidean(x, na.rm = FALSE)
relation_is_Ferrers(x, na.rm = FALSE)
relation_is_acyclic(x)
relation_is_antisymmetric(x, na.rm = FALSE)
relation_is_asymmetric(x, na.rm = FALSE)
relation_is_bijective(x)
relation_is_binary(x)
relation_is_complete(x, na.rm = FALSE)
relation_is_coreflexive(x, na.rm = FALSE)
relation_is_crisp(x, na.rm = FALSE)
relation_is_cyclic(x)
relation_is_endorelation(x)
relation_is_equivalence(x, na.rm = FALSE)
relation_is_functional(x)
relation_is_homogeneous(x)
relation_is_injective(x)
relation_is_interval_order(x, na.rm = FALSE)
relation_is_irreflexive(x, na.rm = FALSE)
relation_is_left_total(x)
relation_is_linear_order(x, na.rm = FALSE)
relation_is_match(x, na.rm = FALSE)
relation_is_negatively_transitive(x, na.rm = FALSE)
relation_is_partial_order(x, na.rm = FALSE)
relation_is_preference(x, na.rm = FALSE)
relation_is_preorder(x, na.rm = FALSE)
relation_is_quasiorder(x, na.rm = FALSE)
relation_is_quasitransitive(x, na.rm = FALSE)
relation_is_quaternary(x)
relation_is_reflexive(x, na.rm = FALSE)
relation_is_right_total(x)
relation_is_semiorder(x, na.rm = FALSE)
relation_is_semitransitive(x, na.rm = FALSE)
relation_is_strict_linear_order(x, na.rm = FALSE)
relation_is_strict_partial_order(x, na.rm = FALSE)
relation_is_strongly_complete(x, na.rm = FALSE)
relation_is_surjective(x)
relation_is_symmetric(x, na.rm = FALSE)
relation_is_ternary(x)
relation_is_tournament(x, na.rm = FALSE)
relation_is_transitive(x, na.rm = FALSE)
relation_is_trichotomous(x, na.rm = FALSE)
relation_is_weak_order(x, na.rm = FALSE)
relation_has_missings(x)
}
\arguments{
  \item{x}{an object inheriting from class \code{\link{relation}}.}
  \item{na.rm}{a logical indicating whether tuples with missing
    memberships are excluded in the predicate computations.}
  \item{predicate}{character vector matching one of the following (see \bold{details}):
\code{"binary"},
\code{"ternary"},
\code{"quaternary"},
\code{"left_total"},
\code{"right_total"},
\code{"surjective"},
\code{"functional"},
\code{"injective"},
\code{"bijective"},
\code{"endorelation"},
\code{"homogeneous"},
\code{"crisp"},
\code{"complete"},
\code{"match"},
\code{"strongly_complete"},
\code{"reflexive"},
\code{"irreflexive"},
\code{"coreflexive"},
\code{"symmetric"},
\code{"asymmetric"},
\code{"antisymmetric"},
\code{"transitive"},
\code{"negatively_transitive"},
\code{"quasitransitive"},
\code{"Ferrers"},
\code{"semitransitive"},
\code{"trichotomous"},
\code{"Euclidean"},
\code{"equivalence"},
\code{"weak_order"},
\code{"preference"},
\code{"preorder"},
\code{"quasiorder"},
\code{"partial_order"},
\code{"linear_order"},
\code{"strict_partial_order"},
\code{"strict_linear_order"},
\code{"tournament"},
\code{"interval_order"},
\code{"semiorder"},
\code{"acyclic"}
\code{"cyclic"}
  }
\item{\dots}{Additional arguments passed to the predicate functions
  (currently, only \code{na.rm} for most predicates).}
}
\details{
  This help page documents the predicates currently available. Note that
  the preferred way is to use the meta-predicate function
  \code{relation_is(x, "FOO")} instead of the individual predicates
  \code{relation_is_FOO(x)} since the latter will become deprecated in
  future releases.
  
  A binary relation is a relation with arity 2.

  A relation \eqn{R} on a set \eqn{X} is called
  \emph{homogeneous} iff \eqn{D(R) = (X, \dots, X)}.

  An \emph{endorelation} is a binary homogeneous relation.

  For a crisp binary relation, let us write \eqn{x R y} iff \eqn{(x, y)}
  is contained in \eqn{R}.

  A crisp binary relation \eqn{R} is called
  \describe{
    \item{left-total:}{for all \eqn{x} there is at least one \eqn{y}
      such that \eqn{x R y}.}
    \item{right-total:}{for all \eqn{y} there is at least one \eqn{x}
      such that \eqn{x R y}.}
    \item{functional:}{for all \eqn{x} there is at most one \eqn{y}
      such that \eqn{x R y}.}
    \item{surjective:}{the same as right-total.}
    \item{injective:}{for all \eqn{y} there is at most one \eqn{x}
      such that \eqn{x R y}.}
    \item{bijective:}{left-total, right-total, functional and
      injective.}
  }

  A crisp endorelation \eqn{R} is called
  \describe{
    \item{reflexive:}{\eqn{x R x} for all \eqn{x}.}
    \item{irreflexive:}{there is no \eqn{x} such that \eqn{x R x}.}
    \item{coreflexive:}{\eqn{x R y} implies \eqn{x = y}.}
    \item{symmetric:}{\eqn{x R y} implies \eqn{y R x}.}
    \item{asymmetric:}{\eqn{x R y} implies that not \eqn{y R x}.}
    \item{antisymmetric:}{\eqn{x R y} and \eqn{y R x} imply that
      \eqn{x = y}.}
    \item{transitive:}{\eqn{x R y} and \eqn{y R z} imply that
      \eqn{x R z}.}
    \item{complete:}{for all distinct \eqn{x} and \eqn{y}, \eqn{x R y}
      or \eqn{y R x}.}
    \item{strongly complete:}{for all \eqn{x} and \eqn{y}, \eqn{x R y}
      or \eqn{y R x} (i.e., complete and reflexive).}
    \item{negatively transitive:}{not \eqn{x R y} and not \eqn{y R z}
      imply that not \eqn{x R z}.}
    \item{Ferrers:}{\eqn{x R y} and \eqn{z R w} imply \eqn{x R w} or
      \eqn{y R z}.}
    \item{semitransitive:}{\eqn{x R y} and \eqn{y R z} imply \eqn{x R w}
      or \eqn{w R z}.}
    \item{quasitransitive:}{
      \eqn{x R y} and not \eqn{y R x} and
      \eqn{y R z} and not \eqn{z R y} imply
      \eqn{x R z} and not \eqn{z R x}
      (i.e., the asymmetric part of \eqn{R} is transitive).}
    \item{trichotomous:}{exactly one of \eqn{x R y}, \eqn{y R x}, or
      \eqn{x = y} holds.}
    \item{Euclidean:}{\eqn{x R y} and \eqn{x R z} imply \eqn{y R z}.}
    \item{acyclic:}{the transitive closure of R is antisymmetric.}
    \item{cyclic:}{R is not acyclic.}
  }

  Some combinations of these basic properties have special names
  because of their widespread use:
  \describe{
    \item{preorder:}{reflexive and transitive.}
    \item{quasiorder:}{the same as preorder.}
    \item{equivalence:}{a symmetric preorder (reflexive, symmetric, and
      transitive).}
    \item{weak order:}{a complete preorder (complete, reflexive, and
      transitive).}
    \item{preference:}{the same as weak order.}
    \item{partial order:}{an antisymmetric preorder (reflexive,
      antisymmetric, and transitive).}
    \item{strict partial order:}{irreflexive, antisymmetric, and
      transitive, or equivalently: asymmetric and transitive).}
    \item{linear order:}{a complete partial order.}
    \item{strict linear order:}{a complete strict partial order.}
    \item{match:}{strongly complete.}
    \item{tournament:}{complete and asymmetric.}
    \item{interval order:}{complete and Ferrers.}
    \item{semiorder:}{a semitransitive interval order.}
  }

  If \eqn{R} is a weak order (\dQuote{(weak) preference relation}),
  \eqn{I = I(R)} defined by \eqn{x I y} iff \eqn{x R y} and \eqn{y R x}
  is an equivalence, the \emph{indifference relation} corresponding to
  \eqn{R}.

  There seem to be no commonly agreed definitions for order relations:
  e.g., \bibcitet{Fishburn:1972} requires these to be irreflexive.

  For a fuzzy binary relation \eqn{R}, let \eqn{R(x, y)} denote the
  membership of \eqn{(x, y)} in the relation.  Write \eqn{T} and \eqn{S}
  for the fuzzy t-norm (intersection) and t-conorm (disjunction),
  respectively (min and max for the \dQuote{standard} Zadeh family).
  Then generalizations of the above basic endorelation predicates are as
  follows.
  \describe{
    \item{reflexive:}{\eqn{R(x, x) = 1} for all \eqn{x}.}
    \item{irreflexive:}{\eqn{R(x, x) = 0} for all \eqn{x}.}
    \item{coreflexive:}{\eqn{R(x, y) > 0} implies \eqn{x = y}.}
    \item{symmetric:}{
      \eqn{R(x, y) = R(y, x)} for all \eqn{x \ne y}{x != y}.}
    \item{asymmetric:}{
      \eqn{T(R(x, y), R(y, x)) = 0} for all \eqn{x, y}.} 
    \item{antisymmetric:}{
      \eqn{T(R(x, y), R(y, x)) = 0} for all \eqn{x \ne y}{x != y}.}
    \item{transitive:}{\eqn{T(R(x, y), R(y, z)) \le R(x, z)} for all
      \eqn{x, y, z}.}
    \item{complete:}{
      \eqn{S(R(x, y), R(y, x)) = 1} for all \eqn{x \ne y}{x != y}.}
    \item{strongly complete:}{
      \eqn{S(R(x, y), R(y, x)) = 1} for all \eqn{x, y}.}
    \item{negatively transitive:}{
      \eqn{R(x, z) \le S(R(x, y), R(y, z))} for all \eqn{x, y, z}.}
    \item{Ferrers:}{
      \eqn{T(R(x, y), R(z, w)) \le S(R(x, w), R(z, y))} for all
      \eqn{x, y, z, w}.}
    \item{semitransitive:}{
      \eqn{T(R(x, w), R(w, y)) \le S(R(x, z), R(z, y))} for all
      \eqn{x, y, z, w}.}
  }
  The combined predicates are obtained by combining the basic predicates
  as for crisp endorelations (see above).

  A relation has missings iff at least one cell in the incidence matrix
  is \code{NA}. In addition to \code{relation_has_missings()}, an
  \code{is.na} method for relations is available, returning a matrix of
  logicals corresponding to the incidences tested for missingness.
}
\references{
  \bibshow{Fishburn:1972, Varian:2002}
}
\examples{
require("sets")
R <- relation(domain = c(1, 2, 3), graph = set(c(1, 2), c(2, 3)))
summary(R)

## Note the possible effects of NA-handling:
relation_incidence(R)
relation_is(R, "transitive") ## clearly FALSE

relation_incidence(R)[1, 2] <- NA
relation_incidence(R)
relation_is(R, "transitive") ## clearly NA

## The following gives TRUE, since NA gets replaced with 0:
relation_is(R, "transitive", na.rm = TRUE)
}
\keyword{math}
