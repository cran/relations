\name{predicates}
\alias{relation_is_antisymmetric}
\alias{relation_is_asymmetric}
\alias{relation_is_bijective}
\alias{relation_is_binary}
\alias{relation_is_complete}
\alias{relation_is_coreflexive}
\alias{relation_is_endorelation}
\alias{relation_is_equivalence}
\alias{relation_is_functional}
\alias{relation_is_injective}
\alias{relation_is_irreflexive}
\alias{relation_is_left_total}
\alias{relation_is_linear_order}
\alias{relation_is_partial_order}
\alias{relation_is_reflexive}
\alias{relation_is_right_total}
\alias{relation_is_surjective}
\alias{relation_is_symmetric}
\alias{relation_is_tournament}
\alias{relation_is_transitive}
\alias{relation_is_weak_order}
\alias{relation_is_preference}
\alias{relation_is_preorder}
\alias{relation_is_quasiorder}
\title{Relation Predicates}
\description{
  Predicate functions for testing for binary relations and
  endorelations, and special kinds thereof.
}
\usage{
relation_is_antisymmetric(x)
relation_is_asymmetric(x)
relation_is_bijective(x)
relation_is_binary(x)
relation_is_complete(x)
relation_is_coreflexive(x)
relation_is_endorelation(x)
relation_is_equivalence(x)
relation_is_functional(x)
relation_is_injective(x)
relation_is_irreflexive(x)
relation_is_left_total(x)
relation_is_linear_order(x)
relation_is_partial_order(x)
relation_is_reflexive(x)
relation_is_right_total(x)
relation_is_surjective(x)
relation_is_symmetric(x)
relation_is_tournament(x)
relation_is_transitive(x)
relation_is_weak_order(x)
relation_is_preference(x)
relation_is_preorder(x)
relation_is_quasiorder(x)
}
\arguments{
  \item{x}{an object inheriting from class \code{\link{relation}}.}
}
\details{
  A binary relation is a relation with arity 2.
  An \emph{endorelation} \eqn{R} on a set \eqn{X} is a relation with
  domain \eqn{D(R) = (X, X)}, i.e., a binary relation on \eqn{X}.

  Let us write \eqn{x R y} iff \eqn{(x, y)} is contained in \eqn{R}.

  A binary relation \eqn{R} is called
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

  An endorelation \eqn{R} is called
  \describe{
    \item{reflexive:}{\eqn{x R x} for all \eqn{x}.}
    \item{irreflexive:}{there is no \eqn{x} such that \eqn{x R x}.}
    \item{coreflexive:}{\eqn{x R y} implies \eqn{x = y}.}
    \item{symmetric:}{\eqn{x R y} implies \eqn{y R x}.}
    \item{asymmetric:}{\eqn{x R y} implies that not \eqn{y R x}.}
    \item{antisymmetric}{\eqn{x R y} and \eqn{y R x} imply that
      \eqn{x = y}.}
    \item{transitive:}{\eqn{x R y} and \eqn{y R z} imply that
      \eqn{x R z}.}
    \item{complete:}{for all \eqn{x} and \eqn{y}, \eqn{x R y} or
      \eqn{y R x}.}
  }

  Some combinations of these basic properties have special names
  because of their widespread use:
  \describe{
    \item{preorder:}{reflexive and transitive.}
    \item{quasiorder:}{the same as preorder.}
    \item{equivalence:}{a symmetric preorder.}
    \item{weak order:}{complete and transitive.}
    \item{preference:}{the same as weak order.}
    \item{partial order:}{an antisymmetric preorder.}
    \item{linear order:}{a complete partial order.}
    \item{tournament:}{complete and antisymmetric.}
  }

  If \eqn{R} is a weak order (\dQuote{weak preference relation}),
  \eqn{I = I(R)} defined by \eqn{x I y} iff \eqn{x R y} and \eqn{y R x}
  is an equivalence, the \emph{indifference relation} corresponding to
  \eqn{R}.

  There seem to be no commonly agreed definitions for order relations:
  e.g., Fishburn (1972) requires these to be irreflexive.
}
\references{
  P. C. Fishburn (1972),
  \emph{Mathematics of decision theory}.
  Methods and Models in the Social Sciences 3.
  Mouton: The Hague.

  H. R. Varian (2002),
  \emph{Intermediate Microeconomics: A Modern Approach}.
  6th Edition. W. W. Norton \& Company.
  
}
\keyword{math}
