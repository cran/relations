\name{relation}
\alias{relation}
\alias{as.relation}
\alias{is.relation}
\title{Relations}
\description{Creation and manipulation of relations.}
\usage{
relation(domain = NULL, incidence = NULL, graph = NULL, charfun = NULL)
as.relation(x)
is.relation(x)
}
\arguments{
  \item{domain}{List (or tuple) of (possibly named) sets (or vectors)
    used as the domain, recycled as needed to fit the arity of the relation.
    If \code{domain} is not a list or tuple, it is converted to a list.}
  \item{incidence}{Array of 0/1 (or \code{TRUE}/\code{FALSE})
    values.  Note that one-dimensional incidences are also accepted.
    The \code{names}/\code{dimnames} attribute of \code{incidence} is
    used as domain if this is not explicitly given using the
    \code{domain} argument.}
  \item{graph}{Either a set of equally sized tuples, or a list of (possibly,
    generic) vectors of same length where each component specifies one relation
    element, or a data frame where each row specifies one relation
    element.  For the latter, the columns correspond to the domain sets, and the
    colnames are used as their labels if specified.}
  \item{charfun}{A characteristic function of the relation, i.e., a
    predicate function taking \eqn{k} arguments, with \eqn{k} equal to
    the arity of the relation.}
  \item{x}{an \R object.}
}
\details{
  Given \eqn{k} sets of objects \eqn{X_1}, \ldots, \eqn{X_k}, a
  \eqn{k}-ary relation \eqn{R} on \eqn{D(R) = (X_1, \ldots, X_k)} is a
  subset \eqn{G(R)} of the Cartesian product \eqn{X_1 \times \cdots
    \times X_k}{X_1 x ... x X_k}.  We refer to \eqn{D(R)} and \eqn{G(R)}
  as the \emph{domain} and the \emph{graph} of the relation,
  respectively (alternative notions are that of \emph{ground} and
  \emph{figure}, respectively).  We also refer to \eqn{s = (s_1, \ldots,
    s_k)}, where each \eqn{s_i} gives the cardinality of \eqn{X_i}, as
  the \emph{size} of the relation.

  Strictly speaking, the relation is the \emph{pair} \eqn{(D(R), G(R))};
  often, relations are identified with their graph.  We say that a
  \eqn{k}-tuple \eqn{t} is \emph{contained} in the relation \eqn{R} iff
  it is an element of \eqn{G(R)}.

  The \emph{characteristic function} \eqn{f_R} (sometimes also referred
  to as indicator function) of a relation \eqn{R} is the predicate
  (boolean-valued) function on the Cartesian product \eqn{X_1 \times
    \cdots \times X_k}{X_1 x ... x X_k} such that \eqn{f_R(t)} is true
  iff the \eqn{k}-tuple \eqn{t} is in \eqn{G(R)}.

  Relations with arity 2, 3, and 4 are typically referred to as
  \emph{binary}, \emph{ternary}, and \emph{quaternary} relations,
  respectively.  An \emph{endorelation} on \eqn{X} (or binary relation
  \emph{over} \eqn{X}) is a binary relation with domain \eqn{(X, X)}.
  See \link[relation_is_endorelation]{predicates} for the most important
  types of endorelations.
  
  Relations with the same domain can naturally be ordered according to
  their graphs.  I.e., \eqn{R_1 \le R_2} iff \eqn{G(R_1)} is a subset of
  \eqn{G(R_2)}, or equivalently, if every \eqn{k}-tuple \eqn{t}
  contained in \eqn{R_1} is also contained in \eqn{R_2}.  This induces a
  lattice structure, with meet (greatest lower bound) and join (least
  upper bound) the intersection and union of the graphs, respectively,
  also known as the \emph{intersection} and \emph{union} of the
  relations.  The least element moves metric on this lattice is the
  \emph{symmetric difference metric}, i.e., the cardinality of the
  symmetric difference of the graphs (the number of tuples in exactly
  one of the relation graphs).

  The \emph{complement} of a relation \eqn{R} is the relation with
  domain \eqn{D(R)} whose graph is the complement of \eqn{G(R)}, i.e.,
  which contains exactly the tuples not contained in \eqn{R}.

  For binary relations \eqn{R}, it is customary to write \eqn{x R y} iff 
  \eqn{(x, y)} is contained in \eqn{R}.  For binary relations \eqn{R}
  and \eqn{S} with domains \eqn{(X, Y)} and \eqn{(Y, Z)}, the
  \emph{composition} of \eqn{R} and \eqn{S} is defined by taking
  \eqn{x S z} iff there is an \eqn{y} such that \eqn{x R y} and
  \eqn{y S z}.  The \emph{dual} (or \emph{converse}) \eqn{R^*} of the
  relation \eqn{R} with domain \eqn{(X, Y)} is the relation with domain
  \eqn{(Y, X)} such that \eqn{x R^* y} iff \eqn{y R x}.

  Package \pkg{relations} implements finite relations as an S3 class
  which allows for a variety of representations (even though currently,
  only dense array representations of the incidences are employed).
  Other than by the generator,
  relations can be obtained by coercion via the generic function
  \code{as.relation}, which has methods for at least logical and numeric
  vectors, unordered and ordered factors, arrays including matrices, and
  data frames.  Unordered factors are coerced to equivalence relations;
  ordered factors and numeric vectors are coerced to order relations.
  Logical vectors give unary relations (predicates).  A (feasible)
  \eqn{k}-dimensional array is taken as the incidence of a \eqn{k}-ary
  relation.  Finally, a data frame is taken as a relation table.  Note
  that missing values will be propagated in the coercion.
  
  Basic relation operations are available as group methods: \code{min}
  and \code{max} give the meet and join, and \code{range} a
  \link[relation_ensemble]{relation ensemble} with these two.
  Comparison operators implement the natural ordering in the relation
  lattice.  Where applicable, \code{!} gives the complement, \code{&}
  and \code{|} intersection and union, and \code{*} composition,
  respectively.  Finally, \code{t} gives the dual.

  There is a \code{\link[plot.relation]{plot}} method for certain
  endorelations provided that package \pkg{Rgraphviz} is available.

  There is also the notion of fuzzy relations, for which each tuple is
  contained in the graph with a certain membership value.  (I.e., the
  graph is a fuzzy subset of the Cartesian product of the elements of
  the domain.)  Basic support for fuzzy relations will be added
  eventually.
}
\seealso{
  \code{\link{relation_incidence}} for obtaining incidences;
  \code{\link{relation_domain}} for determining domain, arity, and
  size;
  \code{\link{relation_graph}} for determining the graph of a relation;
  \code{\link{relation_charfun}} for determining the characteristic
  function;
  \code{\link{algebra}} for further operations defined on relations.
}
\examples{
## A relation created by specifying the graph:
R <- relation(graph = data.frame(A = c(1, 1:3), B = c(2:4, 4)))
relation_incidence(R)
## extract domain
relation_domain(R)
## extract graph
relation_graph(R)
## both ("a pair of domain and graph" ...)
as.tuple(R)

## (Almost) the same using the set specification
## (the domain labels are missing).
R <- relation(graph = set(tuple(1,2), tuple(1,3), tuple(2,4), tuple(3,4)))
## equivalent to:
## relation(graph = list(1:2, c(1,3), c(2,4), c(3,4)))
relation_incidence(R)

## Explicitly specifying the domain:
R <- relation(domain = list(A = letters[1:3], B = LETTERS[1:4]),
              graph = set(tuple("a","B"), tuple("a","C"),
                          tuple("b","D"), tuple("c","D")))
relation_incidence(R)

## Domains can be composed of arbitrary R objects:
R <- relation(domain = set(c, "test"),
              graph = set(tuple(c, c), tuple(c, "test")))
relation_incidence(R)

## Characteristic function ("a divides b"):
R <- relation(domain = list(1 : 10, 1 : 10),
              charfun = function(a, b) b \%\% a == 0)
relation_incidence(R)
## R is a partial order: plot the Hasse diagram provided that
## Rgraphviz is available:
if(require("Rgraphviz")) plot(R)

## conversions and operators
x <- matrix(0, 3, 3)
R1 <- as.relation(row(x) >= col(x))
R2 <- as.relation(row(x) <= col(x))
R3 <- as.relation(row(x) <  col(x))
relation_incidence(max(R1, R2))
relation_incidence(min(R1, R2))
R3 < R2
relation_incidence(R1 * R2)
relation_incidence(! R1)
relation_incidence(t(R2))
}
\keyword{math}
