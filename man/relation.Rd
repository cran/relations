\name{relation}
\alias{relation}
\alias{endorelation}
\alias{homorelation}
\alias{as.relation}
\alias{is.relation}
\alias{codual}
\alias{codual.relation}
\alias{sym}
\alias{sym.relation}
\alias{asy}
\alias{asy.relation}
\title{Relations}
\description{Creation and manipulation of relations.}
\usage{
relation(domain = NULL, incidence = NULL, graph = NULL,
         charfun = NULL)
endorelation(domain = NULL, incidence = NULL, graph = NULL,
             charfun = NULL)
homorelation(domain = NULL, incidence = NULL, graph = NULL,
             charfun = NULL)
as.relation(x, \dots)
is.relation(x)
}
\arguments{
  \item{domain}{List (or tuple) of (possibly named) sets (or vectors)
    used as the domain, recycled as needed to fit the arity of the relation.
    If \code{domain} is not a list or tuple, it is converted to a list.}
  \item{incidence}{A numeric array with values in the unit interval, or
    a logical array.  Note that one-dimensional incidences are also
    accepted.  The \code{names}/\code{dimnames} attribute of
    \code{incidence} is used as domain if this is not explicitly given
    using the \code{domain} argument.}
  \item{graph}{Either a set of equally sized tuples, or a list of
    (possibly, generic) vectors of same length where each component
    specifies one relation element, or a data frame where each row
    specifies one relation element.  For the latter, the columns
    correspond to the domain sets, and the colnames are used as their
    labels if specified.}
  \item{charfun}{A characteristic function of the relation, i.e., a
    predicate function taking \eqn{k} arguments, with \eqn{k} equal to
    the arity of the relation.}
  \item{x}{an \R object.}
  \item{\dots}{Further arguments passed to \code{as.relation()}
    methods (currently not used for those defined in the relations
    package).}
}
\details{
  Given \eqn{k} sets of objects \eqn{X_1}, \ldots, \eqn{X_k}, a
  \eqn{k}-ary relation \eqn{R} on \eqn{D(R) = (X_1, \ldots, X_k)} is a
  (possibly fuzzy) subset \eqn{G(R)} of the Cartesian product \eqn{X_1
    \times \cdots \times X_k}{X_1 x ... x X_k}.  We refer to \eqn{D(R)}
  and \eqn{G(R)} as the \emph{domain} and the \emph{graph} of the
  relation, respectively (alternative notions are that of \emph{ground}
  and \emph{figure}, respectively).  We also refer to \eqn{s = (s_1,
    \ldots, s_k)}, where each \eqn{s_i} gives the cardinality of
  \eqn{X_i}, as the \emph{size} of the relation.

  Strictly speaking, the relation is the \emph{pair} \eqn{(D(R), G(R))};
  often, relations are identified with their graph.  If \eqn{G(R)} is a
  crisp subset of \eqn{D(R)}, \eqn{R} is a \emph{crisp relation}.  In
  this case, we say that a \eqn{k}-tuple \eqn{t} is \emph{contained} in
  the relation \eqn{R} iff it is an element of \eqn{G(R)}.

  The \emph{characteristic function} \eqn{f_R} of a relation \eqn{R} is
  the membership function of \eqn{G(R)}, giving for each \eqn{k}-tuple
  \eqn{t} in \eqn{D(R)} the membership (amount of belongingness) of
  \eqn{t} to \eqn{G(R)}.  In the crisp case, \eqn{f_R} is also referred
  to as the indicator function of the relation, and is a binary (0/1)
  function such that \eqn{f_R(t)} is one iff \eqn{t} is in \eqn{G(R)}.

  Relations with arity 2, 3, and 4 are typically referred to as
  \emph{binary}, \emph{ternary}, and \emph{quaternary} relations,
  respectively.  A \emph{homorelation} on \eqn{X} is a relation with
  homogeneous domain, i.e. \eqn{(X, X, \dots, X)}.
  An \emph{endorelation} on \eqn{X} (or binary relation
  \emph{over} \eqn{X}) is a binary homorelation.
  See \link[=relation_is]{predicates} for the most important
  types of endorelations.
  
  Relations with the same domain can naturally be ordered according to
  their graphs.  I.e., \eqn{R \le S} iff \eqn{G(R)} is a subset of
  \eqn{G(S)}, or equivalently, iff \eqn{f_R(t) \le f_S(t)} for every
  \eqn{k}-tuple \eqn{t} (in the crisp case, iff every tuple contained in
  \eqn{R} is also contained in \eqn{S}).  This induces a lattice
  structure, with meet (greatest lower bound) and join (least upper
  bound) the intersection and union of the graphs, respectively, also
  known as the \emph{intersection} and \emph{union} of the relations.
  The least element moves metric on this lattice is the
  \emph{symmetric difference metric}, i.e., the Manhattan distance
  between the collections of membership values (incidences).  In the
  crisp case, this gives the cardinality of the symmetric difference of
  the graphs (the number of tuples in exactly one of the relation
  graphs).

  The \emph{complement} (or \emph{negation}) \eqn{R^c} of a relation
  \eqn{R} is the relation with domain \eqn{D(R)} whose graph is the
  complement of \eqn{G(R)} (in the crisp case, containing exactly the
  tuples not contained in \eqn{R}).

  For binary crisp relations \eqn{R}, it is customary to write
  \eqn{x R y} iff \eqn{(x, y)} is contained in \eqn{R}.  For binary
  crisp relations \eqn{R_1} and \eqn{R_2} with domains \eqn{(X, Y)} and 
  \eqn{(Y, Z)}, the \emph{composition} \eqn{T = R_1 * R_2} of \eqn{R_1}
  and \eqn{R_2} is defined by taking \eqn{x S z} iff there is an \eqn{y}
  such that \eqn{x R_1 y} and \eqn{y R_2 z}.  The \emph{transpose} (or
  \emph{inverse}) \eqn{R^{t}} of the relation \eqn{R} with domain
  \eqn{(X, Y)} is the relation with domain \eqn{(Y, X)} such that
  \eqn{y R^{t} x} iff \eqn{x R y}.  The 
  \emph{codual} (Clark (1990), also known as the \sQuote{dual} in the
  fuzzy preference literature, e.g., Ovchinnikov (1991)) is the
  complement of the transpose (equivalently, the transpose of the
  complement).

  For binary fuzzy relations \eqn{R}, one often writes \eqn{R(x, y)} for
  the membership of the pair \eqn{(x, y)} in the relation.  The above
  notions need to take the fuzzy logic employed (as described by the
  fuzzy t-norm (intersection) \eqn{T}, t-conorm (disjunction) \eqn{S},
  and negation \eqn{N}) into account.  Let \eqn{R}, \eqn{R_1} and
  \eqn{R_2} be binary relations with appropriate domains.  Then the
  memberships for \eqn{(x, y)} of the complement, transpose and codual of
  \eqn{R} are \eqn{N(R(x, y))}, \eqn{R(y, x)} and \eqn{N(R(y, x))},
  respectively.  The membership of \eqn{(x, y)} for the composition of
  \eqn{R_1} and \eqn{R_2} is \eqn{\max_z T(R_1(x, z), R_2(z, y))}.
  
  Package \pkg{relations} implements finite relations as an S3 class
  which allows for a variety of representations (even though currently,
  typically dense array representations of the incidences are employed).
  Other than by the generator,
  relations can be obtained by coercion via the generic function
  \code{as.relation()}, which has methods for at least logical and numeric
  vectors, unordered and ordered factors, arrays including matrices, and
  data frames.  Unordered factors are coerced to equivalence relations;
  ordered factors and numeric vectors are coerced to order relations.
  Logical vectors give unary relations (predicates).  A (feasible)
  \eqn{k}-dimensional array is taken as the incidence of a \eqn{k}-ary
  relation.  Finally, a data frame is taken as a relation table.  Note
  that missing values will be propagated in the coercion.

  \code{endorelation()} is a wrapper for \code{relation()}, trying to
  guess a suitable domain from its arguments to create an
  endorelation.  If a domain is given, all labels are combined and the
  result (as a list) recycled as needed.
  
  Basic relation operations are available as group methods: \code{min()}
  and \code{max()} give the meet and join, and \code{range()} a
  \link[=relation_ensemble]{relation ensemble} with these two.
  Comparison operators implement the natural ordering in the relation
  lattice.  Where applicable, \code{!} gives the complement (negation),
  \code{&} and \code{|} intersection and union, and \code{*}
  composition, respectively.  Finally, \code{t()} gives the transpose
  and \code{codual()} the codual.

  There is a \code{\link[=plot.relation]{plot}()} method for certain
  crisp endorelations provided that package \pkg{Rgraphviz} is
  available.

  For crisp endorelations \eqn{R}, \code{sym()} and \code{asy()} give
  the symmetric and asymmetric parts of \eqn{R}, defined as the
  intersection of \eqn{R} with its transpose, or, respectively, with its
  codual (the complement of its transpose).

  The \code{summary()} method applies all \link{predicates} available
  and returns a logical vector with the corresponding results.
}
\references{
  S. A. Clark (1990),
  A folk meta-theorem in the foundations of utility theory.
  \emph{Mathematical Social Sciences}, \bold{19}/3, 253--267.
  \doi{10.1016/0165-4896(90)90065-F}.
  
  S. Ovchinnikov (1991),
  Similarity relations, fuzzy partitions, and fuzzy orderings.
  \emph{Fuzzy Sets and Systems}, \bold{40}/1, 107--126.
  \doi{10.1016/0165-0114(91)90048-U}.
}
\seealso{
  \code{\link{relation_incidence}()} for obtaining incidences;
  \code{\link{relation_domain}()} for determining domain, arity, and
  size;
  \code{\link{relation_graph}()} for determining the graph of a relation;
  \code{\link{relation_charfun}()} for determining the characteristic
  function;
  \link{predicates} for available predicate functions; and
  \link{algebra} for further operations defined on relations.
}
\examples{
require("sets")				# set(), tuple() etc.

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
R <- relation(graph = set(tuple(1,2), tuple(1,3),
                          tuple(2,4), tuple(3,4)))
## equivalent to:
## relation(graph = list(c(1,2), c(1,3), c(2,4), c(3,4)))
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

### endorelation
s <- set(pair("a","b"), pair("c","d"))
relation_incidence(relation(graph = s))
relation_incidence(endorelation(graph = s))
}
\keyword{math}
