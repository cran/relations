\name{set}
\alias{set}
\alias{as.set}
\alias{is.set}
\alias{set_is_empty}
\alias{set_is_subset}
\alias{set_is_proper_subset}
\alias{set_is_equal}
\alias{set_is_element}
\alias{set_union}
\alias{set_intersection}
\alias{set_symdiff}
\alias{set_complement}
\alias{set_power}
\alias{set_cartesian}
\alias{set_combn}
\alias{\%D\%}
\alias{\%e\%}
\title{Sets}
\description{Creation and manipulation of sets.}
\usage{
set(\dots)
as.set(x)
is.set(x)

set_is_empty(x)
set_is_subset(a, b)
set_is_proper_subset(a, b)
set_is_equal(a, b)
set_is_element(e, b)

set_union(...)
set_intersection(...)
set_symdiff(...)
set_complement(a, b)
set_power(x)
set_cartesian(\dots)
set_combn(x, m)
}
\arguments{
  \item{x}{For \code{as.set()} and \code{is.set()}:
    an \R object. A set object otherwise.}
  \item{a, b}{Set objects.}
  \item{e}{An \R object.}
  \item{m}{Number of elements to choose.}
  \item{\dots}{For \code{set()}: \R objects, and set objects otherwise.}
}
\details{
  These functions represent basic infrastructure for handling sets
  of general (\R) objects. The \code{set_is_\var{foo}()} predicates
  are vectorized. In addition
  to the methods defined, one can use the following operators:
  \code{+} and \code{|} for the union, \code{-} for
  the complement, \code{&} for the
  intersection, \code{\%D\%} for the symmetric difference,
  \code{*} and \code{^n} for the
  (\eqn{n}-fold) cartesian product, \code{2^} for the power set,
  \code{\%e\%} for the element-of predicate,
  \code{<} and \code{<=} for
  the (proper) subset predicate, \code{>} and \code{>=} for
  the (proper) superset predicate, and \code{==} and \code{!=} for
  (in)equality. The \code{length} method for sets gives the
  cardinality. \code{set_combn} returns the set of all
  subsets of specfied length.
  
}
\seealso{
  \code{\link{set_outer}},
  \code{\link{tuple}} for tuples (\dQuote{vectors}), and
  \code{\link{relation}}.
}
\examples{
## constructor
s <- set(1,2,3)
s

## a more complex set
set(c,"test",list(1,2,3))

## converter
s2 <- as.set(2:5)
s2

## set of sets
set(set(), set(1))

## cartesian product
s * s2
s * s
s ^ 2 # same as above
s ^ 3

## power set
2 ^ s

## tuples
s3 <- set(tuple(1,2,3), tuple(2,3,4))
s3

## Predicates:

## element
1:2 \%e\% s
tuple(1,2,3) \%e\% s3

## subset
s <= s2
s2 >= s # same

## proper subset
s < s

## complement, union, intersection, symmetric difference:
s - 1
s + set("a") # or use: s | set("a")
s & s
s \%D\% s2
set(1,2,3) - set(1,2)
set_intersection(set(1,2,3), set(2,3,4), set(3,4,5))
set_union(set(1,2,3), set(2,3,4), set(3,4,5))
set_symdiff(set(1,2,3), set(2,3,4), set(3,4,5))

## subsets:
set_combn(as.set(1:3),2)
}
\keyword{math}
