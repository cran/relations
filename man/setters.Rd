\name{setters}
\alias{relation_domain<-}
\alias{relation_domain_names<-}
\alias{relation_graph<-}
\alias{relation_incidence<-}
\title{Modify Relations}
\description{
  Modify relations by (re)setting their domain, graph, or incidences.
}
\usage{
relation_domain(x) <- value
relation_domain_names(x) <- value
relation_graph(x) <- value
relation_incidence(x) <- value
}
\arguments{
  \item{x}{an \R object inheriting from class \code{\link{relation}}.}
  \item{value}{for setting the domain, a tuple (or list) as long as the
    arity of the relation \code{x}, with sets of cardinality (for lists:
    numbers of elements) identical to the size of \code{x}.
    
    For setting the graph, either a set of tuples of equal lengths
    (arity of the relation) or a data frame or something coercible to
    this, with the values of the components of the given tuples (rows)
    always elements of the corresponding elements of the domain of
    \code{x}.
    
    For setting incidences, a numeric array with values in the unit
    interval or a logical array, with dimension the size of the relation
    \code{x}.

    For setting the domain names, a character vector as long as the
    arity of the relation \code{x}.
  }
}
\seealso{
  \code{\link{relation_domain}} for getting the domain of a relation;
  \code{\link{relation_domain_names}} for getting the domain names;
  \code{\link{relation_graph}} for getting the graph;
  \code{\link{relation_incidence}} for getting the incidences;
  \code{\link{relation}} for basic information.
}
\examples{
R <- as.relation(1 : 3)
print(R)

relation_domain(R)
## tuple format:
relation_domain(R) <- pair(X = set("a","b","c"), Y = set("A","B","C"))
relation_domain(R)
## the same in list format:
relation_domain(R) <- list(X = letters[1:3], Y = LETTERS[1:3])
relation_domain(R)

relation_domain_names(R) <- c("XX","YY")
relation_domain_names(R)

relation_incidence(R)
relation_incidence(R) <- diag(1, 3, 3)
relation_incidence(R)

relation_graph(R)
## set format:
relation_graph(R) <- set(pair("a","B"), pair("a","C"), pair("b","C"))
relation_graph(R)
## the same in data frame format:
relation_graph(R) <- data.frame(c("a", "a", "b"), c("B", "C", "C"))
relation_graph(R)

}
\keyword{math}
