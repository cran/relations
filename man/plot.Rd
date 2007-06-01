\name{plot}
\alias{plot.relation}
\title{Visualize Relations}
\description{
  Visualize certain endorelations by plotting a Hasse Diagram of their 
  transitive reduction. 
}
\usage{
\method{plot}{relation}(x, attrs = list(edge = list(arrowsize = "0")), ...)
}
\arguments{
  \item{x}{an \R object inheriting from class \code{\link{relation}}.}
  \item{attrs}{argument passed to the plot method for class
    \code{\link[graph:graphNEL-class]{graphNEL}}.}
  \item{\dots}{Other arguments passed to the
    \code{\link[graph:graphNEL-class]{graphNEL}} plot method.}
}
\details{
  Visualization is currently only provided for antisymmetric and
  transitive endorelations, and requires that package \pkg{Rgraphviz} is
  available.
}
\seealso{
  \code{\link{relation}}
}
\examples{
if(require("Rgraphviz")) {
  ## simple example
  plot(as.relation(1 : 5))

  ## inclusion on a power set:
  ps <- 2 ^ set("a", "b", "c")
  inc <- set_outer(ps, set_is_subset)
  plot(relation(incidence = inc))
}
}
\keyword{math}
