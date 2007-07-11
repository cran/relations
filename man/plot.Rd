\name{plot}
\alias{plot.relation}
\title{Visualize Relations}
\description{
  Visualize certain endorelations by plotting a Hasse Diagram of their 
  transitive reduction. 
}
\usage{
\method{plot}{relation}(x, attrs = list(graph = list(rankdir = "BT"),
                                        edge = list(arrowsize = "0"),
                                        node = list(shape = "rectangle",
                                                    fixedsize = FALSE)),
                        ...)

}
\arguments{
  \item{x}{an \R object inheriting from class \code{\link{relation}}.}
  \item{attrs}{argument passed to the plot method for class
    \code{\link[graph:graphNEL-class]{graphNEL}}.}
  \item{\dots}{Other arguments passed to the
    \code{\link[graph:graphNEL-class]{graphNEL}} plot method.}
}
\details{
  Visualization is currently only provided for either antisymmetric or
  complete transitive endorelations and requires that package \pkg{Rgraphviz} is
  available. In case of transitive complete relations (preferences), the
  dual complement is plotted. Note that the default settings create a diagram with nodes
  ordered bottom-up and with no arrows.
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
