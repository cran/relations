\name{plot}
\alias{plot.relation}
\alias{plot.relation_ensemble}
\title{Visualize Relations}
\description{
  Visualize certain crisp endorelations by plotting a Hasse Diagram of
  their transitive reduction. 
}
\usage{
\method{plot}{relation}(x,
     attrs = list(graph = list(rankdir = "BT"),
                               edge = list(arrowsize = "0"),
                               node = list(shape = "rectangle",
                                           fixedsize = FALSE)),
     limit = 6L,
     labels = NULL,
     main = NULL,
     ...)

\method{plot}{relation_ensemble}(x,
     attrs = list(list(graph = list(rankdir = "BT"),
                                    edge = list(arrowsize = "0"),
                                    node = list(shape = "rectangle",
                                                fixedsize = FALSE))),
     ..., layout = NULL, main = NULL)
}
\arguments{
  \item{x}{an \R object inheriting from class \code{\link{relation}} .}
  \item{attrs}{argument passed to the plot method for class
    \code{\link[graph:graphNEL-class]{graphNEL}}. For the
    \code{relation_ensemble} method, it is a \emph{list} of such objects,
    recycled as needed.}
  \item{limit}{Argument passed to the labeling function creating default
    labels for the nodes (see \code{\link[sets]{LABELS}}).}
  \item{labels}{Optional list of character vectors defining unique
    labels for the nodes.}
  \item{layout}{integer vector of length 2 specifying the number of rows
     and columns of the screen layout. If \code{NULL}, the layout is sqare.}
  \item{\dots}{Other arguments passed to the
    \code{\link[graph:graphNEL-class]{graphNEL}} plot method.}
  \item{main}{character vector used for the main title(s). If \code{NULL},
     the title(s) is (are) set to the type of the visualized relation(s).}
}
\details{
  Visualization
%  is currently only provided for either antisymmetric or
%  complete transitive crisp endorelations and
requires that package \pkg{Rgraphviz} is available.
For partial orders, a Hasse diagram is plotted.
In case of transitive complete
relations (weak orders, preferences), the dual is plotted. For all
other relations, the asymmetric part is plotted. Note that the
default settings create a diagram with nodes ordered bottom-up and
with no arrows.
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
