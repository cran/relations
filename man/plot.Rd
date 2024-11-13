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
                  edge = list(arrowsize = NULL),
                  node = list(shape = "rectangle",
                              fixedsize = FALSE)),
     limit = 6L,
     labels = NULL,
     main = NULL,
     type = c("simplified", "raw"),
     ...)

\method{plot}{relation_ensemble}(x,
     attrs = list(list(graph = list(rankdir = "BT"),
                       edge = list(arrowsize = NULL),
                       node = list(shape = "rectangle",
                                   fixedsize = FALSE))),
     type = "simplified",
     limit = 6L,
     labels = NULL,
     ..., 
     layout = NULL, main = NULL)
}
\arguments{
  \item{x}{an \R object inheriting from class \code{\link{relation}} or
     \code{\link{relation_ensemble}}.}
  \item{attrs}{argument passed to the plot method for class
    \code{\link[graph]{graphNEL}}.
    
    Note that for the \code{relation_ensemble} method, it is a
    \emph{list} of such objects, recycled as needed.}
  \item{type}{character vector of either "simplified" or "raw"
     strings, one for each relation plotted. (See details.)}
  \item{limit}{Argument passed to the labeling function creating default
    labels for the nodes (see \code{\link[sets]{LABELS}()}).
    Recycled as needed for relation ensembles.}
  \item{labels}{Optional list of character vectors defining unique
    labels for the nodes. List of such lists for relation ensembles.}
  \item{layout}{integer vector of length 2 specifying the number of rows
     and columns of the screen layout.  If \code{NULL}, the layout is
     square.}
  \item{\dots}{Other arguments passed to the
    \code{\link[graph]{graphNEL}} plot method.}
  \item{main}{character vector used for the main title(s). If \code{NULL},
     the title(s) is (are) set to the type of the visualized relation(s).}
}
\details{
  Visualization requires that package \pkg{Rgraphviz} is available. If
     type is "simplified" (default), the transitive reduction is
     first computed to reduce visual complexity (especially for transitive
     relations).
     For partial orders, a Hasse diagram is plotted.  In case
     of (acyclic) transitive complete relations (i.e., weak orders, preferences), the dual is plotted.
  For all other acyclic relations, the asymmetric part is plotted. (Note that
  the default settings in these cases create a diagram with nodes ordered bottom-up and
     with no arrows.) For cyclic relations, a raw graph (with arrows) of
     the corresponding transitive
     reduction is computed. If type is "raw", a directed graph without
     any transformation is plotted from the relation.
}
\seealso{
  \code{\link{relation}()},
  \code{\link{transitive_reduction}()}
}
\examples{
require("sets")				# set() etc.
if(require("Rgraphviz")) {
  ## simple example
  plot(as.relation(1 : 5))

  ## inclusion on a power set:
  ps <- 2 ^ set("a", "b", "c")
  inc <- set_outer(ps, set_is_subset)
  R <- relation(incidence = inc)
  plot(relation_ensemble(R, R), type = c("simplified", "raw"))
}
}
\keyword{math}
