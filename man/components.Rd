\name{components}
\alias{relation_connected_components}
\alias{relation_condensation}
\alias{relation_component_representation}
\title{Connected components}
\description{
  Computes (strongly or weakly) connected components of an endorelation.
}
\usage{
relation_connected_components(x, type = c("strongly", "weakly"))
relation_condensation(x)
relation_component_representation(x)
}
\arguments{
  \item{x}{an \R object inheriting from class \code{\link{relation}},
    representing a crisp endorelation without missings.}
  \item{type}{character string indicating the kind of components sought.}
}
\value{
  For \code{relation_connected_components()}, an object of class
  \code{relation_classes_of_objects}, i.e., a list of sets giving the
  elements of the corresponding connected components, named by the
  leaders' character representation. The list of leaders is added as a
  \code{leaders} attribute.

  For \code{relation_condensation()}, an (acyclic) endorelation.

  For \code{relation_component_representation()}, an endorelation with
  same domain as \code{x}.
  
}
\details{
  Let \eqn{G} be the graph of an endorelation \eqn{R}.
  
  A \emph{weakly} connected component of some node \eqn{k} in \eqn{G} is
  the set of all nodes reachable from \eqn{k}.  A \emph{strongly}
  connected component of some node \eqn{k} is the set of all nodes
  reachable from \eqn{k}, from which \eqn{k} can be reached. Each
  component is represented by some element, the \emph{leader}.

  The \emph{component representation} graph of a cyclic endorelation
  \eqn{R} is composed of directed cycles, one for each strongly
  connected component of \eqn{R} containing more than one element,
  linking all corresponding elements.

  The \emph{condensation} of \eqn{R} is the graph of all leaders of
  \eqn{R}.
}
\seealso{
  \code{\link{plot.relation}()},
  \code{\link{transitive_reduction}()}
}
\examples{
## example from La Poutre and van Leeuwen:

require("sets")				# set(), pair() etc.

G <- set(pair(1L, 2L), pair(2L, 1L), pair(1L, 3L), pair(3L, 1L),
         pair(3L, 7L), pair(2L, 5L), pair(2L, 6L), pair(6L, 5L),
         pair(5L, 7L), pair(4L, 6L), pair(5L, 4L), pair(4L, 7L))

R <- endorelation(graph = G)

relation_connected_components(R)
relation_graph(relation_condensation(R))
relation_graph(relation_component_representation(R))
}
\references{
  S. Warshall (1962),
  A theorem on Boolean matrices.
  \emph{Journal of the ACM}, \bold{9}/1, 11--12.
  \doi{10.1145/321105.321107}.

  J. A. La Poutré and J. van Leeuwen (1988),
  Maintenance of Transitive Closures and Transitive Reductions of
  Graphs.
  \emph{Proceedings of the International Workshop of Graph-Theoretic
  Concepts in Computer Science}, Springer, London, 106--120.
}
\keyword{math}
