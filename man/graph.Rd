\name{graph}
\alias{relation_graph}
\title{Relation Graph}
\description{
  Determine the graph of a relation.
}
\usage{
relation_graph(x)
}
\arguments{
  \item{x}{an \R object inheriting from class \code{\link{relation}}.}
}
\value{
  An object of class \code{relation_graph}, inheriting from
  \code{\link[sets]{set}}.
}
\seealso{
  \code{\link[sets]{set}()};
  \code{\link{relation}()};
  \code{\link{relation_graph<-}()} for modifying the graph.
}
\examples{
## A simple relation:
R <- as.relation(c(A = 1, B = 2, C = 3))
relation_graph(R)
}
\keyword{math}
