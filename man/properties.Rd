\name{properties}
\alias{relation_properties}
\alias{relation_property}
%\alias{relation_property<-}
\title{Relation Properties}
\description{Retrieve relation properties.}
\usage{
relation_properties(x)
relation_property(x, which)
%relation_property(x, which) <- value
}
\arguments{
  \item{x}{A relation.}
  \item{which}{Property name (character string).}
%  \item{value}{New value for the property.}
}
\details{
  These functions are used for the \emph{extrinsic} properties of
  relations. Others can be retrieved using the predicate functions.
}
\seealso{
  \code{\link{relation}} and, e.g.,
  \code{\link{relation_is_binary}} for all
  predicate functions defined on relations.
}
\examples{
x <- as.relation(1 : 3)
relation_properties(x)
relation_property(x, "is_endorelation")
% relation_property(x, "is_linear_order")
% relation_property(x, "is_linear_order") <- relation_is_linear_order(x)
% relation_property(x, "is_linear_order")
}
\keyword{math}
