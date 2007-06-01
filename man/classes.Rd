\name{classes}
\alias{relation_class_ids}
\title{Relation Equivalence Classes}
\description{
  Provide class ids for an equivalence relation, or the indifference
  relation of a weak order.
}
\usage{
relation_class_ids(x)
}
\arguments{
  \item{x}{an object inheriting from class \code{\link{relation}}.}
}
\value{
  A numeric vector with class ids corresponding to the classes of the
  equivalence relation, or the indifference relation of the weak order
  with ids ordered according to increasing preference.
}
\examples{
## Equivalence.
f <- factor(rep(c("Good", "Bad", "Ugly"), c(3, 2, 1)))
R <- as.relation(f)
relation_is_equivalence(R)
table(ids = relation_class_ids(R), orig = f)
## Weak order ("weak preference").
f <- as.ordered(f)
R <- as.relation(f)
relation_is_weak_order(R)
table(ids = relation_class_ids(R), orig = f)
}
\keyword{math}
