\name{classes}
\alias{relation_class_ids}
\alias{relation_classes}
\title{Relation Equivalence Classes}
\description{
  Provide class ids or classes, respectively, for an equivalence
  relation or the indifference relation of a weak order.
}
\usage{
relation_class_ids(x)
relation_classes(x)
}
\arguments{
  \item{x}{an object inheriting from class \code{\link{relation}}
    representing a crisp endorelation.}
}
\value{
  For \code{relation_class_ids}, a numeric vector with class ids
  corresponding to the classes of the equivalence relation, or the
  indifference relation of the weak order with ids ordered according to
  increasing preference.
  
  For \code{relation_classes}, an object of class
  \code{relation_classes_of_objects}, which is a list of sets giving the
  elements in the corresponding classes, named by the class ids.
}
\examples{
## Equivalence.
f <- factor(rep(c("Good", "Bad", "Ugly"), c(3, 2, 1)))
R <- as.relation(f)
relation_is_equivalence(R)
table(ids = relation_class_ids(R), orig = f)

relation_classes(R)

## Weak order ("weak preference").
f <- ordered(f, levels = c("Ugly", "Bad", "Good"))
R <- as.relation(f)
relation_is_weak_order(R)
table(ids = relation_class_ids(R), orig = f)

relation_classes(R)
}
\keyword{math}
