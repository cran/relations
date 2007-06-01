\name{table}
\alias{relation_table}
\title{Relation Table}
\description{
  Returns a tabular representation of a relation like a
  \dQuote{view} of a relational database table.
}
\usage{
relation_table(x)
}
\arguments{
  \item{x}{an object inheriting from class \code{\link{relation}}.}
}
\value{
  An object of class \code{relation_table}, inheriting from class
  \code{\link{data.frame}}.
}
\seealso{
  \code{\link{relation_join}}
}
\examples{
R <- data.frame(Name = c("David", "John"),
                Age = c(33, 66),
                stringsAsFactors = FALSE)
R <- as.relation(R)
relation_table(R)
}
\keyword{math}
