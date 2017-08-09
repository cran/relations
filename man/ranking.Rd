\name{ranking}
\alias{ranking}
\alias{as.ranking}
\alias{is.ranking}
\title{Rankings}
\description{
  Creates a ranking object.
}
\usage{
ranking(x, domain = NULL, decreasing = TRUE, complete = FALSE)
as.ranking(x, \dots)
is.ranking(x)
}
\arguments{
  \item{x}{For \code{ranking()}: either an atomic vector interpreted as
    labels of the ranked objects, or as their scores in case of a named
    numeric vector (with the names giving the labels), or a list of
    atomic vectors representing equivalence classes of labels.
    For \code{as.ranking()}: an \R object coercible to a ranking object
    (including \code{relation} objects).}
  \item{domain}{object coercible to a set, from which the labels usable
    in \code{x} are derived.  If \code{NULL}, it is created from
    \code{x}.}
  \item{decreasing}{logical indicating whether the ranking orders objects
    from the best to the worst (\code{TRUE}), or the other way round.}
  \item{complete}{logical specifying whether missing values should be
    imputed, if any. Missing elements are those from \code{domain}
    not used in \code{x}. If \code{decreasing} is \code{TRUE}
    (\code{FALSE}), all missings are ranked tied behind (ahead) the
    worst (best) ranked object.}
  \item{\dots}{currently not used.}
}
\value{
  An object of class \code{ranking}.
}
\seealso{
  \code{\link{relation}()}
}
\examples{
## simple rankings
OBJECTS <- c("Apples", "Bananas", "Oranges", "Lemons")
print(R <- ranking(OBJECTS))
ranking(OBJECTS[2:4], domain = OBJECTS)
ranking(OBJECTS[2:4], domain = OBJECTS, complete = TRUE)

## ranking with ties (weak orders)
ranking(list(c("PhD", "MD"), "MSc", c("BSc", "BA")))

## ranking A > B ~ C with D missing:
ranking(c(A = 1, B = 2, C = 2, D = NA))

## coercion functions
identical(as.ranking(as.relation(R)), R)
}
\keyword{math}
