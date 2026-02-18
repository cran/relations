\name{Cetacea}
\alias{Cetacea}
\title{Cetacea Data}
\description{
  A data set with 16 variables on 36 different types of cetacea.
}
\usage{
data("Cetacea")
}
\format{
  A data frame with 36 observations on 16 categorical variables.  The
  first 15 variables relate to morphology, osteology, or behavior, and
  have both self-explanatory names and levels.  The last (\code{CLASS})
  gives a common zoological classification.
}
\source{
  \bibshow{Vescia:1985}
}
\examples{
data("Cetacea")
summary(Cetacea)
## Show the cetacea types grouped by class.
split(rownames(Cetacea), Cetacea$CLASS)
}
\keyword{datasets}
