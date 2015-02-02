\name{Felines}
\alias{Felines}
\encoding{UTF-8}
\title{Felines Data}
\description{
  A data set with 14 variables on 30 different types of felines.
}
\usage{
data("Felines")
}
\format{
  A data frame with 30 observations on 14 categorical variables (the
  first 10 morphological, the last 4 behavioral), with names in French
  and numeric levels as in the reference.  Names, descriptions in French
  and English and the numbers of levels are as follows.
  \describe{
    \item{\code{TYPPEL}:}{Aspect du pelage; coat; 4.}
    \item{\code{LONGPOIL}:}{Fourrure; fur; 2.}
    \item{\code{OREILLES}:}{Oreilles; ears; 2.}
    \item{\code{TAILLE}:}{Taille au garrot; waist; 3.}
    \item{\code{POIDS}:}{Poids; weight; 3.}
    \item{\code{LONGUEUR}:}{Longueur du corps; body length; 3.}
    \item{\code{QUEUE}:}{Longueur de la queue; tail length; 3.}
    \item{\code{DENTS}:}{Canines développées; carnassials; 2.}
    \item{\code{LARYNX}:}{Os hyaoide; larynx; 2.}    
    \item{\code{RETRACT}:}{Griffes rétractiles; retractible claws; 2.}
    \item{\code{COMPORT}:}{Comportement prédateur; predatory behavior;
      3.} 
    \item{\code{TYPPROIE}:}{Type de la proie; type of prey; 3.}
    \item{\code{ARBRES}:}{Monte ou non aux arbres; climbs trees or not; 2.}
    \item{\code{CHASSE}:}{Chasse (courre ou affut); chases via chivy or
      ambush; 2.}
  }
}
\source{
  F. Marcotorchino and P. Michaud (1982),
  Agregation de similarités en classification automatique.
  \emph{Revue de Statistique Appliquée}, \bold{30}(2), 21--44.
}
\examples{
data("Felines")
summary(Felines)
}
\keyword{datasets}
