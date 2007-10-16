\name{choice}
\alias{relation_choice}
\title{Relation-Based Choices}
\description{
  Choose objects based on an ensemble of relations between these.
}
\usage{
relation_choice(x, method = "symdiff", control = list(), ...)
}
\arguments{
  \item{x}{an ensemble of endorelations.}
  \item{method}{a character string specifying one of the built-in
    methods, or a function to be taken as a user-defined method.
    See \bold{Details} for available built-in methods.}
  \item{control}{a list of control parameters.  See \bold{Details}.}
  \item{\dots}{a list of control parameters (overruling those specified
    in \code{control}).}
}
\value{
  A set with the chosen objects, or a list of such sets.
}
\details{
  A social choice function is a rule for \emph{choosing} from a set
  \eqn{D} of objects, i.e., selecting suitable subsets of \eqn{D}.
  Voting rules used in elections are the most prominent example of such
  functions, which typically aggregate individual preferences (e.g., of
  voters).

  Choice method \code{"symdiff"} (currently the only one available)
  chooses a given number \eqn{k} of objects (\dQuote{winners}) by
  determining a relation \eqn{R} which minimizes \eqn{\sum_b d(R_b, R)}
  over all relations for which winners are always strictly preferred to
  losers, without any further constraints on the relations between pairs
  of winners or pairs of losers, where \eqn{d} is symmetric difference
  dissimilarity and the \eqn{R_b} are crisp endorelations.  (Note that
  this is different from computing consensus preference relations.)

  Available control options include:
  \describe{
    \item{\code{k}}{an integer giving the number of objects/winners to
      be chosen.}
    \item{\code{all}}{a logical indicating whether all optimal choices
      should be obtained.  By default, or if \code{all} is false, only a
      single optimal choice is computed.}
  }
}
\examples{
data("SVM_Benchmarking_Classification")
## Determine the three best classification learners in the above sense.
relation_choice(SVM_Benchmarking_Classification, k = 3)
}
\keyword{math}
