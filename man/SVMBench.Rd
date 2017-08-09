\name{SVMBench}
\alias{SVM_Benchmarking_Classification}
\alias{SVM_Benchmarking_Regression}
\alias{SVM_Benchmarking_Classification_Consensus}
\alias{SVM_Benchmarking_Regression_Consensus}
\title{SVM Benchmarking Data and Consensus Relations}
\description{
  \code{SVM_Benchmarking_Classification} and
  \code{SVM_Benchmarking_Regression} represent the
  results of a benchmark study (Meyer, Leisch and Hornik, 2003)
  comparing Support Vector Machines to 
  other predictive methods on real and artificial data sets
  involving classification and regression
  methods, respectively.
  In addition,
  \code{SVM_Benchmarking_Classification_Consensus}
  and \code{SVM_Benchmarking_Regression_Consensus}
  provide consensus rankings derived from these data.
}
\usage{
data("SVM_Benchmarking_Classification")
data("SVM_Benchmarking_Regression")
data("SVM_Benchmarking_Classification_Consensus")
data("SVM_Benchmarking_Regression_Consensus")
}
\format{
  \code{SVM_Benchmarking_Classification}
  (\code{SVM_Benchmarking_Regression}) is an ensemble of 21 (12)
  relations representing pairwise comparisons of 17 classification (10
  regression) methods on 21 (12) data sets.  Each relation of the 
  ensemble summarizes the results for a particular data set.  The
  relations are reflexive endorelations on the set of methods
  employed, with a pair \eqn{(a, b)} of distinct methods contained in a
  relation iff both delivered results on the corresponding data set and
  \eqn{a} did not perform significantly better than \eqn{b} at the 5\%
  level.  Since some methods failed on some data sets, the relations are
  not guaranteed to be complete or transitive.  See Meyer et al. (2003)
  for details on the experimental design of the benchmark study, and
  Hornik and Meyer (2007) for the pairwise comparisons.

  The corresponding consensus objects are lists of ensembles of
  consensus relations fitted to the benchmark results.
  For each of the following three endorelation families:
  \code{SD/L} (\dQuote{linear orders}),
  \code{SD/O} (\dQuote{partial orders}), and \code{SD/W}
  (\dQuote{weak orders}), \emph{all} possible consensus
  relations have been computed (see \code{\link{relation_consensus}}).
  For both classification and regression,
  the three relation ensembles obtained are provided as a named list of
  length 3. See Hornik & Meyer (2007) for details on the meta-analysis.
}
\source{
  D. Meyer, F. Leisch, and K. Hornik (2003),
  The support vector machine under test.
  \emph{Neurocomputing}, \bold{55}, 169--186.
  \doi{10.1016/S0925-2312(03)00431-4}.

  K. Hornik and D. Meyer (2007),
  Deriving consensus rankings from benchmarking experiments.
  In R. Decker and H.-J. Lenz,
  \emph{Advances in Data Analysis}.
  Studies in Classification, Data Analysis, and Knowledge Organization.
  Springer-Verlag: Heidelberg, 163--170.  
}
\examples{
data("SVM_Benchmarking_Classification")

## 21 data sets
names(SVM_Benchmarking_Classification)

## 17 methods
relation_domain(SVM_Benchmarking_Classification)

## select weak orders
weak_orders <-
    Filter(relation_is_weak_order, SVM_Benchmarking_Classification)

## only the artifical data sets yield weak orders
names(weak_orders)

## visualize them using Hasse diagrams
if(require("Rgraphviz")) plot(weak_orders)

## Same for regression:
data("SVM_Benchmarking_Regression")

## 12 data sets
names(SVM_Benchmarking_Regression)

## 10 methods
relation_domain(SVM_Benchmarking_Regression)

## select weak orders
weak_orders <-
    Filter(relation_is_weak_order, SVM_Benchmarking_Regression)

## only two of the artifical data sets yield weak orders
names(weak_orders)

## visualize them using Hasse diagrams
if(require("Rgraphviz")) plot(weak_orders)

## Consensus solutions:

data("SVM_Benchmarking_Classification_Consensus")
data("SVM_Benchmarking_Regression_Consensus")

## The solutions for the three families are not unique
print(SVM_Benchmarking_Classification_Consensus)
print(SVM_Benchmarking_Regression_Consensus)

## visualize the consensus weak orders
classW <- SVM_Benchmarking_Classification_Consensus$W
regrW <- SVM_Benchmarking_Regression_Consensus$W
if(require("Rgraphviz")) {
    plot(classW)
    plot(regrW)
}

## in tabular style:
ranking <- function(x) rev(names(sort(relation_class_ids(x))))
sapply(classW, ranking)
sapply(regrW, ranking)

## (prettier and more informative:)
relation_classes(classW[[1L]])
}
\keyword{datasets}
