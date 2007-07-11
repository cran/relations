\name{SVMBench}
\alias{SVM_Benchmarking_Classification}
\alias{SVM_Benchmarking_Regression}
\alias{SVM_Benchmarking_Classification_Consensus}
\alias{SVM_Benchmarking_Regression_Consensus}
\title{SVM Benchmarking Data and Consensus Relations}
\description{
  \code{SVM_Benchmarking_Classification} and
  \code{SVM_Benchmarking_Regression} represent the
  results of a benchmark study comparing Support Vector Machines to
  other predictive methods on real and artificial data sets
  involving classification and regression
  methods, respectively. \code{SVM_Benchmarking_Classification_Consensus}
  and \code{SVM_Benchmarking_Regression_Consensus}
  are consensus rankings derived from these data.
}
\usage{
data("SVM_Benchmarking_Classification")
data("SVM_Benchmarking_Regression")
data("SVM_Benchmarking_Classification_Consensus")
data("SVM_Benchmarking_Regression_Consensus")
}
\format{
  \code{SVM_Benchmarking_Classification}
  (\code{SVM_Benchmarking_Regression}) is an ensemble of 21
  (12) relations representing pairwise comparisons of 17 classification
  (10 regression) methods on 21 (12) data sets. Each relation of the
  ensemble summarizes the results for a particular data set. The
  relations are endorelations on the set of methods employed. Since some
  methods did fail on some data sets, the relations are not
  guaranteed to be complete, transitive, or even reflexive. See
  Meyer et al. (2003) for details on the experimental design of the benchmark
  study, and Hornik et Meyer (2007) for the pairwise comparisons.

  \code{SVM_Benchmarking_Classification_Consensus}
  and \code{SVM_Benchmarking_Regression_Consensus} are lists of
  ensembles of consensus relations fitted to the benchmark
  results. For each of the following three endorelation families:
  \code{SD/L} (\dQuote{linear orders}),
  \code{SD/O} (\dQuote{partial orders}), and \code{SD/P}
  (\dQuote{preferences}), \emph{all} possible consensus
  relations have been computed (see \code{\link{relation_consensus}}).
  For both classification and regression,
  the three relation ensembles obtained are provided as a named list of
  length 3. See Hornik et Meyer (2007) for details on the meta-analysis.
}
\source{
D. Meyer, F. Leisch, and K. Hornik (2003),
The support vector machine under test.
\emph{Neurocomputing}, \bold{55}:169--186.

K. Hornik and D. Meyer (2007),
Deriving consensus rankings from benchmarking experiments.
In R. Decker and H.-J. Lenz,
\emph{Advances in Data Analysis}.
Studies in Classification, Data Analysis, and Knowledge Organization.
Springer-Verlag.
}
\examples{
data("SVM_Benchmarking_Classification")

## 21 data sets
names(SVM_Benchmarking_Classification)

## 17 methods
relation_domain(SVM_Benchmarking_Classification)

## select preferences
P <- sapply(SVM_Benchmarking_Classification, relation_is_preference)

## only the artifical data sets yield preferences
names(SVM_Benchmarking_Classification)[P]

## visualize them using Hasse diagrams
if (require("Rgraphviz"))
    plot(SVM_Benchmarking_Classification[P])

## Same for regression:
data("SVM_Benchmarking_Regression")

## 12 data sets
names(SVM_Benchmarking_Regression)

## 10 methods
relation_domain(SVM_Benchmarking_Regression)

## select preferences
P <- sapply(SVM_Benchmarking_Regression, relation_is_preference)

## only two of the artifical data sets yield preferences
names(SVM_Benchmarking_Regression)[P]

## visualize them using Hasse diagrams
if (require("Rgraphviz"))
    plot(SVM_Benchmarking_Regression[P])

## Consensus solutions:

data("SVM_Benchmarking_Classification_Consensus")
data("SVM_Benchmarking_Regression_Consensus")

## The solutions for the three families are not unique
print(SVM_Benchmarking_Classification_Consensus)
print(SVM_Benchmarking_Regression_Consensus)

## visualize the consensus preferences
if (require("Rgraphviz")) {
    plot(SVM_Benchmarking_Classification_Consensus$P)
    plot(SVM_Benchmarking_Regression_Consensus$P)
}

## in tabular style:
ranking <- function(x) rev(names(sort(relation_class_ids(x))))
sapply(SVM_Benchmarking_Classification_Consensus$P, ranking)
sapply(SVM_Benchmarking_Regression_Consensus$P, ranking)

## (prettier and more informative:)
relation_classes(SVM_Benchmarking_Classification_Consensus$P[[1]])
}
\keyword{datasets}
