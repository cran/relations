\name{algebra}
\alias{algebra}
\alias{relation_union}
\alias{relation_intersection}
\alias{relation_symdiff}
\alias{relation_complement}
\alias{relation_cartesian}
\alias{relation_projection}
\alias{relation_selection}
\alias{relation_division}
\alias{relation_join}
\alias{relation_semijoin}
\alias{relation_antijoin}
\alias{relation_remainder}
\alias{\%><\%}
\alias{\%=><\%}
\alias{\%><=\%}
\alias{\%=><=\%}
\alias{\%|><\%}
\alias{\%><|\%}
\alias{\%|><|\%}
\alias{\%|>\%}
\alias{\%<|\%}
\alias{\%U\%}
\title{Relational Algebra}
\description{Various \dQuote{relational algebra}-like operations.}
\usage{
relation_projection(x, margin = NULL)
relation_selection(x, subset)
relation_cartesian(x, y, \dots)
relation_union(x, y, \dots)
relation_complement(x, y)
relation_intersection(x, y, \dots)
relation_symdiff(x, y)
relation_division(x, y)
relation_remainder(x, y)
relation_join(x, y, \dots)
relation_semijoin(x, y, \dots)
relation_antijoin(x, y, \dots)
}
\arguments{
  \item{x, y}{Relation objects.}
  \item{margin}{Either a character vector of domain names, or an integer
    vector of domain indices.}
  \item{subset}{Expression resulting in a logical vector of length equal
    to the number of tuples in the graph.}
  \item{\dots}{Relation objects for \code{relation_union},
    \code{relation_intersection} and \code{relation_cartesian};
    otherwise, parameters passed to \code{\link[base]{merge}}.}
}
\details{
  These functions provide functionality similar to the corresponding
  operations defined in relational algebra theory as introduced by Codd
  (1970).  Note, however, that domains in database relations, unlike the
  concept of relations we use here, are unordered.  In fact, a database
  relation (\dQuote{table}) is defined as a set of elements called
  \dQuote{tuples}, where the \dQuote{tuple} componentes are named, but
  unordered.  So in fact, a \dQuote{tuple} in this sense is a set of
  mappings from the attribute names into the union of the attribute
  domains.

  The \emph{projection} of a relation on a specified margin (i.e., a
  vector of domain names or indices) is the relation obtained when all
  tuples are restricted to this margin. As a consequence, duplicate
  tuples are removed.

  The \emph{selection} of a relation is the relation obtained by taking
  a subset of the relation graph, defined by some logical expression.

  The \emph{Cartesian product} of two relations is obtained by basically
  building the Cartesian product of all graph elements, but combining the
  resulting pairs into single tuples.

  The \emph{union} of two relations simply combines the graph elements
  of both relations; the \emph{complement} of two relations \eqn{X} and
  \eqn{Y} removes the tuples of \eqn{Y} from \eqn{X}.

  The \emph{intersection} (\emph{symmetric difference})
  of two relations is the relation with all tuples they have (do not have)
  in common.

  The \emph{division} of relation \eqn{X} by relation \eqn{Y} is the
  reversed Cartesian product. The result is a relation with the domain
  unique to \eqn{X} and containing the maximum number of tuples which,
  multiplied by \eqn{Y}, are contained in \eqn{X}. The \emph{remainder}
  of this operation is the complement of \eqn{X} and the division of
  \eqn{X} by \eqn{Y}. Note that for both operations, the domain of
  \eqn{Y} must be contained in the domain of \eqn{X}.

  The (natural) \emph{join} of two relations is their Cartesian product,
  restricted to the subset where the elements of the common attributes
  do match. The left/right/full outer join of two relations
  \eqn{X} and \eqn{Y} is
  the union of \eqn{X}/\eqn{Y}/\eqn{X} and \eqn{Y},
  and the inner join of \eqn{X} and \eqn{Y}.
  The implementation uses \code{\link{merge}}, and so the
  left/right/full outer joins are obtained by setting
  \code{all.x}/\code{all.y}/\code{all} to \code{TRUE} in
  \code{relation_join}.
  The domains to be matched are specified using \code{by}.

  The left (right) \emph{semijoin} of two relations \eqn{X} and \eqn{Y}
  is the join of these, projected to the attributes of \eqn{X}
  (\eqn{Y}). Thus, it yields all tuples of \eqn{X}
  (\eqn{Y}) participating in the join of \eqn{X} and \eqn{Y}. 

  The left (right) \emph{antijoin} of two relations \eqn{X} and \eqn{Y}
  is the complement of \eqn{X} (\eqn{Y}) and the join of both,
  projected to the attributes of \eqn{X} (\eqn{Y}).
  Thus, it yields all tuples of \eqn{X}
  (\eqn{Y}) \emph{not} participating in the join of \eqn{X} and \eqn{Y}. 

  The operators \code{\%><\%}, \code{\%=><\%}, \code{\%><=\%},
  \code{\%=><=\%}, \code{\%|><\%}, \code{\%><|\%}, \code{\%|><|\%},
  \code{\%|>\%}, \code{\%<|\%}, and \code{\%U\%} can be used for the
  Cartesian product, left outer join, right outer join, full outer join,
  left semi-join, right semi-join, join, left antijoin, right antijoin,
  and union, respectively.
  
}
\references{
  E. F. Codd (1970),
  A relational model of data for large shared data banks.
  \emph{Communications of the ACM}, \bold{13}/6, 377--387. 
}
\seealso{
  \code{\link{relation}}
}
\examples{
## projection
Person <-
    data.frame(Name = c("Harry", "Sally", "George", "Helena", "Peter"),
               Age = c(34, 28, 29, 54, 34),
               Weight = c(80, 64, 70, 54, 80),
               stringsAsFactors = FALSE)
Person <- as.relation(Person)
relation_table(Person)
relation_table(relation_projection(Person, c("Age", "Weight")))

## selection
relation_table(R1 <- relation_selection(Person, Age < 29))
relation_table(R2 <- relation_selection(Person, Age >= 34))
relation_table(R3 <- relation_selection(Person, Age == Weight))

## union
relation_table(R1 \%U\% R2)

## works only for the same domains:
relation_table(R2 | R3)

## complement
relation_table(Person - R2)

## intersection
relation_table(relation_intersection(R2, R3))

## works only for the same domains:
relation_table(R2 & R3)

## symmetric difference
relation_table(relation_symdiff(R2, R3))

## Cartesian product
Employee <-
    data.frame(Name = c("Harry", "Sally", "George", "Harriet", "John"),
               EmpId = c(3415, 2241, 3401, 2202, 3999),
               DeptName = c("Finance", "Sales", "Finance", "Sales", "N.N."),
	       stringsAsFactors = FALSE)
Employee <- as.relation(Employee)
relation_table(Employee)
Dept <- data.frame(DeptName = c("Finance", "Sales", "Production"),
                   Manager = c("George", "Harriet", "Charles"),
                   stringsAsFactors = FALSE)
Dept <- as.relation(Dept)
relation_table(Dept)

relation_table(Employee \%><\% Dept)

## Natural join
relation_table(Employee \%|><|\% Dept)

## left (outer) join
relation_table(Employee \%=><\% Dept)

## right (outer) join
relation_table(Employee \%><=\% Dept)

## full outer join
relation_table(Employee \%=><=\% Dept)

## antijoin
relation_table(Employee \%|>\% Dept)
relation_table(Employee \%<|\% Dept)

## semijoin
relation_table(Employee \%|><\% Dept)
relation_table(Employee \%><|\% Dept)

## division
Completed <-
    data.frame(Student = c("Fred", "Fred", "Fred", "Eugene",
                           "Eugene", "Sara", "Sara"),
               Task = c("Database1", "Database2", "Compiler1",
                        "Database1", "Compiler1", "Database1",
                        "Database2"),
               stringsAsFactors = FALSE)
Completed <- as.relation(Completed)
relation_table(Completed)
DBProject <- data.frame(Task = c("Database1", "Database2"),
                        stringsAsFactors = FALSE)
DBProject <- as.relation(DBProject)
relation_table(DBProject)

relation_table(Completed \%/\% DBProject)

## division remainder
relation_table(Completed \%\% DBProject)

}
\keyword{math}
