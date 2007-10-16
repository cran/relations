## A simple framework for representing and solving mixed integer linear
## programs (MILPs) of the form
##   optimize obj' x
##   such that mat %*% x dir rhs
##             x[int] integer
## (and of course x >= 0).

MILP <-
function(direction, objective, constraints, integers = integer())
{
    ## For the time being, no defaults for direction (min/max).
    direction <- match.arg(direction, c("min", "max"))

    ## In the simples case, 'constraints' is a (not necessarily named)
    ## list with mat, dir and rhs.  Advanced solvers might allow for
    ## more advanced constraints, but let's worry about this later (and
    ## maybe also a little MILP_constraints() wrapper ...).

    structure(list(direction = direction, objective = objective,
                   constraints = constraints, integers = integers),
              class = "MILP")
}

solve_MILP <- function(x, solver = "lp", control = list())
{
    ## Currently, lpSolve::lp is the only available/used solver, but
    ## hopefully we can have glpk and Symphony soon.

    ## One of the key ideas is to use allow using a simple triplet
    ## matrix representation for the constraint matrix, so that we can
    ## eventually take advantage of advanced solvers which can handle
    ## sparse constraint matrices.  For the time being:

    y <- lpSolve::lp(x$direction,
                     x$objective,
                     as.matrix(x$constraints[[1L]]),
                     x$constraints[[2L]],
                     x$constraints[[3L]],
                     int.vec = x$integers)
    z <- y$solution
    list(solution = z, value = sum(z * y$objective), status = y$status)
}
