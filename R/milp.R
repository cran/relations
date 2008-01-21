## A simple framework for representing and solving mixed integer linear
## programs (MILPs) of the form
##   optimize obj' x
##   such that mat %*% x dir rhs
##             x[int] integer
## (and of course x >= 0).

MILP <-
function(objective, constraints, integers = integer(), maximum = FALSE)
{
    ## In the simples case, 'constraints' is a (not necessarily named)
    ## list with mat, dir and rhs.  Advanced solvers might allow for
    ## more advanced constraints, but let's worry about this later (and
    ## maybe also a little MILP_constraints() wrapper ...).

    structure(list(objective = objective, constraints = constraints,
                   integers = integers, maximum = maximum),
              class = "MILP")
}

solve_MILP <-
function(x, solver = c("lpsolve", "glpk", "symphony"), control = list())
{
    ## Ideally, we would use some registration mechanism for solvers.
    
    ## One of the key ideas is to use allow using a simple triplet
    ## matrix representation for the constraint matrix, so that we can
    ## eventually take advantage of advanced solvers which can handle
    ## sparse constraint matrices.  For the time being:

    solver <- match.arg(solver)
    status <- NULL

    ## Currently, all solvers interfaced allow to specify the integer
    ## variables as a numeric vector with their *positions*.  Let us
    ## allow both kinds for the time being: if we interface solvers
    ## which only allow for a logical vectors of integer indicators,
    ## these need something like
    ##   if(is.numeric(integers <- x$integers)) {
    ##       integers <- logical(length(x$objective))
    ##       integers[x$integers] <- TRUE
    ##   }

    if(is.logical(integers <- x$integers))
        integers <- which(integers)

    if(solver == "lpsolve") {
        out <- lpSolve::lp(if(x$maximum) "max" else "min",
                           x$objective,
                           as.matrix(x$constraints[[1L]]),
                           x$constraints[[2L]],
                           x$constraints[[3L]],
                           int.vec = integers)
        solution <- out$solution
        objval <- sum(solution * out$objective)
        status <- out$status
    }
    else if(solver == "glpk") {
        out <- Rglpk::Rglpl_solve_LP(x$objective,
                                     x$constraints[[1L]],
                                     x$constraints[[2L]],
                                     x$constraints[[3L]],
                                     x$integers,
                                     x$maximum)
        solution <- out$solution
        objval <- out$optimum
        status <- out$status

    }
    else if(solver == "symphony") {
        out <- Rsymphony::Rsymphony_solve_LP(x$objective,
                                             x$constraints[[1L]],
                                             x$constraints[[2L]],
                                             x$constraints[[3L]],
                                             x$integers,
                                             x$maximum)
        solution <- out$solution
        objval <- out$optimum
        status <- out$status
    }
    list(solution = solution, objval = objval, status = status)
}
