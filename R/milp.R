## A simple framework for representing and solving mixed integer linear
## programs (MILPs) of the form
##   optimize obj' x
##   such that mat %*% x dir rhs
## with possibly given types (C/I/B for continuous/integer/binary) and
## given additional (lower and upper) bounds on x.
## (Default of course x >= 0).

MILP <-
function(objective, constraints, bounds = NULL, types = NULL,
         maximum = FALSE)
{
    ## In the simples case, 'constraints' is a (not necessarily named)
    ## list with mat, dir and rhs.  Advanced solvers might allow for
    ## more advanced constraints, but let's worry about this later (and
    ## maybe also a little MILP_constraints() wrapper ...).

    structure(list(objective = objective, constraints = constraints,
                   bounds = bounds, types = types, maximum = maximum),
              class = "MILP")
}

solve_MILP <-
function(x, solver = c("lpsolve", "glpk", "symphony", "cplex"),
         control = list())
{
    ## Ideally, we would use some registration mechanism for solvers.
    
    ## One of the key ideas is to use allow using a simple triplet
    ## matrix representation for the constraint matrix, so that we can
    ## eventually take advantage of advanced solvers which can handle
    ## sparse constraint matrices.  For the time being:

    solver <- match.arg(solver)
    status <- NULL

    if(solver == "lpsolve") {
        ## Currently, no direct support for bounds.
        ## <FIXME>
        ## Should rewrite the given bounds into additional constraints.
        if(!is.null(x$bounds))
            stop("Solver currently does not support variable bounds.")
        ## </FIXME>
        ## Version 5.6.1 has added sparse matrix support via formal
        ## 'dense.const' as well as binary variable types.
        mat <- x$constraints[[1L]]
        is_sparse <- inherits(mat, "simple_triplet_matrix")
        types <- .expand_types(x$types, length(x$objective))
        out <- if(is_sparse)
            lpSolve::lp(if(x$maximum) "max" else "min",
                        x$objective,
                        const.dir = x$constraints[[2L]],
                        const.rhs = x$constraints[[3L]],
                        int.vec = which(types == "I"),
                        binary.vec = which(types == "B"),
                        dense.const = cbind(mat$i, mat$j, mat$v))
        else 
            lpSolve::lp(if(x$maximum) "max" else "min",
                        x$objective,
                        as.matrix(mat),
                        x$constraints[[2L]],
                        x$constraints[[3L]],
                        int.vec = which(types == "I"),
                        binary.vec = which(types == "B"))
        solution <- out$solution
        objval <- sum(solution * out$objective)
        status <- out$status
    }
    else if(solver == "glpk") {
        out <- Rglpk::Rglpk_solve_LP(x$objective,
                                     x$constraints[[1L]],
                                     x$constraints[[2L]],
                                     x$constraints[[3L]],
                                     bounds = x$bounds,
                                     types = x$types,
                                     max = x$maximum)
        solution <- out$solution
        objval <- out$optimum
        status <- out$status
    }
    else if(solver == "symphony") {
        out <- Rsymphony::Rsymphony_solve_LP(x$objective,
                                             x$constraints[[1L]],
                                             x$constraints[[2L]],
                                             x$constraints[[3L]],
                                             bounds = x$bounds,
                                             types = x$types,
                                             max = x$maximum)
        solution <- out$solution
        objval <- out$objval
        status <- out$status
    }
    else if(solver == "cplex") {
        ## Currently, no direct support for bounds.
        ## <FIXME>
        ## Should expand the given bounds and map into lb/ub arguments.
        if(!is.null(x$bounds))
            stop("Solver currently does not support variable bounds.")
        ## </FIXME>
        .as_Rcplex_sense <- function(x) {
            TABLE <- c("L", "L", "G", "G", "E")
            names(TABLE) <- c("<", "<=", ">", ">=", "==")
            TABLE[x]
        }
        sense <- .as_Rcplex_sense(x$constraints[[2L]])
        types <- .expand_types(x$types, length(x$objective))        
        mat <- x$constraints[[1L]]
        is_sparse <- inherits(mat, "simple_triplet_matrix")
        out <- if(is_sparse){
            ## reorder indices as CPLEX needs a column major order
            ## representation i.e., column indices j have to be in
            ## ascending order.
            column_major_order <- order(mat$j)
            mat$i <- mat$i[column_major_order]
            mat$j <- mat$j[column_major_order]
            mat$v <- mat$v[column_major_order]
            Rcplex::Rcplex(cvec = x$objective,
                           Amat = mat,
                           sense = sense,
                           bvec = x$constraints[[3L]],
                           vtype = types,
                           objsense = if(x$maximum) "max" else "min",
                           control = list(trace = 0, round = 1)
                           )
          }
        else
            Rcplex::Rcplex(cvec = x$objective,
                           Amat = as.matrix(mat),
                           sense = sense,
                           bvec = x$constraints[[3L]],
                           vtype = types,
                           objsense = if(x$maximum) "max" else "min",
                           control = list(trace = 0, round = 1)
                           )
        solution <- out$xopt
        ## For the time being ...
        ## since Rcplex 0.1-4 integers are rounded (via control argument)
        ## but no new optimal solution based on this values is calculated
        ## solution[(types == "I") | (types == "B")] <-
        ##  round(solution[(types == "I") | (types == "B")])
        objval <- sum(solution * x$objective)
        ## CPLEX returns 101 if the solution is optimal for MILPs
        status <- ifelse(out$status == 101, 0, out$status) 
    }
    list(solution = solution, objval = objval, status = status)
}

.make_types <-
function(n, I = NULL, B = NULL)
{
    ## Create MILP variable types spec from possibly given positions of
    ## integer and binary variables.
    types <- rep.int("C", n)
    if(!is.null(I)) types[I] <- "I"
    if(!is.null(B)) types[B] <- "B"
    types
}

.expand_types <-
function(x, n)
{
    if(is.null(x)) {
        ## Continuous by default.
        rep.int("C", n)
    }
    else {
        if(!is.character(x) || !all(x %in% c("C", "I", "B")))
            stop("Invalid MILP variable types.")
        ## Be nicer than necessary ...
        rep(x, length.out = n)
    }
}

