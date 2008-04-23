## Relation fitters.

### * fit_relation_symdiff

## General purpose fitter for families
.SD_families <- c("E", "L", "O", "P", "T", "C", "A", "S", "R")
## where:
##   E ... Equivalence relations               REF SYM TRA
##   L ... Linear orders                   TOT     ASY TRA
##   O ... Partial order                           ASY TRA
##   P ... complete Preorder ("ordering")  TOT REF     TRA
##   T ... Tournament                      TOT IRR ASY
##   C ... Complete relations              TOT
##   A ... Antisymmetric relations                 ASY
##   S ... Symmetric relations                     SYM
##   M ... Matches                         TOT REF
##   R ... tRansitive relations                        TRA
## and
##   TOT ... total/complete
##   REF ... reflexive
##   ASY ... antisymmetric,
##   IRR ... irreflexive
##   TRA ... transitive
##
## Families which are not REF or IRR by definition will be reflexive iff
## the (weighted) majority of all input relations is.

## Keep this is sync with .relation_meta_db (currently in utilities.R).

## <NOTE>
## Ideally we would move to a database for families which specify their
## parametrizations (currently, upper triangular or off-diagonal) and
## constraints etc., in particular making it more convenient and
## reliable to add support for new families.
## A single classification (e.g., upper_tri vs. offdiag) does not do the
## job.
## </NOTE>

## Number of non-redundant object pairs.
.n_of_pairs <-
function(n, family)
{
    family <- match.arg(family, .SD_families)
    N <- n * (n - 1)                    # Number of distinct pairs.
    switch(EXPR = family,
           E =, L =, S =, T = N / 2,    # upper_tri parametrization
           A =, C =, M =, O =, P =, R = N    # offdiag parametrization
           )
}

## Number of transitivity constraints for the non-redundant pairs.
## (None for tournaments.)
.n_of_transitivity_constraints <-
function(n, family)
{
    family <- match.arg(family, .SD_families)
    N <- n * (n - 1) * (n - 2)          # Number of distinct triples.
    switch(EXPR = family,
           E = N / 2, L = N / 3, O =, P =, R = N,
           A =, C =, M =, S =, T = 0L   # Families w/out TRA.
           )
}

## Make function giving the position of incidence (i, j) in the vector
## of non-redundant incidences used for symdiff fitting (upper.tri() for
## E/L/S/T and .offdiag() for A/C/M/O/P/R, respectively).
.make_pos <-
function(n, family)
{
    family <- match.arg(family, .SD_families)
    if(family %in% c("E", "L", "S", "T"))
        function(i, j) {
            ## Position of x_{ij} in x[upper.tri(x)].
            i + (j - 1) * (j - 2) / 2
        }
    else
        function(i, j) {
            ## Position of x_{ij} in x[.offdiag(x)]
            (n - 1) * (j - 1) + i - (i >= j)
        }
}

fit_relation_symdiff <-
function(x, family = .SD_families, control = list())
{
    sparse <- control$sparse
    if(is.null(sparse)) sparse <- FALSE

    ## Number of objects:
    n <- nrow(x)
    objective_in <- if (family %in% c("L", "T"))
        (x - t(x))[upper.tri(x)]
    else if (family %in% c("E", "S"))
        (x + t(x))[upper.tri(x)]
    else
        x[.offdiag(x)]

    ## Handle constraints implied by the family to be fitted.
    ## Need all variables in { 0 , 1 }, i.e., >= 0 and <= 1 and
    ## integer, plus maybe totality or antisymmetry, plus maybe
    ## transitivity.
    NP <- .n_of_pairs(n, family)
    ## <NOTE>
    ## At least for the time being, these families are exactly the
    ## families using an off-diagonal parametrization (except "R").
    family_is_A_or_C_or_M_or_O_or_P <-
        family %in% c("A", "C", "M", "O", "P")
    ## </NOTE>
    eye <-
        if(!sparse) diag(1, NP) else .simple_triplet_diag_matrix(1, NP)
    constr_mat <-
        rbind(eye,
              eye,
              if(family_is_A_or_C_or_M_or_O_or_P)
              .make_tot_or_asy_constraint_mat(n, sparse),
              .make_transitivity_constraint_mat(n, family, sparse))
    constr_dir <-
        c(rep.int(">=", NP),
          rep.int("<=", NP),
          if(family_is_A_or_C_or_M_or_O_or_P)
              .make_tot_or_asy_constraint_dir(n, family),
          .make_transitivity_constraint_dir(n, family))
    constr_rhs <-
        c(rep.int(0, NP),
          rep.int(1, NP),
          if(family_is_A_or_C_or_M_or_O_or_P)
              .make_tot_or_asy_constraint_rhs(n),
          .make_transitivity_constraint_rhs(n, family))

    ## Handle additional constraints.
    acmaker <-
        .make_additional_constraint_maker_using_incidences(n, family, sparse)
    if(!is.null(A <- control$constraints)) {
        A <- .canonicalize_additional_constraints(A, family)
        add <- acmaker(A)
        constr_mat <- rbind(constr_mat, add$mat)
        constr_dir <- c(constr_dir, add$dir)
        constr_rhs <- c(constr_rhs, add$rhs)
    }

    labels <- dimnames(x)

    ## Compute diagonal:
    ## For families which are reflexive or irreflexive, set all ones or
    ## zeroes, respectively.  For all other families, a diagonal entry
    ## will be set iff it is in the (non-strict weighted) majority of
    ## the input relations.
    ## <FIXME>
    ## What is a reasonable thing to do when using the fitter code to
    ## determine fuzzy euclidean consensus relations?
    ## </FIXME>
    diagonal <- if(family %in% c("E", "M", "O"))
        rep.int(1, n)
    else if(family == "T")
        rep.int(0, n)
    else
        diag(x) >= 0

    if(!is.null(all <- control$all) && identical(all, TRUE)) {
        verbose <- control$verbose
        if(is.null(verbose))
            verbose <- getOption("verbose")
        return(.find_all_relation_symdiff_optima(n, family,
                                                 objective_in,
                                                 constr_mat,
                                                 constr_dir,
                                                 constr_rhs,
                                                 seq_len(NP),
                                                 A,
                                                 acmaker,
                                                 labels,
                                                 verbose,
                                                 diagonal,
                                                 control$solver,
                                                 control$control))
    }

    out <- solve_MILP(MILP(objective_in,
                           list(constr_mat,
                                constr_dir,
                                constr_rhs),
                           seq_len(NP),
                           maximum = TRUE),
                      control$solver,
                      control$control)
    .stop_if_lp_status_is_nonzero(out$status, family)

    ## Turn the solution back into a full incidence matrix.
    fit <- if(family_is_A_or_C_or_M_or_O_or_P || (family == "R"))
        .make_incidence_from_offdiag(round(out$solution),
                                     family, labels, diagonal)
    else
        .make_incidence_from_upper_tri(round(out$solution),
                                       family, labels, diagonal)
    ## For the time being, tack some of the MILP results on so that we
    ## can look at them (but not everything due to size ...)
    attr(fit, ".lp") <- out[c("solution", "objval")]

    fit
}

### * .stop_if_lp_status_is_nonzero

.stop_if_lp_status_is_nonzero <-
function(status, family)
{
    if(status != 0) {
        ## This should really only be possible in case additional
        ## constraints were given.
        stop(gettextf("Given constraints are incompatible with family '%s'.",
                      family))
    }
}

### * .find_all_relation_symdiff_optima

.find_all_relation_symdiff_optima <-
function(n, family, obj, mat, dir, rhs, int, A, acmaker, labels,
         verbose = FALSE, diagonal, solver, control)
{
    ## Start by computing one optimal solution.
    out <- solve_MILP(MILP(obj, list(mat, dir, rhs), int,
                           maximum = TRUE),
                       solver, control)
    ## Check status:
    .stop_if_lp_status_is_nonzero(out$status, family)
    ## Value of the optimum:
    Vopt <- out$objval

    ## And now find all optimal additional constrains (i.e., at the end,
    ## all optimal incidences in 3-column matrix (i, j, x_{ij}) form) by
    ## looping over all non-redundant (i, j) pairs.
    ##
    ## Somewhat tricky: in case the fitting is to include additional
    ## constraints, we need to ensure that we only try to add feasible
    ## additional constraints.

    ## Determine a 2-column matrix of (i, j) index pairs to loop over.
    ind <- if(family %in% c("E", "L", "S", "T")) {
        ## Loop over all pairs 1 <= i < j <= n.
        which(upper.tri(matrix(0, n, n)), arr.ind = TRUE)
    }
    else {
        ## Loop over all pairs 1 <= i != j <= n.
        which(.offdiag(matrix(0, n, n)), arr.ind = TRUE)
    }
    ## But exclude pairs already specified in additional constraints.
    if(!is.null(A)) {
        hpos <- if(family %in% c("E", "L", "S", "T")) {
            function(i, j)
                ifelse(i < j, n * (j - 1) + i, n * (i - 1) + j)
        }
        else {
            function(i, j)
                n * (j - 1) + i
        }
        ind <- ind[apply(outer(hpos(ind[, 1L], ind[, 2L]),
                               hpos(A[, 1L], A[, 2L]),
                               "!="),
                         2, all), ]
    }

    Alist <- list(matrix(0, 0L, 3L))
    for(k in seq_len(nrow(ind))) {
        Alist <- do.call("c",
                         lapply(Alist,
                                .add_constraint_for_single_pair,
                                ind[k, 1L], ind[k, 2L],
                                acmaker, Vopt, obj, mat, dir, rhs, int,
                                solver, control))
        if(verbose)
            message(gettextf("N_of_pairs: %d *** N_of_optimal_branches: %d",
                             k, length(Alist)))
    }

    lapply(Alist, .make_incidence_from_triples, family, labels, diagonal)
}

.add_constraint_for_single_pair <-
function(A, i, j, acmaker, Vopt, obj, mat, dir, rhs, int,
         solver, control)
{
    ## Figure out whether additional constraint x_{ij} = 0 or x_{ij} =
    ## 1, or both, are "optimal".

    optimize_with_valid_additional_constraints <-
        function(A) {
        ## Compute optimum for a valid matrix A of additional
        ## constraints.
        add <- acmaker(A)
        out <- solve_MILP(MILP(obj,
                               list(rbind(mat, add$mat),
                                    c(dir, add$dir),
                                    c(rhs, add$rhs)),
                               int,
                               maximum = TRUE),
                          solver, control)
        if(out$status != 0) return(-Inf)
        out$objval
    }

    Aij0 <- rbind(A, c(i, j, 0))
    Aij1 <- rbind(A, c(i, j, 1))
    ## <NOTE>
    ## At least when using the objval given by lpSolve:lp(), we cannot
    ## reliably test v < Vopt due to numerical precision issues.  This
    ## seems to be less of an issue when using
    ##   sum(out$objective * out$solution)
    ## (as done by solve_MILP()) but then we've still seen cases where
    ## max(v0, v1) < Vopt, hence we compare with some tolerance (could
    ## improve ...).
    tol <- 1e-10
    ## (Smaller than .Machine$double.eps^0.5 as e.g. used in the
    ## all.equal() comparisons.)
    v <- optimize_with_valid_additional_constraints(Aij0)
    if(v < Vopt - tol) return(list(Aij1))
    v <- optimize_with_valid_additional_constraints(Aij1)
    if(v < Vopt - tol) return(list(Aij0))
    ## To debug, use the following.
    ##   v0 <- optimize_with_valid_additional_constraints(Aij0)
    ##   v1 <- optimize_with_valid_additional_constraints(Aij1)
    ##   cat(i, j, v0, v1, Vopt, "\n")
    ##   if(all(c(v0, v1) < Vopt)) browser()
    ##   if(v0 < Vopt) return(list(Aij1))
    ##   if(v1 < Vopt) return(list(Aij0))
    ##   cat("splitting", i, j, "\n")
    list(Aij0, Aij1)
}

### * .find_all_relation_symdiff_optima_E_or_P_k

.find_all_relation_symdiff_optima_E_or_P_k <-
function(n, nc, family, obj, mat, dir, rhs, int, A, acmaker, labels,
         verbose = FALSE, solver, control)
{
    ## <NOTE>
    ## This still has a lot of code duplicated with
    ##   .find_all_relation_symdiff_optima()
    ## </NOTE>

    ## Start by computing one optimal solution.
    out <- solve_MILP(MILP(obj, list(mat, dir, rhs), int,
                           maximum = TRUE),
                       solver, control)
    ## Check status:
    .stop_if_lp_status_is_nonzero(out$status, family)
    ## Value of the optimum:
    Vopt <- out$objval

    ## Need to loop over all pairs of objects i and classes k:
    ind <- cbind(rep(seq_len(n), each = nc),
                 rep.int(seq_len(nc), n))

    ## (Note that we currently do not worry about pairs already
    ## specified in additional constraints: this is not possible
    ## directly anyways ...)

    Alist <- list(matrix(0, 0L, 3L))
    for(k in seq_len(nrow(ind))) {
        Alist <- do.call("c",
                         lapply(Alist,
                                .add_constraint_for_single_pair,
                                ind[k, 1L], ind[k, 2L],
                                acmaker, Vopt, obj, mat, dir, rhs, int,
                                solver, control))
        if(verbose)
            message(gettextf("N_of_pairs: %d *** N_of_optimal_branches: %d",
                             k, length(Alist)))
    }

    lapply(Alist, .make_incidence_from_class_membership_triples,
           family, labels)
}

.make_incidence_from_class_memberships <-
function(M, family, labels)
{
    family <- match.arg(family, c("E", "P"))
    I <- if(family == "E")
        tcrossprod(M)
    else {
        nc <- ncol(M)
        E <- matrix(0, nc, nc)
        E[row(E) <= col(E)] <- 1
        tcrossprod(M %*% E, M)
    }
    dimnames(I) <- labels
    I
}

.make_incidence_from_class_membership_triples <-
function(x, family, labels)
{
    M <- matrix(0, max(x[, 1L]), max(x[, 2L]))
    M[x[, -3L, drop = FALSE]] <- x[, 3L]
    ## Could also (more efficiently?) do:
    ##  ind <- which(x[, 3L] > 0)
    ##  M[x[ind, -3L, drop = FALSE]] <- 1
    .make_incidence_from_class_memberships(M, family, labels)
}

### * Constraint generators: transitivity.

### ** .make_transitivity_constraint_mat

.make_transitivity_constraint_mat <-
function(n, family, sparse = FALSE)
{
    family <- match.arg(family, .SD_families)

    NP <- .n_of_pairs(n, family)
    if ((n <= 2L) || (family %in% c("A", "C", "M", "S", "T"))) {
        if(!sparse)
            return(matrix(0, 0, NP))
        else
            return(.simple_triplet_zero_matrix(0, NP))
    }

    NC <- .n_of_transitivity_constraints(n, family)
    pos <- .make_pos(n, family)

    if(family %in% c("E", "L")) {
        ## Create a matrix with all combinations of triples i < j < k in
        ## the rows.
        ind <- seq_len(n)
        z <- as.matrix(expand.grid(ind, ind, ind))[, c(3L, 2L, 1L)]
        z <- z[(z[, 1L] < z[, 2L]) & (z[, 2L] < z[, 3L]), , drop = FALSE]

        p_ij <- pos(z[, 1L], z[, 2L])
        p_ik <- pos(z[, 1L], z[, 3L])
        p_jk <- pos(z[, 2L], z[, 3L])

        ind <- seq_len(NC)

        if(!sparse) {
            out <- matrix(0, NC, NP)
            if(family == "E") {
                ## For equivalence relations, we have 3 transitivity
                ## constraints for each such triple.
                out[cbind(ind, c(p_ij, p_ij, p_ik))] <- 1
                out[cbind(ind, c(p_jk, p_ik, p_jk))] <- 1
                out[cbind(ind, c(p_ik, p_jk, p_ij))] <- -1
            }
            else if(family == "L") {
                ## For linear orders, we only get two constraints.
                NT <- NC / 2
                out[cbind(ind, c(p_ij, p_ij))] <-
                    rep(c(1, -1), each = NT)
                out[cbind(ind, c(p_jk, p_jk))] <-
                    rep(c(1, -1), each = NT)
                out[cbind(ind, c(p_ik, p_ik))] <-
                    rep(c(-1, 1), each = NT)
            }
        } else {
            out <- if(family == "E")
                simple_triplet_matrix(rep.int(ind, 3L),
                                      c(p_ij, p_ij, p_ik,
                                        p_jk, p_ik, p_jk,
                                        p_ik, p_jk, p_ij),
                                      c(rep.int(1, 2 * NC),
                                        rep.int(-1, NC)),
                                      NC, NP)
            else if(family == "L") {
                NT <- NC / 2
                simple_triplet_matrix(rep.int(ind, 3L),
                                      c(p_ij, p_ij,
                                        p_jk, p_jk,
                                        p_ik, p_ik),
                                      c(rep(c(1, -1), each = NT),
                                        rep(c(1, -1), each = NT),
                                        rep(c(-1, 1), each = NT)),
                                      NC, NP)
            }
        }
    }
    else {
        ## Create a matrix with all combinations of distinct triples i,
        ## j, k in the rows.
        ind <- seq_len(n)
        z <- as.matrix(expand.grid(ind, ind, ind))[, c(3L, 2L, 1L)]
        z <- z[(z[, 1L] != z[, 2L])
               & (z[, 2L] != z[, 3L])
               & (z[, 3L] != z[, 1L]), ]

        ind <- seq_len(NC)
        if(!sparse) {
            out <- matrix(0, NC, NP)
            out[cbind(ind, pos(z[, 1L], z[, 2L]))] <- 1
            out[cbind(ind, pos(z[, 2L], z[, 3L]))] <- 1
            out[cbind(ind, pos(z[, 1L], z[, 3L]))] <- -1
        } else {
            out <- simple_triplet_matrix(rep.int(ind, 3L),
                                         c(pos(z[, 1L], z[, 2L]),
                                           pos(z[, 2L], z[, 3L]),
                                           pos(z[, 1L], z[, 3L])),
                                         c(rep.int(1, 2 * NC),
                                           rep.int(-1, NC)),
                                         NC, NP)
        }
    }

    out
}

### ** .make_transitivity_constraint_dir

.make_transitivity_constraint_dir <-
function(n, family)
{
    if(n <= 2) return(character())
    family <- match.arg(family, .SD_families)
    rep.int("<=", .n_of_transitivity_constraints(n, family))
}

### ** .make_transitivity_constraint_rhs

.make_transitivity_constraint_rhs <-
function(n, family)
{
    if(n <= 2) return(double())

    family <- match.arg(family, .SD_families)

    NC <- .n_of_transitivity_constraints(n, family)

    if(family == "L")
        rep(c(1, 0), each = NC / 2)
    else
        rep.int(1, NC)
}

### * Constraint generators: completeness or antisymmetry.

## This translates into
##
##    x_{ij} + x_{ji} >= 1      [completeness: P]
##    x_{ij} + x_{ji} <= 1      [antisymmetry: O]
##
## For tournaments, we have both, and hence x_{ij} + x_{ji} = 1 as for
## linear orders, such that we use the non-redundant upper diagonal
## pairs representation, and get no additional explicit constraints.
##
## As these constraints only differ by direction, we handle them
## together.

### ** .make_tot_or_asy_constraint_mat

.make_tot_or_asy_constraint_mat <-
function(n, sparse = FALSE)
{
    if(n <= 1) {
        if(!sparse)
            return(matrix(0, 0L, 0L))   # :-)
        else
            return(.simple_triplet_zero_matrix(0L, 0L))
    }

    ## Position of x_{ij} in x[.offdiag(x)]
    pos <- .make_pos(n, "P")

    ## Number of non-redundant pairs.
    NP <- n * (n - 1)
    ## Number of constraints.
    NC <- NP / 2

    ind <- seq_len(n)
    z <- as.matrix(expand.grid(ind, ind))[, c(2L, 1L)]
    z <- z[z[, 1L] < z[, 2L], , drop = FALSE]

    ind <- seq_len(NC)

    if(!sparse) {
        out <- matrix(0, NC, NP)
        out[cbind(ind, pos(z[, 1L], z[, 2L]))] <- 1
        out[cbind(ind, pos(z[, 2L], z[, 1L]))] <- 1
    }
    else {
        out <- simple_triplet_matrix(c(ind, ind),
                                     c(pos(z[, 1L], z[, 2L]),
                                       pos(z[, 2L], z[, 1L])),
                                     rep.int(1, 2 * NC),
                                     NC, NP)
    }

    out
}

### ** .make_tot_or_asy_constraint_dir

.make_tot_or_asy_constraint_dir <-
function(n, family)
{
    NC <- n * (n - 1) / 2
    rep.int(switch(EXPR = family,
                   A =, O = "<=",       # ASY
                   C =, M =, P = ">="   # TOT
                   ),
            NC)
}

### ** .make_tot_or_asy_constraint_rhs

.make_tot_or_asy_constraint_rhs <-
function(n)
{
    NC <- n * (n - 1) / 2
    rep.int(1, NC)
}

### * Constraint generators: explicitly given constraints.

### .make_additional_constraint_maker_using_incidences

.make_additional_constraint_maker_using_incidences <-
function(n, family, sparse = FALSE)
    function(A) {
        ## (This assumes A has already been validated and
        ## canonicalized.)
        na <- nrow(A)
        nc <- .n_of_pairs(n, family)
        pos <- .make_pos(n, family)
        if(!sparse) {
            mat <- matrix(0, na, nc)
            mat[cbind(seq_len(na), pos(A[, 1L], A[, 2L]))] <- 1
        } else {
            mat <- simple_triplet_matrix(seq_len(na),
                                         pos(A[, 1L], A[, 2L]),
                                         rep.int(1, na),
                                         na, nc)
        }
        list(mat = mat,
             dir = rep.int("==", na),
             rhs = A[, 3L])
    }

### .make_additional_constraint_maker_using_memberships_E

.make_additional_constraint_maker_using_memberships_E <-
function(pos, n_of_variables, nc, sparse = FALSE)
    function(A) {
        ## (This assumes A has already been validated and
        ## canonicalized.)
        na <- nrow(A)
        ## For equivalences (E):
        ## A constraint of the form (i, j, 1) implies that
        ## objects i and j are in the same class, i.e.:
        ##   m_{ik} - m_{jk} == 0   for all k.
        ## A constraint of the form (i, j, 0) implies that
        ## objects i and j are in different classes, i.e.:
        ##   m_{ik} + m_{jk} <= 1   for all k.
        ind <- seq_len(na)
        if(!sparse) {
            mat <- matrix(0, nc * na, n_of_variables)
            for(k in seq_len(nc)) {
                mat[cbind(ind, pos(A[, 1L], k))] <- 1
                mat[cbind(ind, pos(A[, 2L], k))] <-
                    ifelse(A[, 3L] == 1, -1, 1)
                ind <- ind + na
            }
        } else {
            i <- rep.int(ind, 2 * nc) +
                rep(seq(from = 0, by = na, length.out = nc),
                    each = 2 * na)
            j <- sapply(seq_len(nc),
                        function(k)
                        c(pos(A[, 1L], k), pos(A[, 2L], k)))
            v <- rep.int(c(rep.int(1, na),
                           ifelse(A[, 3L] == 1, -1, 1)),
                         nc)
            mat <- simple_triplet_matrix(i, j, v,
                                         nc * na, n_of_variables)
        }
        list(mat = mat,
             dir = rep.int(ifelse(A[, 3L] == 1, "==", "<="), nc),
             rhs = rep.int(ifelse(A[, 3L] == 1, 0, 1), nc))
    }

### .make_additional_constraint_maker_using_memberships_P

.make_additional_constraint_maker_using_memberships_P <-
function(pos, n_of_variables, nc, sparse = FALSE)
    function(A) {
        ## (This assumes A has already been validated and
        ## canonicalized.)
        na <- nrow(A)
        ## For preferences (P):
        ## A constraint of the form (i, j, 1) means that i <= j,
        ## or class(i) <= class(j).
        ## I.e., can only have m_{jk} = 1 if \sum_{l <= k} m_{il} = 1,
        ## giving
        ##   m_{jk} <= \sum_{l <= k} m_{il}   for all k.
        ## A constraint of the form (i, j, 0) means that i > j,
        ## or class(i) > class(j).
        ## I.e., can only have m_{ik} = 1 if \sum_{l < k} m_{jl} = 1,
        ## giving
        ##   m_{ik} <= \sum_{l < k} m_{jl}    for all k,
        ## (with the empty sum for k = 1 taken as zero).
        if(!sparse) {
            mat <- matrix(0, nc * na, n_of_variables)
            ind <- i1 <- which(A[, 3L] == 1)
            if(any(i1)) {
                for(k in seq_len(nc)) {
                    mat[cbind(ind, pos(A[i1, 2L], k))] <- 1
                    for(l in seq_len(k))
                        mat[cbind(ind, pos(A[i1, 1L], l))] <- -1
                    ind <- ind + na
                }
            }
            ind <- i0 <- which(A[, 3L] == 0)
            if(any(i0)) {
                for(k in seq_len(nc)) {
                    mat[cbind(ind, pos(A[i0, 1L], k))] <- 1
                    for(l in seq_len(k - 1L))
                        mat[cbind(ind, pos(A[i0, 2L], l))] <- -1
                    ind <- ind + na
                }
            }
        } else {
            i0 <- i1 <- j0 <- j1 <- integer()
            v0 <- v1 <- double()
            ## This is somewhat messy to compute explicitly in one pass,
            ## so let's simply use a loop for building things up (but
            ## avoid looping over j <- c(j, something) constructs for
            ## performance reasons.
            ind <- which(A[, 3L] == 0)
            if(any(ind)) {
                len <- length(ind)
                i0 <- unlist(lapply(seq_len(nc),
                                    function(k)
                                    rep.int(ind + (k - 1L) * na, k)
                                    ))
                j0 <- unlist(lapply(seq_len(nc),
                                    function(k)
                                    c(pos(A[i0, 1L], k),
                                      unlist(lapply(seq_len(k - 1L),
                                                    function(l)
                                                    pos(A[i0, 2L], l))))
                                    ))
                v0 <- unlist(lapply(seq_len(nc),
                                    function(k)
                                    c(rep.int(1, len),
                                      rep.int(-1, (k - 1L) * len))))
            }
            ind <- which(A[, 3L] == 1)
            if(any(ind)) {
                len <- length(ind)
                i1 <- unlist(lapply(seq_len(nc),
                                    function(k)
                                    rep.int(ind + (k - 1L) * na, k + 1L)
                                    ))
                j1 <- unlist(lapply(seq_len(nc),
                                    function(k)
                                    c(pos(A[ind, 2L], k),
                                      unlist(lapply(seq_len(k),
                                                    function(l)
                                                    pos(A[ind, 1L], l))))
                                    ))
                v1 <- unlist(lapply(seq_len(nc),
                                    function(k)
                                    c(rep.int(1, len),
                                      rep.int(-1, k * len))))
            }
            mat <- simple_triplet_matrix(c(i0, i1),
                                         c(j0, j1),
                                         c(v0, v1),
                                         nc * na, n_of_variables)
        }

        list(mat = mat,
             dir = rep.int("<=", nc * na),
             rhs = double(nc * na))
    }

### ** .make_additional_constraint_maker_using_o_c_pairs

.make_additional_constraint_maker_using_o_c_pairs <-
function(pos, n_of_variables, nc, sparse = FALSE)
    function(A) {
        ## (This assumes A has already been validated and
        ## canonicalized.)
        na <- nrow(A)
        if(!sparse) {
            mat <- matrix(0, na, n_of_variables)
            mat[cbind(seq_len(na), pos(A[, 1L], A[, 2L]))] <- 1
        } else {
            mat <- simple_triplet_matrix(seq_len(na),
                                         pos(A[, 1L], A[, 2L]),
                                         rep.int(1, na),
                                         na, n_of_variables)
        }
        list(mat = mat, dir = rep.int("==", na), rhs = A[, 3L])
    }

### ** .canonicalize_additional_constraints

.canonicalize_additional_constraints <-
function(A, family)
{
    ## Additional constraints should be given as a 3-column matrix with
    ## rows (i, j, x) meaning that the incidence of i and j should be
    ## equal to x.
    if(!is.matrix(A) || (ncol(A) != 3L))
        stop("Invalid additional incidence constraints.")
    ## For all families, we only use off-diagonal terms, so drop
    ## diagonal ones (which should all be one, of course).
    A <- A[A[, 1L] != A[, 2L], , drop = FALSE]
    if(nrow(A) == 0L) return(matrix(0, 0L, 3L))
    ## For E/L/S/T, we use only pairs i < j, so swap and possibly flip
    ## (L/T) if needed.
    if(family %in% c("E", "L", "S", "T")) {
        ind <- A[, 1L] > A[, 2L]
        if(any(ind))
            A[ind, ] <- cbind(A[ind, c(2L, 1L), drop = FALSE],
                              if(family %in% c("E", "S")) A[ind, 3L]
                              else 1 - A[ind, 3L])
    }
    ## Now validate.
    pos <- max(A) * (A[, 2L] - 1) + A[, 1L]
    ind <- duplicated(pos)
    if(any(ind)) {
        ## Check if the duplicated entries all have the same
        ## incidences.
        lens <- tapply(A[, 3L], pos, function(t) length(unique(t)))
        if(any(lens > 1))
            stop("Incompatible constraints.")
        ## Drop duplicated entries.
        A <- A[!ind, , drop = FALSE]
    }
    A
}

### * Incidence generators

### ** .make_incidence_from_upper_tri

.make_incidence_from_upper_tri <-
function(x, family = c("E", "L", "S", "T"),
         labels = NULL, diagonal)
{
    ## Compute the indicences of a binary relation from the given family
    ## from its upper triangular part (provided this is possible, of
    ## course).

    family <- match.arg(family)

    if(family %in% c("E", "S")) {
        ## Equivalence or symmetric relation.
        y <- diag(diagonal / 2)
        y[upper.tri(y)] <- x
        y <- y + t(y)
    }
    else {
        ## Linear order or tournament.
        y <- diag(diagonal)
        y[upper.tri(y)] <- x
        y[lower.tri(y)] <- 1 - t(y)[lower.tri(y)]
    }

    dimnames(y) <- labels
    y
}

### ** .make_incidence_from_offdiag

.make_incidence_from_offdiag <-
function(x, family = c("A", "C", "M", "O", "P", "R"),
         labels = NULL, diagonal)
{
    family <- match.arg(family)

    ## Compute the indicences of a binary relation from its off-diagonal
    ## part (provided this is possible, of course).

    ## <NOTE>
    ## Diagonal entries for the incidences are a mess.
    ## For antisymmetric relations, it really does/should not matter.
    ## We use the standard family definitions to infer reflexivity or
    ## irreflexivity; where neither is implied, we use a majority vote
    ## for reflexivity to determine the diagonal entries.
    ## </NOTE>

    y <- diag(diagonal)
    y[.offdiag(y)] <- x

    dimnames(y) <- labels
    y
}

### ** .make_incidence_from_triples

.make_incidence_from_triples <-
function(x, family, labels = NULL, diagonal)
{
    ## Compute the indicences of a binary relation from a 3-column
    ## matrix the rows of which are triples (i, j, x_{ij}) with x_{ij}
    ## the incidence at position (i, j).

    ## Set up incidences for the non-redundant pairs.
    I <- diag(diagonal)
    I[x[, -3L, drop = FALSE]] <- x[, 3L]
    ## And complete according to family.
    if(family %in% c("E", "L", "S", "T")) {
        ind <- lower.tri(I)
        I[ind] <- if(family %in% c("E", "S"))
            t(I)[ind]
        else
            1 - t(I)[ind]
    }

    dimnames(I) <- labels
    I
}

### * fit_relation_symdiff_E_k

fit_relation_symdiff_E_k <-
function(C, nc, control = list())
{
    ## Fit equivalence with at most nc classes.

    sparse <- control$sparse
    if(is.null(sparse)) sparse <- FALSE

    ## Using the membership matrix M = [m_{ik}] corresponding to the
    ## equivalence relation (partition), we have to maximize
    ##   \sum_{i,j,k} c_{ij} m_{ik} m_{jk}
    ## over all binary stochastic matrices M.  This is done by a simple
    ## linearization of the quadratic integer program.
    ##
    ## Let
    no <- nrow(C)
    ## be the number of objects.  Then the membership matrix M is a
    ## binary stochastic no x nc matrix.  There are no^2 x nc products
    ## y_{ijk} = m_{ik} m_{jk}.  For each we have the constraints
    ##    y_{ijk} - m_{ik} - m_{jk} >= -1,
    ##    y_{ijk} - m_{ik}          <=  0
    ##    y_{ijk}          - m_{jk} <=  0
    ## and each y_{ijk} gets coefficient c_{ij} in the objective
    ## function.

    ## We put all decision variables into a matrix [c(Y), c(M)] (so that
    ## the Y part has k "blocks" because k varies the slowest).

    n_of_y_variables <- no ^ 2 * nc
    n_of_m_variables <- no * nc
    n_of_variables <- n_of_y_variables + n_of_m_variables

    pos_m <- function(i, k) {
        ## Position of variable m_{ik}.
        n_of_y_variables + i + (k - 1L) * no
    }

    ind_o <- seq_len(no)
    ind_c <- seq_len(nc)
    ind_y <- seq_len(n_of_y_variables)
    z <- as.matrix(expand.grid(ind_o, ind_o, ind_c))

    ## Build the three constraint objects:
    ind <- ind_y
    if(!sparse) {
        constraint_mat <- matrix(0, 3 * n_of_y_variables, n_of_variables)
        constraint_mat[cbind(ind, ind_y)] <- 1
        constraint_mat[cbind(ind, pos_m(z[, 1L], z[, 3L]))] <- -1
        constraint_mat[cbind(ind, pos_m(z[, 2L], z[, 3L]))] <- -1
        ind <- ind + n_of_y_variables
        constraint_mat[cbind(ind, ind_y)] <- 1
        constraint_mat[cbind(ind, pos_m(z[, 1L], z[, 3L]))] <- -1
        ind <- ind + n_of_y_variables
        constraint_mat[cbind(ind, ind_y)] <- 1
        constraint_mat[cbind(ind, pos_m(z[, 2L], z[, 3L]))] <- -1
    } else {
        len <- length(ind)
        constraint_mat <-
            simple_triplet_matrix(c(rep.int(ind, 3L),
                                    rep.int(ind + n_of_y_variables, 2L),
                                    rep.int(ind + 2 * n_of_y_variables,
                                            2L)),
                                  c(ind,
                                    pos_m(z[, 1L], z[, 3L]),
                                    pos_m(z[, 2L], z[, 3L]),
                                    ind,
                                    pos_m(z[, 1L], z[, 3L]),
                                    ind,
                                    pos_m(z[, 2L], z[, 3L])),
                                  c(rep.int(1, len),
                                    rep.int(-1, 2 * len),
                                    rep.int(1, len),
                                    rep.int(-1, len),
                                    rep.int(1, len),
                                    rep.int(-1, len)),
                                  3 * n_of_y_variables, n_of_variables)
    }

    constraint_dir <- c(rep.int(">=", n_of_y_variables),
                        rep.int("<=", 2 * n_of_y_variables))
    constraint_rhs <- c(rep.int(-1, n_of_y_variables),
                        rep.int(0, 2 * n_of_y_variables))
    ## (And don't forget that we need a binary stochastic matrix M.)
    constraint_mat <-
        rbind(constraint_mat,
              cbind(matrix(0, no, n_of_y_variables),
                    kronecker(rbind(rep.int(1, nc)), diag(1, no))))
    constraint_dir <- c(constraint_dir, rep.int("=", no))
    constraint_rhs <- c(constraint_rhs, rep.int(1, no))
    ## Handle possibly additional explicit constrains.
    ## (Which specify that pairs of objects are in relation or not.)
    acmaker <-
        .make_additional_constraint_maker_using_memberships_E(pos_m,
                                                              n_of_variables,
                                                              nc, sparse)
    if(!is.null(A <- control$constraints)) {
        A <- .canonicalize_additional_constraints(A, "E")
        add <- acmaker(A)
        constraint_mat <- rbind(constraint_mat, add$mat)
        constraint_dir <- c(constraint_dir, add$dir)
        constraint_rhs <- c(constraint_rhs, add$rhs)
    }

    labels <- dimnames(C)

    ## Set up augmented target function.
    objective_in <- c(rep.int(c(C), nc), double(n_of_m_variables))

    integer_positions <- seq_len(n_of_m_variables) + n_of_y_variables

    if(!is.null(all <- control$all) && identical(all, TRUE)) {
        verbose <- control$verbose
        if(is.null(verbose))
            verbose <- getOption("verbose")
        ## For finding all solutions, we loop over combinations of
        ## objects and classes:
        acmaker <-
            .make_additional_constraint_maker_using_o_c_pairs(pos_m,
                                                              n_of_variables,
                                                              nc, sparse)
        return(.find_all_relation_symdiff_optima_E_or_P_k(no,
                                                          nc,
                                                          "E",
                                                          objective_in,
                                                          constraint_mat,
                                                          constraint_dir,
                                                          constraint_rhs,
                                                          integer_positions,
                                                          A,
                                                          acmaker,
                                                          labels,
                                                          verbose,
                                                          control$solver,
                                                          control$control))
    }

    y <- solve_MILP(MILP(objective_in,
                         list(constraint_mat,
                              constraint_dir,
                              constraint_rhs),
                         integer_positions,
                         maximum = TRUE),
                    control$solver,
                    control$control)
    .stop_if_lp_status_is_nonzero(y$status, "E")

    M <- matrix(y$solution[integer_positions], ncol = nc)

    .make_incidence_from_class_memberships(M, "E", labels)
}

### * fit_relation_symdiff_P_k

fit_relation_symdiff_P_k <-
function(C, nc, control = list())
{
    ## Fit preference with at most nc classes.

    sparse <- control$sparse
    if(is.null(sparse)) sparse <- FALSE

    ## Using the membership matrix M = [m_{ik}] corresponding to the
    ## indifference relation (partition) with columns sorted according
    ## to increasing preference, we have to maximize
    ##   \sum_{i,j,k,l} c_{ij} I(k <= l) m_{ik} m_{jl}
    ## over all binary stochastic matrices M.  This is done by a simple
    ## linearization of the quadratic integer program.

    ## See above, but we need y_{ijkl} = m_{ik} m_{jl} and corresponding
    ## constraints.

    no <- nrow(C)

    n_of_y_variables <- no ^ 2 * nc ^ 2
    n_of_m_variables <- no * nc
    n_of_variables <- n_of_y_variables + n_of_m_variables

    pos_m <- function(i, k) {
        ## Position of variable m_{ik}.
        n_of_y_variables + i + (k - 1L) * no
    }

    ind_o <- seq_len(no)
    ind_c <- seq_len(nc)
    ind_y <- seq_len(n_of_y_variables)
    z <- as.matrix(expand.grid(ind_o, ind_o, ind_c, ind_c))

    ## Build the three constraint objects:
    ind <- ind_y
    if(!sparse) {
        constraint_mat <- matrix(0, 3 * n_of_y_variables, n_of_variables)
        ind <- ind_y
        constraint_mat[cbind(ind, ind_y)] <- 1
        constraint_mat[cbind(ind, pos_m(z[, 1L], z[, 3L]))] <- -1
        constraint_mat[cbind(ind, pos_m(z[, 2L], z[, 4L]))] <- -1
        ind <- ind + n_of_y_variables
        constraint_mat[cbind(ind, ind_y)] <- 1
        constraint_mat[cbind(ind, pos_m(z[, 1L], z[, 3L]))] <- -1
        ind <- ind + n_of_y_variables
        constraint_mat[cbind(ind, ind_y)] <- 1
        constraint_mat[cbind(ind, pos_m(z[, 2L], z[, 4L]))] <- -1
    }
    else {
        len <- length(ind)
        constraint_mat <-
            simple_triplet_matrix(c(rep.int(ind, 3L),
                                    rep.int(ind + n_of_y_variables, 2L),
                                    rep.int(ind + 2 * n_of_y_variables,
                                            2L)),
                                  c(ind,
                                    pos_m(z[, 1L], z[, 3L]),
                                    pos_m(z[, 2L], z[, 4L]),
                                    ind,
                                    pos_m(z[, 1L], z[, 3L]),
                                    ind,
                                    pos_m(z[, 2L], z[, 4L])),
                                  c(rep.int(1, len),
                                    rep.int(-1, 2 * len),
                                    rep.int(1, len),
                                    rep.int(-1, len),
                                    rep.int(1, len),
                                    rep.int(-1, len)),
                                  3 * n_of_y_variables, n_of_variables)
    }
    constraint_dir <- c(rep.int(">=", n_of_y_variables),
                        rep.int("<=", 2 * n_of_y_variables))
    constraint_rhs <- c(rep.int(-1, n_of_y_variables),
                        rep.int(0, 2 * n_of_y_variables))
    ## (And don't forget that we need a binary stochastic matrix M.)
    constraint_mat <-
        rbind(constraint_mat,
              cbind(matrix(0, no, n_of_y_variables),
                    kronecker(rbind(rep.int(1, nc)), diag(1, no))))
    constraint_dir <- c(constraint_dir, rep.int("=", no))
    constraint_rhs <- c(constraint_rhs, rep.int(1, no))
    ## Handle constraints on the class sizes.
    if(!is.null(l <- control$l)) {
        ## Verify that l is feasible.
        ## It would be nice to be nice, in particular by rescaling the
        ## sizes to sum to the number of objects (so that one can
        ## e.g. specify proportions).  But using something like
        ##   l <- round(l / sum(l) * nc)
        ## does not work (e.g., c(2.4, 2.4, 5.2) sums to 10 but after
        ## rounding only to 9).
        if((length(l) != nc) || (sum(l) != no))
            stop("Invalid specification of class sizes.")
        ## And now add constraints \sum_j m_{jk} = l_k, k = 1, ..., nc.
        constraint_mat <-
            rbind(constraint_mat,
                  cbind(matrix(0, nc, n_of_y_variables),
                        kronecker(diag(1, nc), rbind(rep.int(1, no)))))
        constraint_dir <- c(constraint_dir, rep.int("=", nc))
        constraint_rhs <- c(constraint_rhs, l)
    }
    ## Handle possibly additional explicit constrains.
    ## (Which specify that pairs of objects are in relation or not.)
    acmaker <-
        .make_additional_constraint_maker_using_memberships_P(pos_m,
                                                              n_of_variables,
                                                              nc, sparse)
    if(!is.null(A <- control$constraints)) {
        A <- .canonicalize_additional_constraints(A, "P")
        add <- acmaker(A)
        constraint_mat <- rbind(constraint_mat, add$mat)
        constraint_dir <- c(constraint_dir, add$dir)
        constraint_rhs <- c(constraint_rhs, add$rhs)
    }

    labels <- dimnames(C)

    ## Set up augmented target function, silly way.
    C <- array(C, c(dim(C), nc, nc))
    for(k in ind_c)
        for(l in seq_len(k - 1L))
            C[, , k, l] <- 0

    objective_in <- c(c(C), double(n_of_m_variables))

    integer_positions <- seq_len(n_of_m_variables) + n_of_y_variables

    if(!is.null(all <- control$all) && identical(all, TRUE)) {
        verbose <- control$verbose
        if(is.null(verbose))
            verbose <- getOption("verbose")
        ## For finding all solutions, we loop over combinations of
        ## objects and classes:
        acmaker <-
            .make_additional_constraint_maker_using_o_c_pairs(pos_m,
                                                              n_of_variables,
                                                              nc, sparse)
        return(.find_all_relation_symdiff_optima_E_or_P_k(no,
                                                          nc,
                                                          "P",
                                                          objective_in,
                                                          constraint_mat,
                                                          constraint_dir,
                                                          constraint_rhs,
                                                          integer_positions,
                                                          A,
                                                          acmaker,
                                                          labels,
                                                          verbose,
                                                          control$solver,
                                                          control$control))
    }

    y <- solve_MILP(MILP(objective_in,
                         list(constraint_mat,
                              constraint_dir,
                              constraint_rhs),
                         integer_positions,
                         maximum = TRUE),
                    control$solver,
                    control$control)
    .stop_if_lp_status_is_nonzero(y$status, "P")

    M <- matrix(y$solution[integer_positions], ncol = nc)

    .make_incidence_from_class_memberships(M, "P", labels)
}

### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "### [*]+" ***
### End: ***
