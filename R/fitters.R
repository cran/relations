### Relation fitters.

### * fit_relation_symdiff

## General purpose fitter for families
.SD_families <- c("E", "L", "O", "P", "T")
## where:
##   E ... equivalence relations               REF SYM TRA
##   L ... linear orders                   TOT  *  ASY TRA
##   O ... partial order                       REF ASY TRA
##   P ... complete preorder ("ordering")  TOT  *      TRA
##   T ... tournament                      TOT  *  ASY
## and TOT <-> total/complete, ASY <-> antisymmetric.
## Note that completeness implies reflexivity (indicated by '*').

## Number of non-redundant object pairs.
.n_of_pairs <-
function(n, family)
{
    family <- match.arg(family, .SD_families)
    N <- n * (n - 1)                    # Number of distinct pairs.
    switch(EXPR = family, E = , L = , T = N / 2, O = , P = N)
}

## Number of transitivity constraints for the non-redundant pairs.
## (None for tournaments.)
.n_of_transitivity_constraints <-
function(n, family)
{
    family <- match.arg(family, .SD_families)
    N <- n * (n - 1) * (n - 2)          # Number of distinct triples.
    switch(EXPR = family, E = N / 2, L = N / 3, O = , P = N, T = 0L)
}

## Make function giving the position of incidence (i, j) in the vector
## of non-redundant incidences used for symdiff fitting (upper.tri() for
## E/L/T and .offdiag() for P/O, respectively).
.make_pos <-
function(n, family)
{
    family <- match.arg(family, .SD_families)
    if(family %in% c("E", "L", "T"))
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
    ## We'll need this eventually, so let's try now.
    require("lpSolve")
    
    ## Number of objects:
    n <- nrow(x)
    
    objective_in <- if(family %in% c("L", "T"))
        (x - t(x))[upper.tri(x)]
    else if(family == "E")
        (x + t(x))[upper.tri(x)]
    else
        x[.offdiag(x)]

    ## Handle constraints implied by the family to be fitted.
    ## Need all variables in { 0 , 1 }, i.e., >= 0 and <= 1 and
    ## integer, plus maybe totality or antisymmetry, plus maybe
    ## transitivity.
    NP <- .n_of_pairs(n, family)
    family_is_O_or_P <- family %in% c("O", "P")
    constr_mat <-
        rbind(diag(1, NP),
              diag(1, NP),
              if(family_is_O_or_P)
              .make_tot_or_asy_constraint_mat(n),              
              .make_transitivity_constraint_mat(n, family))
    constr_dir <- 
        c(rep.int(">=", NP),
          rep.int("<=", NP),
          if(family_is_O_or_P)
          .make_tot_or_asy_constraint_dir(n, family),
          .make_transitivity_constraint_dir(n, family))
    constr_rhs <-
        c(rep.int(0, NP),
          rep.int(1, NP),
          if(family_is_O_or_P)          
          .make_tot_or_asy_constraint_rhs(n),
          .make_transitivity_constraint_rhs(n, family))

    ## Handle additional constraints.
    acmaker <-
        .make_additional_constraint_maker_using_incidences(n, family)
    if(!is.null(A <- control$constraints)) {
        A <- .canonicalize_additional_constraints(A, "P")
        add <- acmaker(A)
        constr_mat <- rbind(constr_mat, add$mat)
        constr_dir <- c(constr_dir, add$dir)
        constr_rhs <- c(constr_rhs, add$rhs)
    }

    labels <- dimnames(x)

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
                                                 verbose))
    }
    
    out <- lpSolve::lp("max", 
                       objective_in,
                       constr_mat,
                       constr_dir,
                       constr_rhs,
                       int.vec = seq_len(NP))
    .stop_if_lp_status_is_nonzero(out$status, family)

    ## Turn the solution back into a full incidence matrix.
    fit <- if(family_is_O_or_P)
        .make_incidence_from_offdiag(round(out$solution), family, labels)
    else
        .make_incidence_from_upper_tri(round(out$solution), family, labels)
    ## For the time being, tack some of the lp() results on so that we
    ## can look at them (but not everything due to size ...)
    attr(fit, ".lp") <- out[c("objval", "solution")]

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
         verbose = FALSE)
{
    ## Start by computing one optimal solution.
    out <- lpSolve::lp("max", obj, mat, dir, rhs, int.vec = int)
    ## Check status:
    .stop_if_lp_status_is_nonzero(out$status, family)    
    ## Value of the optimum:
    ## (Seems that out$objval is less precise than this.)
    Vopt <- sum(out$objective * out$solution)

    ## Helper function one.
    ## Compute optimum for a valid matrix A of additional constraints.
    optimize_with_valid_additional_constraints <- function(A) {
        add <- acmaker(A)
        out <- lpSolve::lp("max", obj,
                           rbind(mat, add$mat),
                           c(dir, add$dir),
                           c(rhs, add$rhs),
                           int.vec = int)
        if(out$status != 0) return(-Inf)
        ## Seems that out$objval is less precise than this:
        sum(out$objective * out$solution)
    }

    ## Helper function two.
    ## Figure out whether constraint x_{ij} = 0 or x_{ij} = 1, or both,
    ## are "optimal".
    add_constraint_for_single_pair <- function(A, i, j) {
        Aij0 <- rbind(A, c(i, j, 0))
        Aij1 <- rbind(A, c(i, j, 1))
        ## <NOTE>
        ## At least when using the objval given by lpSolve:lp(), we
        ## cannot reliably test v < Vopt due to numerical precision
        ## issues.  This seems to be less of an issue when using
        ##   sum(out$objective * out$solution)
        ## but then we've still seen cases where max(v0, v1) < Vopt,
        ## hence we compare with some tolerance (could improve ...).
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

    ## And now find all optimal additional constrains (i.e., at the end,
    ## all optimal incidences in 3-column matrix (i, j, x_{ij}) form) by
    ## looping over all non-redundant (i, j) pairs.
    ##
    ## Somewhat tricky: in case the fitting is to include additional
    ## constraints, we need to ensure that we only try to add feasible
    ## additional constraints.

    ## Determine a 2-column matrix of (i, j) index pairs to loop over.
    ind <- if(family %in% c("E", "L", "T")) {
        ## Loop over all pairs 1 <= i < j <= n.
        which(upper.tri(matrix(0, n, n)), arr.ind = TRUE)
    }
    else {
        ## Loop over all pairs 1 <= i != j <= n.
        which(.offdiag(matrix(0, n, n)), arr.ind = TRUE)
    }
    ## But exclude pairs already specified in additional constraints.
    if(!is.null(A)) {
        hpos <- if(family %in% c("E", "L", "T")) {
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
                                add_constraint_for_single_pair,
                                ind[k, 1L], ind[k, 2L]))
        if(verbose)
            cat("N_of_pairs:", k,
                "***",
                "N_of_optimal_branches:", length(Alist),
                "\n")
    }

    lapply(Alist, .make_incidence_from_triples, family, labels, n)
}

### * Constraint generators: transitivity.

### ** .make_transitivity_constraint_mat

.make_transitivity_constraint_mat <-
function(n, family)
{
    family <- match.arg(family, .SD_families)

    NP <- .n_of_pairs(n, family)
    if((n <= 2L) || (family == "T")) return(matrix(0, 0, NP))

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

        out <- matrix(0, NC, NP)
        ind <- seq_len(NC)        

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
            out[cbind(ind, c(p_ij, p_ij))] <- rep(c(1, -1), each = NT)
            out[cbind(ind, c(p_jk, p_jk))] <- rep(c(1, -1), each = NT)
            out[cbind(ind, c(p_ik, p_ik))] <- rep(c(-1, 1), each = NT)
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

        out <- matrix(0, NC, NP)
        ind <- seq_len(NC)
        out[cbind(ind, pos(z[, 1L], z[, 2L]))] <- 1
        out[cbind(ind, pos(z[, 2L], z[, 3L]))] <- 1
        out[cbind(ind, pos(z[, 1L], z[, 3L]))] <- -1
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
function(n)
{
    if(n <= 1) return(matrix(0, 0L, 0L))  # :-)

    ## Position of x_{ij} in x[.offdiag(x)]
    pos <- .make_pos(n, "P")

    ## Number of non-redundant pairs.
    NP <- n * (n - 1)
    ## Number of constraints.
    NC <- NP / 2

    ind <- seq_len(n)
    z <- as.matrix(expand.grid(ind, ind))[, c(2L, 1L)]
    z <- z[z[, 1L] < z[, 2L], ]
    out <- matrix(0, NC, NP)
    ind <- seq_len(NC)
    out[cbind(ind, pos(z[, 1L], z[, 2L]))] <- 1
    out[cbind(ind, pos(z[, 2L], z[, 1L]))] <- 1
    out
}

### ** .make_tot_or_asy_constraint_dir

.make_tot_or_asy_constraint_dir <-
function(n, family)
{
    NC <- n * (n - 1) / 2
    rep.int(switch(EXPR = family, O = "<=", P = ">="), NC)
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
function(n, family)
    function(A) {
        ## (This assumes A has already been validated and
        ## canonicalized.)
        na <- nrow(A)
        mat <- matrix(0, na, .n_of_pairs(n, family))
        pos <- .make_pos(n, family)
        mat[cbind(seq_len(na), pos(A[, 1L], A[, 2L]))] <- 1
        list(mat = mat,
             dir = rep.int("==", na),
             rhs = A[, 3L])
    }

### .make_additional_constraint_maker_using_memberships_E

.make_additional_constraint_maker_using_memberships_E <-
function(pos, n_of_variables, nc)
    function(A) {
        ## (This assumes A has already been validated and
        ## canonicalized.)
        na <- nrow(A)
        mat <- matrix(0, nc * na, n_of_variables)
        ## For equivalences (E):
        ## A constraint of the form (i, j, 1) implies that
        ## objects i and j are in the same class, i.e.:
        ##   m_{ik} - m_{jk} == 0   for all k.
        ## A constraint of the form (i, j, 0) implies that
        ## objects i and j are in different classes, i.e.:
        ##   m_{ik} + m_{jk} <= 1   for all k.
        ind <- seq_len(na)
        for(k in seq_len(nc)) {
            mat[cbind(ind, pos(A[, 1L], k))] <- 1
            mat[cbind(ind, pos(A[, 2L], k))] <-
                ifelse(A[, 3L] == 1, -1, 1)
            ind <- ind + na
        }
        list(mat = mat,
             dir = rep.int(ifelse(A[, 3L] == 1, "==", "<="), nc),
             rhs = rep.int(ifelse(A[, 3L] == 1, 0, 1), nc))
    }

### .make_additional_constraint_maker_using_memberships_P

.make_additional_constraint_maker_using_memberships_P <-
function(pos, n_of_variables, nc)
    function(A) {
        ## (This assumes A has already been validated and
        ## canonicalized.)
        na <- nrow(A)
        mat <- matrix(0, nc * na, n_of_variables)
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
                for(l in seq_len(k - 1))
                    mat[cbind(ind, pos(A[i0, 2L], l))] <- -1
                ind <- ind + na
            }
        }

        list(mat = mat,
             dir = rep.int("<=", nc * na),
             rhs = double(nc * na))
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
    ## For E/L/T, we use only pairs i < j, so swap and possibly flip
    ## (L/T) if needed.
    if(family %in% c("E", "L", "T")) {
        ind <- A[, 1L] > A[, 2L]
        if(any(ind))
            A[ind, ] <- cbind(A[ind, c(2L, 1L), drop = FALSE],
                              if(family == "E") A[ind, 3L]
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
function(x, family = c("E", "L", "T"), labels = NULL, n = NULL)
{
    ## Compute the indicences of a binary relation from the given family
    ## from its upper triangular part (provided this is possible, of
    ## course).

    family <- match.arg(family)

    ## Number of objects.
    n_from_tri_len <- function(n)
        (1 + sqrt(8 * n + 1)) / 2
    if(is.null(n))
        n <- n_from_tri_len(length(x))

    if(family == "E") {
        ## Equivalence relation.
        y <- diag(1 / 2, n, n)
        y[upper.tri(y)] <- x
        y <- y + t(y)
    }
    else {
        ## Linear order or tournament.
        y <- diag(1, n, n)
        y[upper.tri(y)] <- x
        y[lower.tri(y)] <- 1 - t(y)[lower.tri(y)]
    }

    dimnames(y) <- labels
    y
}

### ** .make_incidence_from_offdiag

.make_incidence_from_offdiag <-
function(x, family = c("O", "P"), labels = NULL, n = NULL)
{
    ## Compute the indicences of a binary relation from its off-diagonal
    ## part (provided this is possible, of course).

    if(is.null(n))
        n <- (1 + sqrt(4 * length(x) + 1)) / 2

    ## <NOTE>
    ## Diagonal entries for the incidences are a mess.
    ## For antisymmetric relations, it really does/should not matter.
    ## Nevertheless, we now (unlike Fishburn) include reflexivity in the
    ## definition of families L, O, P, and T. (For all but O, this is in
    ## fact implied by the standard definition of completeness.)
    ## </NOTE>

    y <- diag(1, n, n)
    y[.offdiag(y)] <- x
    
    dimnames(y) <- labels
    y
}

### ** .make_incidence_from_triples

.make_incidence_from_triples <-
function(x, family, labels = NULL, n = NULL)
{
    ## Compute the indicences of a binary relation from a 3-column
    ## matrix the rows of which are triples (i, j, x_{ij}) with x_{ij}
    ## the incidence at position (i, j).
    
    if(is.null(n))
        n <- max(x)
    ## Set up incidences for the non-redundant pairs.
    I <- matrix(0, n, n)
    I[x[, -3L, drop = FALSE]] <- x[, 3L]
    ## And complete according to family.
    diag(I) <- 1
    if(family %in% c("E", "L", "T")) {
        ind <- lower.tri(I)
        I[ind] <- if(family == "E")
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
        n_of_y_variables + i + (k - 1) * no
    }

    ind_o <- seq_len(no)
    ind_c <- seq_len(nc)
    ind_y <- seq_len(n_of_y_variables)
    z <- as.matrix(expand.grid(ind_o, ind_o, ind_c))
    
    ## Build the three constraint objects:
    constraint_mat <- matrix(0, 3 * n_of_y_variables, n_of_variables)
    ind <- ind_y
    constraint_mat[cbind(ind, ind_y)] <- 1    
    constraint_mat[cbind(ind, pos_m(z[, 1L], z[, 3L]))] <- -1
    constraint_mat[cbind(ind, pos_m(z[, 2L], z[, 3L]))] <- -1
    ind <- ind + n_of_y_variables
    constraint_mat[cbind(ind, ind_y)] <- 1
    constraint_mat[cbind(ind, pos_m(z[, 1L], z[, 3L]))] <- -1
    ind <- ind + n_of_y_variables
    constraint_mat[cbind(ind, ind_y)] <- 1
    constraint_mat[cbind(ind, pos_m(z[, 2L], z[, 3L]))] <- -1
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
    acmaker <-
        .make_additional_constraint_maker_using_memberships_E(pos_m,
                                                              n_of_variables,
                                                              nc)
    if(!is.null(A <- control$constraints)) {
        A <- .canonicalize_additional_constraints(A, "E")
        add <- acmaker(A)
        constraint_mat <- rbind(constraint_mat, add$mat)
        constraint_dir <- c(constraint_dir, add$dir)
        constraint_rhs <- c(constraint_rhs, add$rhs)
    }

    ## Set up augmented target function.
    objective_in <- c(rep.int(c(C), nc), double(n_of_m_variables))

    integer_positions <- seq_len(n_of_m_variables) + n_of_y_variables

    if(!is.null(all <- control$all) && identical(all, TRUE)) {
        verbose <- control$verbose
        if(is.null(verbose))
            verbose <- getOption("verbose")
        return(.find_all_relation_symdiff_optima(no, "E",
                                                 objective_in,
                                                 constraint_mat,
                                                 constraint_dir,
                                                 constraint_rhs,
                                                 integer_positions,
                                                 A,
                                                 acmaker,
                                                 dimnames(C),
                                                 verbose))
    }

    y <- lpSolve::lp("max",
                     objective_in,
                     constraint_mat,
                     constraint_dir,
                     constraint_rhs,
                     int.vec = integer_positions)
    .stop_if_lp_status_is_nonzero(y$status, "E")

    M <- matrix(y$solution[integer_positions], nc = nc)

    ## Return incidences.
    tcrossprod(M)
    
}

### * fit_relation_symdiff_P_k

fit_relation_symdiff_P_k <-
function(C, nc, control = list())
{
    ## Fit preference with at most nc classes.

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
        n_of_y_variables + i + (k - 1) * no
    }

    ind_o <- seq_len(no)
    ind_c <- seq_len(nc)
    ind_y <- seq_len(n_of_y_variables)
    z <- as.matrix(expand.grid(ind_o, ind_o, ind_c, ind_c))
    
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
    acmaker <-
        .make_additional_constraint_maker_using_memberships_P(pos_m,
                                                              n_of_variables,
                                                              nc)
    if(!is.null(A <- control$constraints)) {
        A <- .canonicalize_additional_constraints(A, "P")
        add <- acmaker(A)
        constraint_mat <- rbind(constraint_mat, add$mat)
        constraint_dir <- c(constraint_dir, add$dir)
        constraint_rhs <- c(constraint_rhs, add$rhs)
    }

    ## Set up augmented target function, silly way.
    C <- array(C, c(dim(C), nc, nc))
    for(k in ind_c)
        for(l in seq(from = 1, length = k - 1))
            C[, , k, l] <- 0
    
    objective_in <- c(c(C), double(n_of_m_variables))

    integer_positions <- seq_len(n_of_m_variables) + n_of_y_variables

    if(!is.null(all <- control$all) && identical(all, TRUE)) {
        verbose <- control$verbose
        if(is.null(verbose))
            verbose <- getOption("verbose")
        return(.find_all_relation_symdiff_optima(no, "P",
                                                 objective_in,
                                                 constraint_mat,
                                                 constraint_dir,
                                                 constraint_rhs,
                                                 integer_positions,
                                                 A,
                                                 acmaker,
                                                 dimnames(C),
                                                 verbose))
    }

    y <- lpSolve::lp("max",
                     objective_in,
                     constraint_mat,
                     constraint_dir,
                     constraint_rhs,
                     int.vec = integer_positions)
    .stop_if_lp_status_is_nonzero(y$status, "P")

    M <- matrix(y$solution[integer_positions], nc = nc)

    ## Return incidences.
    E <- matrix(0, nc, nc)
    E[row(E) <= col(E)] <- 1
    tcrossprod(M %*% E, M)
}

### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "### [*]+" ***
### End: ***
