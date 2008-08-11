relation_choice <-
function(x, method = "symdiff", control = list(), ...)
{
    dots <- list(...)
    control[names(dots)] <- dots

    relations <- as.relation_ensemble(x)

    if(!length(relations))
        stop("Cannot compute choice from empty ensemble.")

    known_methods <-
        list("symdiff" = ".relation_choice_symdiff")
    if(is.character(method)) {
        ## Hopefully of length one, add some tests eventually ...
        if(is.na(ind <- pmatch(method, names(known_methods))))
            stop(gettextf("Method '%s' is not a valid choice method."))
        method <- get(known_methods[[ind]][1L])
    }
    else if(!is.function(method))
        stop("Invalid argument 'method'.")

    method(relations, control)
}

## <FIXME>
## Add information on which problem is solved, and how this is done.
## </FIXME>

.relation_choice_symdiff <-
function(relations, control)
{
    ## <FIXME>
    ## Shouldn't this be true for all choice functions?
    ## (As we choose from given preferences on a set of objects.)
    ## </FIXME>
    if(!.is_ensemble_of_endorelations(relations))
        stop("Need an ensemble of endorelations.")
    if(!.is_ensemble_of_crisp_relations(relations))
        stop("Need an ensemble of crisp relations.")

    ## Argument handling.
    k <- control$k
    if(is.null(k)) k <- 1L              # Single winner by default.
    all <- control$all
    if(is.null(all)) all <- FALSE
    sparse <- control$sparse
    if(is.null(sparse)) sparse <- FALSE
    solver <- control$solver
    control <- control$control
    ## And of course, more sanity checking would be nice ...

    x <- .make_fit_relation_symdiff_M(relations, 1)
    n <- nrow(x)
    ## Set diagonal to 0.
    diag(x) <- 0
    ## Compute coefficients of target function.
    beta <- colSums(pmin(x, 0)) - rowSums(pmax(x, 0))
    x <- abs(x)
    alpha <- (x + t(x))[upper.tri(x)]
    obj <- c(alpha, beta)
    ## Constraints.
    C <- .make_symdiff_winners_constraints(n, k, sparse = sparse)
    ## Underlying set of objects to choose from.
    ## (Domain of the choice problem.)
    D <- as.list(.domain(relations)[[1L]])
    ## Solve.
    if(all)
        return(.find_all_relation_symdiff_winners(alpha, beta, C, n, D,
                                                  solver, control))
    binary_positions <- .n_of_symdiff_choice_v_variables(n) + seq_len(n)
    out <- solve_MILP(MILP(obj, list(C$mat, C$dir, C$rhs),
                           types = .make_types(length(obj),
                                               B = binary_positions),
                           maximum = TRUE),
                      solver, control)
    u <- out$solution[binary_positions]
    as.set(D[u == 1])
}

## Number of v variables for symdiff winners.
.n_of_symdiff_choice_v_variables <- function(n) n * (n - 1) / 2
## Number of u variables for symdiff winners.
.n_of_symdiff_choice_u_variables <- function(n) n

.make_symdiff_winners_constraints <-
function(n, k, sparse = FALSE)
{
    ## Constraints:
    ##   v_{ij} - u_i       <= 0
    ##   v_{ij}       - u_j <= 0
    ##   v_{ij} - u_i - u_j >= -1
    ##            u_i       <= 1
    ##            sum(u_i)  == k

    n_u <- n
    n_v <- .n_of_symdiff_choice_v_variables(n)

    ## Create a matrix with all combinations of triples i < j in the
    ## rows (note that we use the same order as in upper.tri()).
    ind <- seq_len(n)
    z <- as.matrix(expand.grid(ind, ind))
    z <- z[z[, 1L] < z[, 2L], , drop = FALSE]
    ## Constraint mat.
    ## Note that we v_ij are arranged in the same order as z, so that
    ## the v_ij constraint parts are blocks of identity matrices.
    ind <- p_v_ij <- seq_len(n_v)
    p_u_i <- n_v + z[, 1L]
    p_u_j <- n_v + z[, 2L]
    mat_nr <- 3L * n_v + n_u + 1L
    mat_nc <- n_v + n_u
    if(!sparse) {
        mat <- matrix(0, mat_nr, mat_nc)
        ##   v_{ij} - u_i       <= 0
        mat[cbind(ind, p_v_ij)] <- 1
        mat[cbind(ind, p_u_i)] <- -1
        ##   v_{ij}       - u_j <= 0
        ind <- ind + n_v
        mat[cbind(ind, p_v_ij)] <- 1
        mat[cbind(ind, p_u_j)] <- -1
        ##   v_{ij} - u_i - u_j >= -1
        ind <- ind + n_v
        mat[cbind(ind, p_v_ij)] <- 1
        mat[cbind(ind, p_u_i)] <- -1
        mat[cbind(ind, p_u_j)] <- -1
        ##            u_i       <= 1
        mat[cbind(3L * n_v + seq_len(n_u), n_v + seq_len(n_u))] <- 1
        ##            sum(u_i)  == k
        mat[3L * n_v + n_u + 1L, n_v + seq_len(n_u)] <- 1
    }
    else {
        mat_i <- c(rep.int(ind, 2L),
                   rep.int(ind + n_v, 2L),
                   rep.int(ind + 2L * n_v, 3L),
                   3L * n_v + seq_len(n_u),
                   rep.int(3L * n_v + n_u + 1L, n_u))
        mat_j <- c(p_v_ij, p_u_i,
                   p_v_ij, p_u_j,
                   p_v_ij, p_u_i, p_u_j,
                   rep.int(n_v + seq_len(n_u), 2L))
        mat_v <- c(rep.int(1, n_v), rep.int(-1, n_v),
                   rep.int(1, n_v), rep.int(-1, n_v),
                   rep.int(1, n_v), rep.int(-1, 2L * n_v),
                   rep.int(1, 2L * n_u))
        mat <- simple_triplet_matrix(mat_i, mat_j, mat_v,
                                     mat_nr, mat_nc)
    }
    ## Constraint dir.
    dir <- c(rep.int("<=", 2L * n_v),
             rep.int(">=", n_v),
             rep.int("<=", n_u),
             "==")
    ## Constraint rhs.
    rhs <- c(rep.int(0, 2L * n_v),
             rep.int(-1, n_v),
             rep.int(1, n_u),
             k)
    list(mat = mat, dir = dir, rhs = rhs)
}

.find_all_relation_symdiff_winners <-
function(alpha, beta, C, n, domain, solver, control)
{
    ## Find all relation symdiff "winners" by a simple branch-and-cut
    ## strategy, using successive reductions.

    node <- function(u, beta, value = 0, done = FALSE)
        list(u = u, beta = beta, value = value, done = done)

    bin_pos <- function(n)
        .n_of_symdiff_choice_v_variables(n) + seq_len(n)

    ## Constraints for the full problem.
    mat <- C$mat
    dir <- C$dir
    rhs <- C$rhs

    ## Value function for given problem instance.
    V <- function(obj, mat, dir, rhs, pos) {
        out <- solve_MILP(MILP(obj, list(mat, dir, rhs),
                               types = .make_types(length(obj), B = pos),
                               maximum = TRUE),
                          solver, control)
        if(out$status != 0) return(-Inf)
        out$objval
    }
    ## Determine the optimal value.
    Vopt <- V(c(alpha, beta), mat, dir, rhs, bin_pos(n))
    ## <FIXME>
    ## Check status.
    ## </FIXME>

    k <- rhs[length(rhs)]

    splitter <- function(node, n, alpha, delta, mat, dir, rhs) {
        ## This to be called for n > 1.
        if(node$done) return(node)
        tol <- 1e-10
        len <- length(rhs)
        pos <- bin_pos(n - 1L)
        u <- node$u
        beta <- node$beta
        value <- node$value
        ## Solution with last u = 0.
        beta_0 <- beta[-n]
        rhs_0 <- rhs
        rhs_0[len] <- k - sum(u)
        value_0 <- value
        node_0 <- node(c(u, 0), beta_0, value_0)
        ## Solution with last u = 1.
        beta_1 <- beta_0 + delta
        rhs_1 <- rhs_0
        rhs_1[len] <- rhs_0[len] - 1
        value_1 <- value_0 + beta[n]
        node_1 <- node(c(u, 1), beta_1, value_1)
        ## Try u = 0.
        v <- V(c(alpha, beta_0), mat, dir, rhs_0, pos)
        if(value_0 + v < Vopt - tol) return(list(node_1))
        ## Try u = 1.
        v <- V(c(alpha, beta_1), mat, dir, rhs_1, pos)
        if(value_1 + v < Vopt - tol) return(list(node_0))
        ## If both are optimal, branch:
        list(node_0, node_1)
    }

    finisher <- function(node) {
        ## This to be called for n = 1.
        ## Remeber we proceed from right to left and leave one off.
        u <- node$u
        u <- rev(c(u, k - sum(u)))
        as.set(domain[u == 1])
    }

    ## Main loop.
    nodes <- list(node(integer(), beta))
    while(n > 1L) {
        ## Move from n to n - 1.
        ## Reduce constraints.
        n_v <- .n_of_symdiff_choice_v_variables(n)
        ind <- seq(from = n_v - n + 2L, to = n_v)
        ## Drop elements from objective.
        delta <- alpha[ind]
        alpha <- alpha[-ind]
        ## Drop elements from constraints.
        rind <- c(ind, ind + n_v, ind + 2L * n_v, 3L * n_v + n)
        mat <- mat[-rind, -c(ind, ncol(mat)), drop = FALSE]
        dir <- dir[-rind]
        rhs <- rhs[-rind]
        ## Branch and cut.
        nodes <- do.call("c",
                         lapply(nodes, splitter, n, alpha, delta, mat,
                                dir, rhs))
        ## cat("n:", n, "n of nodes", length(nodes), "\n")
        n <- n - 1L
    }

    nodes <- lapply(nodes, finisher)

    nodes
}

## And now use along the lines of
##   data("SVM_Benchmarking_Regression")
##   relation_choice(SVM_Benchmarking_Regression, k = 1)
##   relation_choice(SVM_Benchmarking_Regression, k = 2)
## etc.
