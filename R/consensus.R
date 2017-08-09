### Consensus relations.

### * relation_consensus

relation_consensus <-
function(x, method = NULL, weights = 1, control = list(), ...)
{
    dots <- list(...)
    control[names(dots)] <- dots

    if(inherits(x, "gset")) {
        relations <- relation_ensemble(list = gset_support(x))
        weights <- gset_memberships(x)
    } else {
        relations <- as.relation_ensemble(x)
    }

    if(!length(relations))
        stop("Cannot compute consensus of empty ensemble.")

    weights <- rep_len(weights, length(relations))
    if(any(weights < 0))
        stop("Argument 'weights' has negative elements.")
    if(!any(weights > 0))
        stop("Argument 'weights' has no positive elements.")

    if(!is.function(method)) {
        if(!inherits(method, "relation_consensus_method")) {
            ## Get the method definition from the registry.
            if(!is.character(method) || (length(method) != 1L))
                stop("Invalid 'method' argument.")
            entry <- get_relation_consensus_method(method)
            if(is.null(entry))
                stop(gettextf("Method '%s' is not a valid consensus method.",
                              method),
                     domain = NA)
            method <- entry
        }
        method <- method$definition
    }

    method(relations, weights, control)
}

### * Relation consensus "methods"

### ** .relation_consensus_Borda

## Kendall/Borda method
.relation_consensus_Borda <-
function(relations, weights, control)
    .relation_consensus_score(relations, weights, control,
                              relation_scores, "Borda")

### ** .relation_consensus_Copeland

## Copeland method
.relation_consensus_Copeland <-
function(relations, weights, control)
    .relation_consensus_score(relations, weights, control,
                              relation_scores, "Copeland")

### ** .relation_consensus_score

.relation_consensus_score <-
function(relations, weights, control, FUN, ...)
{
    ## Several sanity checks could be done here.
    ## In particular, check whether all relations are in fact complete
    ## preferences (or whatever is really necessary).

    if(!.is_ensemble_of_endorelations(relations))
        stop("Need an ensemble of endorelations.")

    ## extract options.
    n <- .n_from_control_list(control)
    f_l_o <- isTRUE(control$L) ## force linear order?

    ## Get the scores.
    scores <- lapply(relations, FUN, ...)

    ## Multiply by the weights and compute the total scores.
    scores <- rowSums(mapply("*", scores, weights))

    ## break ties directly if a single linear order is enforced.
    out <-
        rank(scores,
             ties.method = if (f_l_o && n == 1L) "first" else "average")

    INC <- function(S) outer(S, S, "<=")
    I <- if (f_l_o && n > 1L) {
        ## find all groupwise permutations & combine
        l <- expand.grid(lapply(split(seq_along(out), out), .permute))

        ## Only use up to n combinations
        l <- l[seq_len(min(n, nrow(l))),,drop = FALSE]

        ## create incidences
        unlist(apply(l, 1, function(i) list(INC(unlist(i)))), recursive = FALSE)
    } else INC(out)

    meta <- list(is_endorelation = TRUE,
                 is_complete = TRUE,
                 is_reflexive = TRUE,
                 is_antisymmetric = f_l_o || !any(duplicated(out)),
                 is_transitive = TRUE,
                 scores = scores)
    .make_consensus_from_incidences(.domain(relations), I, meta)
}

### ** .relation_consensus_majority

.relation_consensus_majority <-
function(relations, weights, control)
{
    p <- control$p
    if(is.null(p)) p <- 1 / 2
    ## Add some sanity checking for p eventually.

    incidences <- lapply(relations, relation_incidence)
    weights <- rep_len(weights, length(relations))

    I <- (.weighted_sum_of_arrays(incidences, weights, na.rm = TRUE) /
          .weighted_sum_of_arrays(lapply(incidences, is.finite),
                                  weights))
    I <- if(p == 1)
        I == 1
    else
        I > p

    ## (Could add a tie-breaking mechanism a la Condorcet, with an
    ## option for finding all solutions.)

    .make_relation_from_domain_and_incidence(.domain(relations), I)
}

### ** .relation_consensus_CS

.relation_consensus_CS <-
function(relations, weights, control)
{
    ## Cook and Seiford, Management Science (1978).
    ## Determine a linear order minimizing the aggregate Cook-Seiford
    ## dissimilarity to a given ensemble of relations (originally:
    ## complete rankings [i.e., preferences]).
    ## This can be done by solving a linear sum assignment problem: the
    ## sought linear order is uniquely characterized by its ranks (r_i),
    ## and the target function is
    ##   \sum_b w_b \sum_i |r_i(b) - r_i| = \sum_{i,k} x_{ik} c_{ik}
    ## where
    ##   c_{ik} = \sum_b w_b | r_i(b) - k |
    ## and x_{ik} is one iff r_i is k.
    ## Clearly, this can be generalized to arbitrary score-based
    ## dissimilarities based on a score function which gives the same
    ## range of values for arbitrary linear orders.

    if(!.is_ensemble_of_endorelations(relations))
        stop("Need an ensemble of endorelations.")

    nos <- .n_from_control_list(control)

    n <- relation_size(relations)[1L]
    C <- matrix(0, n, n)
    ## Note that Cook and Seiford use Kendall-style "ranks" which are
    ## sorted in decreasing preference, whereas our default scores work
    ## in the opposite direction.
    incidences <- lapply(relations, relation_incidence)
    scores <- sapply(incidences, .incidence_scores_ranks)
    for(k in seq_len(n))
        C[, k] <- rowSums(sweep(abs(scores - k), 2L, weights, "*"))
    .compare <- function(u) outer(u, u, ">=")
    I <- if(nos > 1L)
        lapply(.find_up_to_n_LSAP_solutions(C, nos), .compare)
    else
        .compare(clue::solve_LSAP(C))
    objval <- .relation_consensus_CS_objval(I, incidences, weights)
    meta <- c(.relation_meta_db[["L"]], list(objval = objval))

    .make_consensus_from_incidences(.domain(relations), I, meta)
}

## Consensus methods for central relations using symdiff, Manhattan or
## Euclidean distance.  Note that
##
## * Symdiff dissimilarity only applies to crisp relation, and agrees
##   with Manhattan dissimilarity for these.
##
## * The restricted symdiff fitters are thus only for crisp ensembles.
##
## * We have restricted Manhattan fitters only for crisp ensembles (the
##   symdiff case).  Fitters for fuzzy ensembles could be added by the
##   usual means of turning an l_1 problem with linear/integer
##   constraints into a MILP (see e.g. CLUE).
##
## * The restricted symdiff fitters can also be used for determining
##   restricted *Euclidean* consensus relations for arbitrary (fuzzy)
##   ensembles.  We accomodate for this by calling the internal work
##   horses with a parameter indicating the (symdiff or Euclidean)
##   "context".
##
## * The restricted Euclidean fitters always give crisp relations.
##   Adding fitters for restricted fuzzy consensus would be possible via
##   SUMT approaches (or maybe these could be formulated as mixed
##   integer quadratic programs, but would this help?).

### ** .relation_consensus_Condorcet

## Symdiff/Manhattan crisp consensus relation for an ensemble of crisp
## relations.

.relation_consensus_Condorcet <-
function(relations, weights, control)
{
    if(!.is_ensemble_of_crisp_relations(relations))
        stop("Need an ensemble of crisp relations.")

    nos <- .n_from_control_list(control)

    incidences <- lapply(relations, relation_incidence)
    M <- .make_fit_relation_symdiff_B(incidences, weights, TRUE)
    I <- M >= 0                         # We do not break ties (>=).
    objval <- .relation_consensus_symdiff_objval(I, incidences, weights)
    meta <- list(is_endorelation = TRUE,
                 is_complete = TRUE,
                 objval = objval)
    if(nos == 1L) {
        meta <- c(meta,
                  list(is_reflexive = all(diag(M) >= 0),
                       is_antisymmetric = all(M != 0)))
        ## According to the way the default solution is defined.
        ## We do not know about transitivity without computing ...
    } else {
        I <- list(I)
        ind <- which(M == 0, arr.ind = TRUE)
        ## Recursively generate up to nos solutions by using 0 or 1 for
        ## the zero entries of M.
        splitter <- function(x, i, j) {
            y <- x
            ## By default we use 1 for the zero entries of M.
            y[i, j] <- 0
            list(x, y)
        }
        k <- 1L
        nr <- nrow(ind)
        while((length(I) < nos) && (k <= nr)) {
            I <- do.call("c",
                         lapply(I, splitter, ind[k, 1L], ind[k, 2L]))
            k <- k + 1L
        }
        if(length(I) > nos)
            I <- I[seq_len(nos)]
    }

    .make_consensus_from_incidences(.domain(relations), I, meta)
}

### ** .relation_consensus_symdiff_A

.relation_consensus_symdiff_A <-
function(relations, weights, control)
    .relation_consensus_symdiff(relations, "A", weights, control)

### ** .relation_consensus_symdiff_C

.relation_consensus_symdiff_C <-
function(relations, weights, control)
    .relation_consensus_symdiff(relations, "C", weights, control)

### ** .relation_consensus_symdiff_E

.relation_consensus_symdiff_E <-
function(relations, weights, control)
{
    k <- control$k
    if(!is.null(k))
        .relation_consensus_symdiff_E_k(relations, weights, k, control)
    else
        .relation_consensus_symdiff(relations, "E", weights, control)
}

### ** .relation_consensus_symdiff_L

.relation_consensus_symdiff_L <-
function(relations, weights, control)
    .relation_consensus_symdiff(relations, "L", weights, control)

### ** .relation_consensus_symdiff_M

.relation_consensus_symdiff_M <-
function(relations, weights, control)
    .relation_consensus_symdiff(relations, "M", weights, control)

### ** .relation_consensus_symdiff_O

.relation_consensus_symdiff_O <-
function(relations, weights, control)
    .relation_consensus_symdiff(relations, "O", weights, control)

### ** .relation_consensus_symdiff_W

.relation_consensus_symdiff_W <-
function(relations, weights, control)
{
    k <- control$k
    if(!is.null(k))
        .relation_consensus_symdiff_W_k(relations, weights, k, control)
    else
        .relation_consensus_symdiff(relations, "W", weights, control)
}

### ** .relation_consensus_symdiff_S

.relation_consensus_symdiff_S <-
function(relations, weights, control)
    .relation_consensus_symdiff(relations, "S", weights, control)

### ** .relation_consensus_symdiff_T

.relation_consensus_symdiff_T <-
function(relations, weights, control)
    .relation_consensus_symdiff(relations, "T", weights, control)

### ** .relation_consensus_symdiff_preorder

.relation_consensus_symdiff_preorder <-
function(relations, weights, control)
    .relation_consensus_symdiff(relations, "preorder", weights, control)

### ** .relation_consensus_symdiff_transitive

.relation_consensus_symdiff_transitive <-
function(relations, weights, control)
    .relation_consensus_symdiff(relations, "transitive", weights, control)

### ** .relation_consensus_manhattan

## Manhattan valued consensus relation for an ensemble of valued
## relations.

.relation_consensus_manhattan <-
function(relations, weights, control)
{
    incidences <- lapply(relations, relation_incidence)
    weights <- rep_len(weights, length(incidences))
    ## Incidences of the consensus relation are the weighted medians.
    I <- array(apply(do.call("cbind", lapply(incidences, c)),
                     1L, clue:::weighted_median, weights),
               dim = .size(relations))
    meta <-
        list(objval =
             .relation_consensus_manhattan_objval(I, incidences, weights))
    .make_relation_from_domain_and_incidence(.domain(relations), I, meta)
}

### ** .relation_consensus_euclidean

## Euclidean valued consensus relation for an ensemble of valued
## relations.

.relation_consensus_euclidean <-
function(relations, weights, control)
{
    weights <- rep_len(weights, length(relations))
    weights <- weights / sum(weights)
    incidences <- lapply(relations, relation_incidence)
    ## Incidences of the consensus relation are the weighted means.
    I <- .weighted_sum_of_arrays(incidences, weights)
    meta <-
        list(objval =
             .relation_consensus_euclidean_objval(I, incidences, weights))
    .make_relation_from_domain_and_incidence(.domain(relations), I, meta)
}

### ** .relation_consensus_euclidean_A

.relation_consensus_euclidean_A <-
function(relations, weights, control)
    .relation_consensus_symdiff(relations, "A", weights, control, TRUE)

### ** .relation_consensus_euclidean_C

.relation_consensus_euclidean_C <-
function(relations, weights, control)
    .relation_consensus_symdiff(relations, "C", weights, control, TRUE)

### ** .relation_consensus_euclidean_E

.relation_consensus_euclidean_E <-
function(relations, weights, control)
{
    k <- control$k
    if(!is.null(k))
        .relation_consensus_symdiff_E_k(relations, weights, k, control, TRUE)
    else
        .relation_consensus_symdiff(relations, "E", weights, control, TRUE)
}

### ** .relation_consensus_euclidean_L

.relation_consensus_euclidean_L <-
function(relations, weights, control)
    .relation_consensus_symdiff(relations, "L", weights, control, TRUE)

### ** .relation_consensus_euclidean_M

.relation_consensus_euclidean_M <-
function(relations, weights, control)
    .relation_consensus_symdiff(relations, "M", weights, control, TRUE)

### ** .relation_consensus_euclidean_O

.relation_consensus_euclidean_O <-
function(relations, weights, control)
    .relation_consensus_symdiff(relations, "O", weights, control, TRUE)

### ** .relation_consensus_euclidean_W

.relation_consensus_euclidean_W <-
function(relations, weights, control)
{
    k <- control$k
    if(!is.null(k))
        .relation_consensus_symdiff_W_k(relations, weights, k, control, TRUE)
    else
        .relation_consensus_symdiff(relations, "W", weights, control, TRUE)
}

### ** .relation_consensus_euclidean_S

.relation_consensus_euclidean_S <-
function(relations, weights, control)
    .relation_consensus_symdiff(relations, "S", weights, control, TRUE)

### ** .relation_consensus_euclidean_T

.relation_consensus_euclidean_T <-
function(relations, weights, control)
    .relation_consensus_symdiff(relations, "T", weights, control, TRUE)

### ** .relation_consensus_euclidean_preorder

.relation_consensus_euclidean_preorder <-
function(relations, weights, control)
    .relation_consensus_symdiff(relations, "preorder",
                                weights, control, TRUE)

### ** .relation_consensus_euclidean_transitive

.relation_consensus_euclidean_transitive <-
function(relations, weights, control)
    .relation_consensus_symdiff(relations, "transitive",
                                weights, control, TRUE)

## Consensus methods for central relations using CKS distance.

### ** .relation_consensus_CKS_A

.relation_consensus_CKS_A <-
function(relations, weights, control)
    .relation_consensus_CKS(relations, "A", weights, control)

### ** .relation_consensus_CKS_C

.relation_consensus_CKS_C <-
function(relations, weights, control)
    .relation_consensus_CKS(relations, "C", weights, control)

### ** .relation_consensus_CKS_E

.relation_consensus_CKS_E <-
function(relations, weights, control)
{
    k <- control$k
    if(!is.null(k))
        .relation_consensus_CKS_E_k(relations, weights, k, control)
    else
        .relation_consensus_CKS(relations, "E", weights, control)
}

## ** .relation_consensus_CKS_G

.relation_consensus_CKS_G <-
function(relations, weights, control)
{
    ## Generalized (Cook-Kress-Seiford) majority.

    incidences <- lapply(relations, relation_incidence)

    M <- .make_fit_relation_symdiff_B(incidences, weights, TRUE)
    Q <- .make_fit_relation_CKS_Q(incidences, weights)
    w <- sum(weights)

    ## See the 'Relations: Issues' writeup for the theory.
    ## With
    ##    m_{ij} = [M]_{ij} with M as above
    ##    q_{ij} = [Q]_{ij} with Q as above 
    ## we have
    ##    m_{ij} = 2 w (\bar{x}_{ij} - 1/2
    ##    q_{ij} = 2 w (\bar{q}_{ij} - 1/2
    ## or equivalently
    ##   \bar{x}_{ij} = 1/2 + m_{ij} / (2 w)
    ##   \bar{q}_{ij} = 1/2 + q_{ij} / (2 w)
    ## and hence
    ##    \bar{x}_{ij} \ge 1/2
    ##      .iff. m_{ij} \ge 0
    ##    \bar{q}_{ij} < 1/2
    ##      .iff. q_{ij} < 0
    ##    \bar{x}_{ij} \ge 1/2 - \bar{q}_{ij}
    ##      .iff. m_{ij} + q_{ij} \ge -w
    I <- (M >= 0) | ((Q < 0) & (M + Q >= - w))
    objval <- .relation_consensus_CKS_objval(I, incidences, weights)
    meta <- list(is_endorelation = TRUE, objval = objval)

    nos <- .n_from_control_list(control)

    if(nos > 1L) {
        I <- list(I)

        ## Split diagonal terms when m_{ii} = 0.
        for(i in which(diag(M) == 0)) {
            I <- do.call("c",
                         lapply(I,
                                function(x, i) {
                                    y <- x
                                    y[i, i] <- 0
                                    list(x, y)
                                },
                                i))
            if(length(I) >= nos) break
        }
        len <- length(I)
        if(len > nos)
            I <- I[seq_len(nos)]
        else if(len < nos) {
            ## Splitting non-diagonal terms is somewhat tricky as we
            ## cannot independently split x_{ij}/x_{ji} for split
            ## positions.
            ## Theory shows that the optimum is characterized via
            ##   \max(q_{ij}, p_{ij}, p_{ji}, 0)
            ## corresponding to 0/0, 1/0, 0/1 and 1/1 for x_{ij}/x_{ji}.
            ## Splits are characterized by
            ##   \max(q_{ij}, p_{ij}, p_{ji}) = 0.
            ## <FIXME>
            ## Check whether the above generalized majority rule always
            ## takes indicidence as one when possible.
            ## </CHECK>
            P <- .make_fit_relation_CKS_P(incidences, weights)
            ## Only need to know where P and Q are zero.
            P <- (P == 0)
            Q <- (Q == 0)
            ## And maybe whether i < j.
            U <- row(P) < col(P)
            ## Helper.
            do_split <- function(I, fun, ind) {
                k <- 1L
                nr <- nrow(ind)
                while((length(I) < nos) && (k <= nr)) {
                    I <- do.call("c",
                                 lapply(I, fun, ind[k, 1L], ind[k, 2L]))
                    k <- k + 1L
                }
                I
            }
            ## If q_{ij} = p_{ij} = 0, can take 0/0 1/0 1/1.
            I <- do_split(I,
                          function(x, i, j) {
                              z <- y <- x
                              x[i, j] <- x[j, i] <- 0
                              y[i, j] <- 1; y[j, i] <- 0
                              z[i, j] <- z[j, i] <- 1
                              list(x, y, z)
                          },
                          which(P & Q, arr.ind = TRUE))
            ## If only q_{ij} = 0, can take 0/0 1/1.
            I <- do_split(I,
                          function(x, i, j) {
                              y <- x
                              x[i, j] <- x[j, i] <- 0
                              y[i, j] <- y[j, i] <- 1
                              list(x, y)
                          },
                          which((Q & !P) & !t(P), arr.ind = TRUE))
            ## If p_{ij} = p_{ji} = 0, can take 1/0 0/1 1/1.
            I <- do_split(I,
                          function(x, i, j) {
                              z <- y <- x
                              x[i, j] <- 1; x[j, i] <- 0
                              y[i, j] <- 0; y[j, i] <- 1
                              z[i, j] <- z[j, i] <- 1
                              list(x, y, z)
                          },
                          which((P & t(P)) & U, arr.ind = TRUE))
            ## If only p_{ij} = 0, can take 1/0 1/1.
            I <- do_split(I,
                          function(x, i, j) {
                              y <- x
                              x[i, j] <- 1; x[j, i] <- 0
                              y[i, j] <- y[j, i] <- 1
                              list(x, y)
                          },
                          which((P & !Q) & !t(P), arr.ind = TRUE))
            ## Do this only once and not inside do_split.
            if(length(I) > nos)
                I <- I[seq_len(nos)]
        }
    }

    .make_consensus_from_incidences(.domain(relations), I, meta)
}

### ** .relation_consensus_CKS_L

.relation_consensus_CKS_L <-
function(relations, weights, control)
    .relation_consensus_CKS(relations, "L", weights, control)

### ** .relation_consensus_CKS_M

.relation_consensus_CKS_M <-
function(relations, weights, control)
    .relation_consensus_CKS(relations, "M", weights, control)

### ** .relation_consensus_CKS_O

.relation_consensus_CKS_O <-
function(relations, weights, control)
    .relation_consensus_CKS(relations, "O", weights, control)

### ** .relation_consensus_CKS_S

.relation_consensus_CKS_S <-
function(relations, weights, control)
    .relation_consensus_CKS(relations, "S", weights, control)

### ** .relation_consensus_CKS_T

.relation_consensus_CKS_T <-
function(relations, weights, control)
    .relation_consensus_CKS(relations, "T", weights, control)

### ** .relation_consensus_CKS_W

.relation_consensus_CKS_W <-
function(relations, weights, control)
{
    k <- control$k
    if(!is.null(k))
        .relation_consensus_CKS_W_k(relations, weights, k, control)
    else
        .relation_consensus_CKS(relations, "W", weights, control)
}

### ** .relation_consensus_CKS_preorder

.relation_consensus_CKS_preorder <-
function(relations, weights, control)
    .relation_consensus_CKS(relations, "preorder",
                            weights, control)

### ** .relation_consensus_CKS_transitive

.relation_consensus_CKS_transitive <-
function(relations, weights, control)
    .relation_consensus_CKS(relations, "transitive",
                            weights, control)

## Consensus methods for central relations using PC distance.

### ** .relation_consensus_PC_A

.relation_consensus_PC_A <-
function(relations, weights, control)
    .relation_consensus_PC(relations, "A", weights, control)

### ** .relation_consensus_PC_C

.relation_consensus_PC_C <-
function(relations, weights, control)
    .relation_consensus_PC(relations, "C", weights, control)

### ** .relation_consensus_PC_E

.relation_consensus_PC_E <-
function(relations, weights, control)
{
    k <- control$k
    if(!is.null(k))
        .relation_consensus_PC_E_k(relations, weights, k, control)
    else
        .relation_consensus_PC(relations, "E", weights, control)
}

### ** .relation_consensus_PC_G

## <FIXME>
## Not yet ...
## </FIXME>
.relation_consensus_PC_G <-
function(relations, weights, control)
{
    ## Generalized paired comparison majority.
    ## The consensus problem is to minimize
    ##   \sum_{i,j} L_{ij} x_{ij} + \sum_{i,j:i<j} Q_{ij} x_{ij}x_{ji}
    ## over all binary x_{ij}.
    ## For the x_{ii}, this is archieved upon taking x_{ii} = 0 if
    ## L_{ii} > 0, and x_{ii} = 1 if L_{ii} < 0 (if L_{ii} = 0, both
    ## 0 and 1 can be taken).
    ## For i < j, the values for x_{ij}/x_{ji} 0/0, 1/0, 0/1 and 1/1 are
    ##   0, L_{ij}, L_{ji}, L_{ij} + L_{ji} + Q_{ij}
    ## and the incidences need to be taken to minimize the above.

    delta <- control$delta
    gamma <- control$gamma
    incidences <- lapply(relations, relation_incidence)
    AB <- .make_fit_relation_PC_AB(incidences, weights, delta, gamma)

    B <- AB$B
    n <- nrow(B)
    vd <- diag(B)
    I <- diag(1 - (vd >= 0), n, n)
    ind <- row(B) < col(B)
    pos_ij <- which(ind, arr.ind = TRUE)
    pos_ji <- pos_ij[, c(2L, 1L)]
    vo1 <- B[ind]
    vo2 <- t(B)[ind]
    vo3 <- vo1 + vo2 + AB$A[ind]
    Mo <- which(cbind(0, vo1, vo2, vo3) == pmin(0, vo1, vo2, vo3),
                arr.ind = TRUE)
    Mo <- Mo[order(Mo[, 1L]), , drop = FALSE]
    ## Rows of Mo corresponding to first minima in each i/j row:
    ro <- 1 + c(0, which(diff(Mo[, 1L]) > 0))
    k <- Mo[ro, 1L]
    l <- Mo[ro, 2L] - 1
    ## Values for q are from 0 to 3 corresponding to x_{ij}/x_{ji} pairs
    ## 0/0, 1/0, 0/1 and 1/1: we can get these from q as remainder and
    ## quotient of integer division by 2.
    I[pos_ij[k, , drop = FALSE]] <- l %% 2
    I[pos_ji[k, , drop = FALSE]] <- l %/% 2

    objval <- .relation_consensus_PC_objval(I, incidences, weights,
                                            delta, gamma)
    meta <- list(is_endorelation = TRUE, objval = objval)

    ## If more than one solution is sought, we need to "split" at all
    ## positions where the minimum is not unique.
    nos <- .n_from_control_list(control)

    if(nos > 1L) {
        I <- list(I)
        ## Split diagonal terms when vd == 0.
        for(k in which(vd == 0)) {
            if(length(I) >= nos) break
            I <- do.call(c,
                         lapply(I,
                                function(x) {
                                    y <- x
                                    ## x[k, k] should be 0.
                                    y[k, k] <- 1
                                    list(x, y)
                                }))
        }
        ## Split off-diagonal terms.
        ## Drop the ones used for the first solution.
        Mo <- Mo[-ro, , drop = FALSE]
        for(e in split(Mo, row(Mo))) {
            if(length(I) >= nos) break
            I <- do.call(c,
                         lapply(I,
                                function(x) {
                                    k <- e[1L]
                                    l <- e[2L] - 1
                                    y <- x
                                    y[pos_ij[k, , drop = FALSE]] <- l %% 2
                                    y[pos_ji[k, , drop = FALSE]] <- l %/% 2
                                    list(x, y)
                                }))
        }
        if(length(I) > nos)
            I <- I[seq_len(nos)]
    }
                                    
    .make_consensus_from_incidences(.domain(relations), I, meta)
}

### ** .relation_consensus_PC_L

.relation_consensus_PC_L <-
function(relations, weights, control)
    .relation_consensus_PC(relations, "L", weights, control)

### ** .relation_consensus_PC_M

.relation_consensus_PC_M <-
function(relations, weights, control)
    .relation_consensus_PC(relations, "M", weights, control)

### ** .relation_consensus_PC_O

.relation_consensus_PC_O <-
function(relations, weights, control)
    .relation_consensus_PC(relations, "O", weights, control)

### ** .relation_consensus_PC_S

.relation_consensus_PC_S <-
function(relations, weights, control)
    .relation_consensus_PC(relations, "S", weights, control)

### ** .relation_consensus_PC_T

.relation_consensus_PC_T <-
function(relations, weights, control)
    .relation_consensus_PC(relations, "T", weights, control)

### ** .relation_consensus_PC_W

.relation_consensus_PC_W <-
function(relations, weights, control)
{
    k <- control$k
    if(!is.null(k))
        .relation_consensus_PC_W_k(relations, weights, k, control)
    else
        .relation_consensus_PC(relations, "W", weights, control)
}

### ** .relation_consensus_PC_preorder

.relation_consensus_PC_preorder <-
function(relations, weights, control)
    .relation_consensus_PC(relations, "preorder", weights, control)

### ** .relation_consensus_PC_transitive

.relation_consensus_PC_transitive <-
function(relations, weights, control)
    .relation_consensus_PC(relations, "transitive", weights, control)

### * Relation consensus method registration

## Note that things are simpler here as for CLUE, where we have "typed"
## db's (partition or hierarchy type).

## We currently do without explicit registry getters and setters, which
## in the non-typed case could simplify to the following:
##   get_methods_from_db <-
##   function(db)
##       objects(db)
##   get_method_from_db <-
##   function(db, name)
##       db[[name]]
##   put_method_into_db <-
##   function(db, name, value)
##       db[[name]] <- value
## However, we provide a getter which allows for partial matching.

relation_consensus_methods_db <- new.env()
get_relation_consensus_method <-
function(name)
{
    keys <- objects(relation_consensus_methods_db)
    ind <- pmatch(name, keys)
    if(is.na(ind))
        stop(gettextf("Invalid consensus method '%s'.", name),
             domain = NA)
    relation_consensus_methods_db[[keys[ind]]]
}
set_relation_consensus_method <-
function(name, definition, ...)
{
    ## Note that consensus methods are not necessarily optimization
    ## based (and hence do not necessarily have associated dissimilarity
    ## and exponent).
    value <- c(list(definition = definition), list(...))
    class(value) <- "relation_consensus_method"
    relation_consensus_methods_db[[name]] <- value
}

set_relation_consensus_method("Borda",
                              .relation_consensus_Borda)

set_relation_consensus_method("Copeland",
                              .relation_consensus_Copeland)

set_relation_consensus_method("majority",
                              .relation_consensus_majority)

## Note that constructive methods do not necessarily give central
## relations.

set_relation_consensus_method("CKS/A",
                              .relation_consensus_CKS_A,
                              dissimilarity = "CKS",
                              exponent = 1)
set_relation_consensus_method("CKS/C",
                              .relation_consensus_CKS_C,
                              dissimilarity = "CKS",
                              exponent = 1)
set_relation_consensus_method("CKS/E",
                              .relation_consensus_CKS_E,
                              dissimilarity = "CKS",
                              exponent = 1)
set_relation_consensus_method("CKS/G",
                              .relation_consensus_CKS_G,
                              dissimilarity = "CKS",
                              exponent = 1)
set_relation_consensus_method("CKS/L",
                              .relation_consensus_CKS_L,
                              dissimilarity = "CKS",
                              exponent = 1)
set_relation_consensus_method("CKS/M",
                              .relation_consensus_CKS_M,
                              dissimilarity = "CKS",
                              exponent = 1)
set_relation_consensus_method("CKS/O",
                              .relation_consensus_CKS_O,
                              dissimilarity = "CKS",
                              exponent = 1)
set_relation_consensus_method("CKS/S",
                              .relation_consensus_CKS_S,
                              dissimilarity = "CKS",
                              exponent = 1)
set_relation_consensus_method("CKS/T",
                              .relation_consensus_CKS_T,
                              dissimilarity = "CKS",
                              exponent = 1)
set_relation_consensus_method("CKS/W",
                              .relation_consensus_CKS_W,
                              dissimilarity = "CKS",
                              exponent = 1)
set_relation_consensus_method("CKS/preorder",
                              .relation_consensus_CKS_preorder,
                              dissimilarity = "CKS",
                              exponent = 1)
set_relation_consensus_method("CKS/transitive",
                              .relation_consensus_CKS_transitive,
                              dissimilarity = "CKS",
                              exponent = 1)

set_relation_consensus_method("CS",
                              .relation_consensus_CS,
                              dissimilarity = "CS",
                              exponent = 1)

set_relation_consensus_method("Condorcet",
                              .relation_consensus_Condorcet,
                              dissimilarity = "symdiff",
                              exponent = 1)

set_relation_consensus_method("symdiff/A",
                              .relation_consensus_symdiff_A,
                              dissimilarity = "symdiff",
                              exponent = 1)
set_relation_consensus_method("symdiff/C",
                              .relation_consensus_symdiff_C,
                              dissimilarity = "symdiff",
                              exponent = 1)
set_relation_consensus_method("symdiff/E",
                              .relation_consensus_symdiff_E,
                              dissimilarity = "symdiff",
                              exponent = 1)
set_relation_consensus_method("symdiff/G",
                              .relation_consensus_Condorcet,
                              dissimilarity = "symdiff",
                              exponent = 1)
set_relation_consensus_method("symdiff/L",
                              .relation_consensus_symdiff_L,
                              dissimilarity = "symdiff",
                              exponent = 1)
set_relation_consensus_method("symdiff/M",
                              .relation_consensus_symdiff_M,
                              dissimilarity = "symdiff",
                              exponent = 1)
set_relation_consensus_method("symdiff/O",
                              .relation_consensus_symdiff_O,
                              dissimilarity = "symdiff",
                              exponent = 1)
## <NOTE>
## Keep this for back-compatibility.
set_relation_consensus_method("symdiff/P",
                              .relation_consensus_symdiff_W,
                              dissimilarity = "symdiff",
                              exponent = 1)
## </NOTE>
set_relation_consensus_method("symdiff/S",
                              .relation_consensus_symdiff_S,
                              dissimilarity = "symdiff",
                              exponent = 1)
set_relation_consensus_method("symdiff/T",
                              .relation_consensus_symdiff_T,
                              dissimilarity = "symdiff",
                              exponent = 1)
set_relation_consensus_method("symdiff/W",
                              .relation_consensus_symdiff_W,
                              dissimilarity = "symdiff",
                              exponent = 1)
set_relation_consensus_method("symdiff/preorder",
                              .relation_consensus_symdiff_preorder,
                              dissimilarity = "symdiff",
                              exponent = 1)
set_relation_consensus_method("symdiff/transitive",
                              .relation_consensus_symdiff_transitive,
                              dissimilarity = "symdiff",
                              exponent = 1)
set_relation_consensus_method("SD/A",
                              .relation_consensus_symdiff_A,
                              dissimilarity = "symdiff",
                              exponent = 1)
set_relation_consensus_method("SD/C",
                              .relation_consensus_symdiff_C,
                              dissimilarity = "symdiff",
                              exponent = 1)
set_relation_consensus_method("SD/E",
                              .relation_consensus_symdiff_E,
                              dissimilarity = "symdiff",
                              exponent = 1)
set_relation_consensus_method("SD/G",
                              .relation_consensus_Condorcet,
                              dissimilarity = "symdiff",
                              exponent = 1)
set_relation_consensus_method("SD/L",
                              .relation_consensus_symdiff_L,
                              dissimilarity = "symdiff",
                              exponent = 1)
set_relation_consensus_method("SD/M",
                              .relation_consensus_symdiff_M,
                              dissimilarity = "symdiff",
                              exponent = 1)
set_relation_consensus_method("SD/O",
                              .relation_consensus_symdiff_O,
                              dissimilarity = "symdiff",
                              exponent = 1)
## <NOTE>
## Keep this for back-compatibility.
set_relation_consensus_method("SD/P",
                              .relation_consensus_symdiff_W,
                              dissimilarity = "symdiff",
                              exponent = 1)
## </NOTE>
set_relation_consensus_method("SD/S",
                              .relation_consensus_symdiff_S,
                              dissimilarity = "symdiff",
                              exponent = 1)
set_relation_consensus_method("SD/T",
                              .relation_consensus_symdiff_T,
                              dissimilarity = "symdiff",
                              exponent = 1)
set_relation_consensus_method("SD/W",
                              .relation_consensus_symdiff_W,
                              dissimilarity = "symdiff",
                              exponent = 1)
set_relation_consensus_method("SD/preorder",
                              .relation_consensus_symdiff_preorder,
                              dissimilarity = "symdiff",
                              exponent = 1)
set_relation_consensus_method("SD/transitive",
                              .relation_consensus_symdiff_transitive,
                              dissimilarity = "symdiff",
                              exponent = 1)

set_relation_consensus_method("manhattan",
                              .relation_consensus_manhattan,
                              dissimilarity = "manhattan",
                              exponent = 1)

set_relation_consensus_method("euclidean",
                              .relation_consensus_euclidean,
                              dissimilarity = "euclidean",
                              exponent = 2)
set_relation_consensus_method("euclidean/A",
                              .relation_consensus_euclidean_A,
                              dissimilarity = "euclidean",
                              exponent = 2)
set_relation_consensus_method("euclidean/C",
                              .relation_consensus_euclidean_C,
                              dissimilarity = "euclidean",
                              exponent = 2)
set_relation_consensus_method("euclidean/E",
                              .relation_consensus_euclidean_E,
                              dissimilarity = "euclidean",
                              exponent = 2)
set_relation_consensus_method("euclidean/G",
                              .relation_consensus_Condorcet,
                              dissimilarity = "euclidean",
                              exponent = 2)
set_relation_consensus_method("euclidean/L",
                              .relation_consensus_euclidean_L,
                              dissimilarity = "euclidean",
                              exponent = 2)
set_relation_consensus_method("euclidean/M",
                              .relation_consensus_euclidean_M,
                              dissimilarity = "euclidean",
                              exponent = 2)
set_relation_consensus_method("euclidean/O",
                              .relation_consensus_euclidean_O,
                              dissimilarity = "euclidean",
                              exponent = 2)
set_relation_consensus_method("euclidean/S",
                              .relation_consensus_euclidean_S,
                              dissimilarity = "euclidean",
                              exponent = 2)
set_relation_consensus_method("euclidean/T",
                              .relation_consensus_euclidean_T,
                              dissimilarity = "euclidean",
                              exponent = 2)
set_relation_consensus_method("euclidean/W",
                              .relation_consensus_euclidean_W,
                              dissimilarity = "euclidean",
                              exponent = 2)
set_relation_consensus_method("euclidean/preorder",
                              .relation_consensus_euclidean_preorder,
                              dissimilarity = "euclidean",
                              exponent = 2)
set_relation_consensus_method("euclidean/transitive",
                              .relation_consensus_euclidean_transitive,
                              dissimilarity = "euclidean",
                              exponent = 2)

set_relation_consensus_method("PC/A",
                              .relation_consensus_PC_A,
                              dissimilarity = "PC",
                              exponent = 1)
set_relation_consensus_method("PC/C",
                              .relation_consensus_PC_C,
                              dissimilarity = "PC",
                              exponent = 1)
set_relation_consensus_method("PC/E",
                              .relation_consensus_PC_E,
                              dissimilarity = "PC",
                              exponent = 1)
set_relation_consensus_method("PC/G",
                              .relation_consensus_PC_G,
                              dissimilarity = "PC",
                              exponent = 1)
set_relation_consensus_method("PC/L",
                              .relation_consensus_PC_L,
                              dissimilarity = "PC",
                              exponent = 1)
set_relation_consensus_method("PC/M",
                              .relation_consensus_PC_M,
                              dissimilarity = "PC",
                              exponent = 1)
set_relation_consensus_method("PC/O",
                              .relation_consensus_PC_O,
                              dissimilarity = "PC",
                              exponent = 1)
set_relation_consensus_method("PC/S",
                              .relation_consensus_PC_S,
                              dissimilarity = "PC",
                              exponent = 1)
set_relation_consensus_method("PC/T",
                              .relation_consensus_PC_T,
                              dissimilarity = "PC",
                              exponent = 1)
set_relation_consensus_method("PC/W",
                              .relation_consensus_PC_W,
                              dissimilarity = "PC",
                              exponent = 1)
set_relation_consensus_method("PC/preorder",
                              .relation_consensus_PC_preorder,
                              dissimilarity = "PC",
                              exponent = 1)
set_relation_consensus_method("PC/transitive",
                              .relation_consensus_PC_transitive,
                              dissimilarity = "PC",
                              exponent = 1)


### * Relation consensus workers

### ** .relation_consensus_symdiff

.relation_consensus_symdiff <-
function(relations, family, weights, control, euclidean = FALSE)
{
    if(!.is_ensemble_of_endorelations(relations))
        stop("Need an ensemble of endorelations.")
    if(!euclidean && !.is_ensemble_of_crisp_relations(relations))
        stop("Need an ensemble of crisp relations.")

    incidences <- lapply(relations, relation_incidence)
    B <- .make_fit_relation_symdiff_B(incidences, weights)
    I <- fit_relation_LP(B, family, control)
    objval <- if(euclidean)
        .relation_consensus_euclidean_objval(I, incidences, weights)
    else
        .relation_consensus_symdiff_objval(I, incidences, weights)
    meta <- c(.relation_meta_db[[family]], list(objval = objval))

    .make_consensus_from_incidences(.domain(relations), I, meta)
}

### ** .relation_consensus_symdiff_E_k

.relation_consensus_symdiff_E_k <-
function(relations, weights, k, control, euclidean = FALSE)
{
    if(!.is_ensemble_of_endorelations(relations))
        stop("Need an ensemble of endorelations.")
    if(!euclidean && !.is_ensemble_of_crisp_relations(relations))
        stop("Need an ensemble of crisp relations.")

    incidences <- lapply(relations, relation_incidence)
    B <- .make_fit_relation_symdiff_B(incidences, weights)
    I <- fit_relation_LP_E_k(B, k, control)
    objval <- if(euclidean)
        .relation_consensus_euclidean_objval(I, incidences, weights)
    else
        .relation_consensus_symdiff_objval(I, incidences, weights)
    meta <- c(.relation_meta_db[["E"]], list(objval = objval))

    .make_consensus_from_incidences(.domain(relations), I, meta)
}

### ** .relation_consensus_symdiff_W_k

.relation_consensus_symdiff_W_k <-
function(relations, weights, k, control, euclidean = FALSE)
{
    if(!.is_ensemble_of_endorelations(relations))
        stop("Need an ensemble of endorelations.")
    if(!euclidean && !.is_ensemble_of_crisp_relations(relations))
        stop("Need an ensemble of crisp relations.")

    incidences <- lapply(relations, relation_incidence)
    B <- .make_fit_relation_symdiff_B(incidences, weights)
    I <- fit_relation_LP_W_k(B, k, control)
    objval <- if(euclidean)
        .relation_consensus_euclidean_objval(I, incidences, weights)
    else
        .relation_consensus_symdiff_objval(I, incidences, weights)
    meta <- c(.relation_meta_db[["W"]], list(objval = objval))

    .make_consensus_from_incidences(.domain(relations), I, meta)
}

## ** .relation_consensus_CKS

.relation_consensus_CKS <-
function(relations, family, weights, control)
{
    if(!.is_ensemble_of_endorelations(relations))
        stop("Need an ensemble of endorelations.")
    if(!.is_ensemble_of_crisp_relations(relations))
        stop("Need an ensemble of crisp relations.")

    incidences <- lapply(relations, relation_incidence)
    AB <- .make_fit_relation_CKS_AB(incidences, weights)
    if(identical(control$simplify, FALSE) ||
       family %in% c("preorder", "transitive")) {
        ## Solve via QP formulation if requested explicitly or
        ## necessary.
        ## Mostly useful for testing purposes.
        I <- fit_relation_QP(AB$A, AB$B, family, control)
    } else {
        if(family %in% c("A", "L", "O", "T")) {
            B <- AB$B
        } else if(family %in% c("E", "S")) {
            B <- AB$B + AB$A
        } else if(family %in% c("C", "M", "W")) {
            B <- AB$B + AB$A + t(AB$A)
        } else {
            stop("Not implemented.")
        }
        I <- fit_relation_LP(B, family, control)
    }
    objval <- .relation_consensus_CKS_objval(I, incidences, weights)
    meta <- c(.relation_meta_db[[family]], list(objval = objval))

    .make_consensus_from_incidences(.domain(relations), I, meta)
}

## ** .relation_consensus_CKS_E_k

.relation_consensus_CKS_E_k <-
function(relations, weights, k, control)
{
    if(!.is_ensemble_of_endorelations(relations))
        stop("Need an ensemble of endorelations.")
    if(!.is_ensemble_of_crisp_relations(relations))
        stop("Need an ensemble of crisp relations.")

    incidences <- lapply(relations, relation_incidence)
    ## The coefficients we need are B + A.
    ## By symmetry, we can also use M = (B + B' + A + A') / 2, and
    ## diagonal terms are irrelevant.
    ## Now for i \ne j,
    ##   [B + B']_{ij} = 2 q_{ij} - p_{ij} - p_{ji}
    ##   [A + A']_{ij} = - q_{ij} + p_{ij} + p_{ji}
    ## Hence we can simply take Q for the coefficients.
    ## Alternatively, we can proceed as for general PC_E_k.
    
    B <- .make_fit_relation_CKS_Q(incidences, weights)
    I <- fit_relation_LP_E_k(B, k, control)
    objval <- .relation_consensus_CKS_objval(I, incidences, weights)
    meta <- c(.relation_meta_db[["E"]], list(objval = objval))

    .make_consensus_from_incidences(.domain(relations), I, meta)
}

## ** .relation_consensus_CKS_W_k

.relation_consensus_CKS_W_k <-
function(relations, weights, k, control)
{
    if(!.is_ensemble_of_endorelations(relations))
        stop("Need an ensemble of endorelations.")
    if(!.is_ensemble_of_crisp_relations(relations))
        stop("Need an ensemble of crisp relations.")

    incidences <- lapply(relations, relation_incidence)
    AB <- .make_fit_relation_CKS_AB(incidences, weights)
    B <- AB$B + AB$A + t(AB$A)
    I <- fit_relation_LP_W_k(B, k, control)    
    objval <- .relation_consensus_CKS_objval(I, incidences, weights)
    meta <- c(.relation_meta_db[["W"]], list(objval = objval))

    .make_consensus_from_incidences(.domain(relations), I, meta)
}

## ** .relation_consensus_PC

.relation_consensus_PC <-
function(relations, family, weights, control)
{
    if(!.is_ensemble_of_endorelations(relations))
        stop("Need an ensemble of endorelations.")
    if(!.is_ensemble_of_crisp_relations(relations))
        stop("Need an ensemble of crisp relations.")

    delta <- control$delta
    gamma <- control$gamma
    incidences <- lapply(relations, relation_incidence)
    AB <- .make_fit_relation_PC_AB(incidences, weights, delta, gamma)
    if(identical(control$simplify, FALSE) ||
       family %in% c("preorder", "transitive")) {
        ## Solve via QP formulation if requested explicitly or
        ## necessary.
        ## Mostly useful for testing purposes.
        I <- fit_relation_QP(AB$A, AB$B, family, control)
    } else {
        if(family %in% c("A", "L", "O", "T")) {
            B <- AB$B
        } else if(family %in% c("E", "S")) {
            B <- AB$B + AB$A
        } else if(family %in% c("C", "M", "W")) {
            B <- AB$B + AB$A + t(AB$A)
        } else {
            stop("Not implemented.")
        }
        I <- fit_relation_LP(B, family, control)
    }
    objval <- .relation_consensus_PC_objval(I, incidences, weights,
                                            delta, gamma)
    meta <- c(.relation_meta_db[[family]], list(objval = objval))

    .make_consensus_from_incidences(.domain(relations), I, meta)
}

## ** .relation_consensus_PC_E_k

.relation_consensus_PC_E_k <-
function(relations, weights, k, control)
{
    if(!.is_ensemble_of_endorelations(relations))
        stop("Need an ensemble of endorelations.")
    if(!.is_ensemble_of_crisp_relations(relations))
        stop("Need an ensemble of crisp relations.")

    delta <- control$delta
    gamma <- control$gamma
    incidences <- lapply(relations, relation_incidence)
    AB <- .make_fit_relation_PC_AB(incidences, weights, delta, gamma)
    B <- AB$B + AB$A
    I <- fit_relation_LP_E_k(B, k, control)
    objval <- .relation_consensus_PC_objval(I, incidences, weights,
                                            delta, gamma)
    meta <- c(.relation_meta_db[["E"]], list(objval = objval))

    .make_consensus_from_incidences(.domain(relations), I, meta)
}

## ** .relation_consensus_PC_W_k

.relation_consensus_PC_W_k <-
function(relations, weights, k, control)
{
    if(!.is_ensemble_of_endorelations(relations))
        stop("Need an ensemble of endorelations.")
    if(!.is_ensemble_of_crisp_relations(relations))
        stop("Need an ensemble of crisp relations.")

    delta <- control$delta
    gamma <- control$gamma
    incidences <- lapply(relations, relation_incidence)
    AB <- .make_fit_relation_PC_AB(incidences, weights, delta, gamma)
    B <- AB$B + AB$A + t(AB$A)
    I <- fit_relation_LP_W_k(B, k, control)
    objval <- .relation_consensus_PC_objval(I, incidences, weights,
                                            delta, gamma)
    meta <- c(.relation_meta_db[["W"]], list(objval = objval))

    .make_consensus_from_incidences(.domain(relations), I, meta)
}

### * Utilities

### ** .make_fit_relation_symdiff_B

.make_fit_relation_symdiff_B <-
function(incidences, weights, flip = FALSE)
{
    ## Compute the array
    ##   \sum_b w_b (1 - 2 incidence(b))
    ## or
    ##   \sum_b w_b (2 incidence(b) - 1)
    ## if flip is true.

    w <- rep_len(weights, length(incidences))
    if(is.relation_ensemble(incidences))
        incidences <- lapply(incidences, relation_incidence)
    M <- .weighted_sum_of_arrays(incidences, w)
    if(flip)
        2 * M - sum(w)
    else
        sum(w) - 2 * M
}

### ** .make_fit_relation_CKS_P

.make_fit_relation_CKS_P <-
function(incidences, weights)
    .make_fit_relation_symdiff_B(lapply(incidences,
                                        function(I) pmin(I, 1 - t(I))),
                                 weights,
                                 TRUE)

### ** .make_fit_relation_CKS_Q

.make_fit_relation_CKS_Q <-
function(incidences, weights)
    .make_fit_relation_symdiff_B(lapply(incidences,
                                        function(I) 1 - pmax(I, t(I))),
                                 weights,
                                 TRUE)

### ** .make_fit_relation_CKS_AB

.make_fit_relation_CKS_AB <-
function(incidences, weights, P = NULL, Q = NULL)
{
    if(is.null(P))
        P <- .make_fit_relation_CKS_P(incidences, weights)
    if(is.null(Q))
        Q <- .make_fit_relation_CKS_Q(incidences, weights)
    B <- Q - P
    diag(B) <- diag(Q) / 2
    A <- P + t(P) - Q
    A[row(A) >= col(A)] <- 0
    list(A = A, B = B)
}
        
### ** .make_fit_relation_PC_AB

.make_fit_relation_PC_AB <-
function(incidences, weights, delta, gamma)
{
    w <- rep_len(weights, length(incidences))
    if(is.relation_ensemble(incidences)) 
        incidences <- lapply(incidences, relation_incidence)
    D <- .relation_dissimilarity_PC_Delta(delta)
    M <- .relation_dissimilarity_PC_M(D)
    ABC <- lapply(incidences, .relation_dissimilarity_PC_ABC,
                  delta, gamma, D, M)
    A <- .weighted_sum_of_arrays(lapply(ABC, `[[`, "A"), w)
    B <- .weighted_sum_of_arrays(lapply(ABC, `[[`, "B"), w)
    list(A = A, B = B)
}

### ** .make_consensus_from_incidences

.make_consensus_from_incidences <-
function(D, I, meta = NULL)
{
    if(is.list(I)) {
        ## In case *all* consensus solutions were sought ...
        relations <-
            lapply(I,
                   function(e)
                   .make_relation_from_domain_and_incidence(D, e, meta))
        relation_ensemble(list = relations)
    }
    else
        .make_relation_from_domain_and_incidence(D, I, meta)
}

## Utilities for computing the values of the objective functions for the
## optimization-based consensus methods.

### ** .relation_consensus_symdiff_objval

.relation_consensus_symdiff_objval <-
function(I, incidences, weights)
{
    ## Be nice.
    if(is.relation_ensemble(incidences))
        incidences <- lapply(incidences, relation_incidence)
    if(is.list(I)) I <- I[[1L]]

    sum(weights *
        sapply(incidences, .incidence_dissimilarity_symdiff, I))
}

### ** .relation_consensus_manhattan_objval

.relation_consensus_manhattan_objval <-
function(I, incidences, weights)
{
    ## Be nice.
    if(is.relation_ensemble(incidences))
        incidences <- lapply(incidences, relation_incidence)
    if(is.list(I)) I <- I[[1L]]

    sum(weights *
        sapply(incidences, .incidence_dissimilarity_manhattan, I))
}

### ** .relation_consensus_euclidean_objval

.relation_consensus_euclidean_objval <-
function(I, incidences, weights)
{
    ## Be nice.
    if(is.relation_ensemble(incidences))
        incidences <- lapply(incidences, relation_incidence)
    if(is.list(I)) I <- I[[1L]]

    sum(weights *
        sapply(incidences, .incidence_dissimilarity_euclidean, I) ^ 2)
}

### ** .relation_consensus_CS_objval

.relation_consensus_CS_objval <-
function(I, incidences, weights)
{
    ## Be nice.
    if(is.relation_ensemble(incidences))
        incidences <- lapply(incidences, relation_incidence)
    if(is.list(I)) I <- I[[1L]]

    sum(weights *
        sapply(incidences, .incidence_dissimilarity_CS, I))
}

### ** .relation_consensus_CKS_objval

.relation_consensus_CKS_objval <-
function(I, incidences, weights)
{
    ## Be nice.
    if(is.relation_ensemble(incidences))
        incidences <- lapply(incidences, relation_incidence)
    if(is.list(I)) I <- I[[1L]]

    sum(weights *
        sapply(incidences, .incidence_dissimilarity_CKS, I))
}

## ** .relation_consensus_PC_objval

.relation_consensus_PC_objval <-
function(I, incidences, weights, delta, gamma)
{
    ## Be nice.
    if(is.relation_ensemble(incidences))
        incidences <- lapply(incidences, relation_incidence)
    if(is.list(I)) I <- I[[1L]]

    sum(weights *
        sapply(incidences,
               .incidence_dissimilarity_PC, I, delta, gamma))
}



### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "### [*]+" ***
### End: ***
