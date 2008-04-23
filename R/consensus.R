### Consensus relations.

### * relation_consensus

relation_consensus <-
function(x, method = NULL, weights = 1, control = list(), ...)
{
    dots <- list(...)
    control[names(dots)] <- dots
    relations <- as.relation_ensemble(x)

    if(!length(relations))
        stop("Cannot compute consensus of empty ensemble.")

    weights <- rep(weights, length.out = length(relations))
    if(any(weights < 0))
        stop("Argument 'weights' has negative elements.")
    if(!any(weights > 0))
        stop("Argument 'weights' has no positive elements.")

    if(!is.function(method)) {
        if(!inherits(method, "relation_consensus_method")) {
            ## Get the method definition from the registry.
            if(!is.character(method) || (length(method) != 1))
                stop("Invalid 'method' argument.")
            entry <- get_relation_consensus_method(method)
            if(is.null(entry))
                stop(gettextf("Method '%s' is not a valid consensus method.",
                              method))
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
                              relation_scores, "differential")

### ** .relation_consensus_score

.relation_consensus_score <-
function(relations, weights, control, FUN, ...)
{
    ## Several sanity checks could be done here.
    ## In particular, check whether all relations are in fact complete
    ## preferences (or whatever is really necessary).

    if(!.is_ensemble_of_endorelations(relations))
        stop("Need an ensemble of endorelations.")

    ## Get the scores.
    scores <- lapply(relations, FUN, ...)
    ## Multiply by the weights and compute the total scores.
    scores <- rowSums(mapply("*", scores, weights))

    ind <- seq_along(scores)
    out <- integer(length = length(ind))
    out[order(scores)] <- ind

    I <- outer(out, out, "<=")
    meta <- list(is_endorelation = TRUE,
                 is_complete = TRUE,
                 is_reflexive = TRUE,
                 is_antisymmetric = !any(duplicated(out)),
                 is_transitive = TRUE,
                 scores = scores)
    .make_relation_from_domain_and_incidence(.domain(relations), I, meta)
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

    all <- control$all
    if(is.null(all)) all <- FALSE

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
    I <- if(all)
        lapply(.find_all_LSAP_solutions(C), .compare)
    else
        .compare(clue::solve_LSAP(C))
    objval <- .relation_consensus_CS_objval(I, incidences, weights)
    meta <- c(.relation_meta_db[["L"]], list(objval = objval))

    .make_consensus_from_incidences(.domain(relations), I, meta)
}

## Consensus methods for central relations using symdiff/Manhattan or
## Euclidean distance.  Note that
##
## * We have restricted symdiff fitters only for crisp ensembles.
##   Fitters for fuzzy ensembles could be added by the usual means of
##   turning an l_1 problem with linear/integer constraints into a
##   MILP (see e.g. CLUE).
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

    incidences <- lapply(relations, relation_incidence)
    M <- .make_fit_relation_symdiff_M(incidences, weights)
    diag(M) <- 1                        # Ensure reflexivity.
    I <- M >= 0                         # We do not break ties (>=).
    objval <- .relation_consensus_symdiff_objval(I, incidences, weights)
    meta <- list(is_endorelation = TRUE,
                 is_complete = TRUE,
                 is_reflexive = TRUE,
                 is_antisymmetric = all(M != 0),
                 objval = objval)
    ## (We do not know about transitivity without computing ...) 
    .make_relation_from_domain_and_incidence(.domain(relations), I, meta)
}

### ** .relation_consensus_symdiff_A

.relation_consensus_symdiff_A <-
function(relations, weights, control)
    .relation_consensus_symdiff(relations, "A", weights, control)

### ** .relation_consensus_symdiff_C

.relation_consensus_symdiff_C <-
function(relations, weights, control)
    .relation_consensus_symdiff(relations, "C", weights, control)

### ** .relation_consensus_symdiff_R

.relation_consensus_symdiff_R <-
function(relations, weights, control)
    .relation_consensus_symdiff(relations, "R", weights, control)

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

### ** .relation_consensus_symdiff_P

.relation_consensus_symdiff_P <-
function(relations, weights, control)
{
    k <- control$k
    if(!is.null(k))
        .relation_consensus_symdiff_P_k(relations, weights, k, control)
    else
        .relation_consensus_symdiff(relations, "P", weights, control)
}

### ** .relation_consensus_symdiff_S

.relation_consensus_symdiff_S <-
function(relations, weights, control)
    .relation_consensus_symdiff(relations, "S", weights, control)

### ** .relation_consensus_symdiff_T

.relation_consensus_symdiff_T <-
function(relations, weights, control)
    .relation_consensus_symdiff(relations, "T", weights, control)

### ** .relation_consensus_manhattan

## Symdiff/Manhattan valued consensus relation for an ensemble of
## valued relations.

.relation_consensus_manhattan <-
function(relations, weights, control)
{
    incidences <- lapply(relations, relation_incidence)
    weights <- rep(weights, length.out = length(incidences))
    ## Incidences of the consensus relation are the weighted medians.
    I <- array(apply(do.call("cbind", lapply(incidences, c)),
                     1L, clue:::weighted_median, weights),
               dim = .size(relations))
    meta <-
        list(objval =
             .relation_consensus_symdiff_objval(I, incidences, weights))
    .make_relation_from_domain_and_incidence(.domain(relations), I, meta)
}

### ** .relation_consensus_euclidean

## Euclidean valued consensus relation for an ensemble of valued
## relations.

.relation_consensus_euclidean <-
function(relations, weights, control)
{
    weights <- rep(weights, length.out = length(relations))
    weights <- weights / sum(weights)
    incidences <- lapply(relations, relation_incidence)
    ## Incidences of the consensus relation are the weighted means.
    I <- array(rowSums(mapply("*", incidences, weights)),
               dim = dim(incidences[[1L]]))
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

### ** .relation_consensus_euclidean_P

.relation_consensus_euclidean_P <-
function(relations, weights, control)
{
    k <- control$k
    if(!is.null(k))
        .relation_consensus_symdiff_P_k(relations, weights, k, control, TRUE)
    else
        .relation_consensus_symdiff(relations, "P", weights, control, TRUE)
}

### ** .relation_consensus_euclidean_S

.relation_consensus_euclidean_S <-
function(relations, weights, control)
    .relation_consensus_symdiff(relations, "S", weights, control, TRUE)

### ** .relation_consensus_euclidean_T

.relation_consensus_euclidean_T <-
function(relations, weights, control)
    .relation_consensus_symdiff(relations, "T", weights, control, TRUE)


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
##   {
##       db[[name]]
##   }
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
        stop(gettextf("Invalid consensus method '%s'.", name))
    relation_consensus_methods_db[[keys[ind]]]
}
set_relation_consensus_method <-
function(name, definition, ...)
{
    ## Note that consensus methods are not necessarily optimization
    ## based (and hence do not necessarily have associated dissimilarity
    ## and exponent).
    value <- structure(c(list(definition = definition), list(...)),
                       class = "relation_consensus_method")
    relation_consensus_methods_db[[name]] <- value
}


set_relation_consensus_method("Borda",
                              .relation_consensus_Borda)
set_relation_consensus_method("Copeland",
                              .relation_consensus_Copeland)
## Note that constructive methods do not necessarily give central
## relations.
set_relation_consensus_method("CS",
                              .relation_consensus_CS,
                              dissimilarity = "CS",
                              exponent = 1)
set_relation_consensus_method("Condorcet",
                              .relation_consensus_Condorcet,
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
set_relation_consensus_method("SD/P",
                              .relation_consensus_symdiff_P,
                              dissimilarity = "symdiff",
                              exponent = 1)
set_relation_consensus_method("SD/S",
                              .relation_consensus_symdiff_S,
                              dissimilarity = "symdiff",
                              exponent = 1)
set_relation_consensus_method("SD/T",
                              .relation_consensus_symdiff_T,
                              dissimilarity = "symdiff",
                              exponent = 1)
set_relation_consensus_method("SD/R",
                              .relation_consensus_symdiff_R,
                              dissimilarity = "symdiff",
                              exponent = 1)
set_relation_consensus_method("manhattan",
                              .relation_consensus_manhattan,
                              dissimilarity = "symdiff",
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
set_relation_consensus_method("euclidean/P",
                              .relation_consensus_euclidean_P,
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
    M <- .make_fit_relation_symdiff_M(incidences, weights)
    I <- fit_relation_symdiff(M, family, control)
    objval <- .relation_consensus_symdiff_objval(I, incidences, weights)
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
    M <- .make_fit_relation_symdiff_M(incidences, weights)
    I <- fit_relation_symdiff_E_k(M, k, control)
    objval <- .relation_consensus_symdiff_objval(I, incidences, weights)
    meta <- c(.relation_meta_db[["E"]], list(objval = objval))

    .make_consensus_from_incidences(.domain(relations), I, meta)
}

### ** .relation_consensus_symdiff_P_k

.relation_consensus_symdiff_P_k <-
function(relations, weights, k, control, euclidean = FALSE)
{
    if(!.is_ensemble_of_endorelations(relations))
        stop("Need an ensemble of endorelations.")
    if(!euclidean && !.is_ensemble_of_crisp_relations(relations))
        stop("Need an ensemble of crisp relations.")

    incidences <- lapply(relations, relation_incidence)
    M <- .make_fit_relation_symdiff_M(incidences, weights)
    I <- fit_relation_symdiff_P_k(M, k, control)
    objval <- .relation_consensus_symdiff_objval(I, incidences, weights)
    meta <- c(.relation_meta_db[["P"]], list(objval = objval))

    .make_consensus_from_incidences(.domain(relations), I, meta)
}

### * Utilities

### ** .make_fit_relation_symdiff_M

.make_fit_relation_symdiff_M <-
function(incidences, weights)
{
    ## Compute the array
    ##   \sum_b w_b (2 incidence(b) - 1)
    ## used in fit_relation_symdiff() and also for the Condorcet
    ## consensus solution.

    w <- rep(weights, length.out = length(incidences))
    if(is.relation_ensemble(incidences))
        incidences <- lapply(incidences, relation_incidence)
    M <- array(rowSums(mapply("*", incidences, w)),
               dim = dim(incidences[[1L]]),
               dimnames = dimnames(incidences[[1L]]))
    2 * M - sum(w)
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

    sum(weights * sapply(incidences, .incidence_dissimilarity_CS, I))
}

### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "### [*]+" ***
### End: ***
