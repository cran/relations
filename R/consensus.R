### Consensus relations.

### * relation_consensus

relation_consensus <-
function(x, method = NULL, weights = 1, control = list())
{
    relations <- as.relation_ensemble(x)

    if(!length(relations))
        stop("Cannot compute consensus of empty ensemble.")

    weights <- rep(weights, length = length(relations))
    if(any(weights < 0))
        stop("Argument 'weights' has negative elements.")
    if(!any(weights > 0))
        stop("Argument 'weights' has no positive elements.")

    ## Eventually, more sophisticated handling of 'method' as in CLUE's
    ## cl_consensus().  For the time being, have an internal look up
    ## table here.
    known_methods <- c(Borda = ".relation_consensus_Borda",
                       Condorcet = ".relation_consensus_Condorcet",
                       "SD/E" = ".relation_consensus_symdiff_E",
                       "SD/L" = ".relation_consensus_symdiff_L",
                       "SD/O" = ".relation_consensus_symdiff_O",
                       "SD/P" = ".relation_consensus_symdiff_P",
                       "SD/T" = ".relation_consensus_symdiff_T"
                       )

    if(is.character(method)) {
        ## Hopefully of length one, add some tests eventually ...
        if(is.na(ind <- pmatch(method, names(known_methods))))
            stop(gettextf("Method '%s' is not a valid consensus method.",
                          method))
        method <- get(known_methods[ind])
    }
    else if(!is.function(method))
        stop("Argument 'method' must be a function or character string.")
    
    method(relations, weights, control)
}

### * Relation consensus "methods"

### ** .relation_consensus_Borda

.relation_consensus_Borda <-
function(relations, weights, control)
{
    ## Several sanity checks would be necessary here.
    ## Check whether all relations are in fact complete preferences (or
    ## whatever is really necessary).
    ## We currently cannot handle non-identical weights.
    ## Etc etc ...

    ## Note that we currently do things the other way round.

    if(!.is_ensemble_of_endorelations(relations))
        stop("Need an ensemble of endorelations.")

    ## First, get the incidences.
    incidences <- lapply(relations, relation_incidence)
    ## From this, get the Borda/Kendall scores.
    scores <- lapply(incidences, rowSums)
    ## Multiply by the weights and compute the total scores.
    scores <- rowSums(mapply("*", scores, weights))

    ind <- seq_along(scores)
    out <- integer(length = length(ind))
    out[order(scores, decreasing = TRUE)] <- ind
    names(out) <- rownames(incidences[[1L]])

    I <- outer(out, out, "<=")
    meta <- list(is_endorelation = TRUE,
                 is_complete = TRUE,
                 is_reflexive = TRUE,
                 is_antisymmetric = !any(duplicated(out)),
                 is_transitive = TRUE,
                 scores = scores)
    .make_relation_from_domain_and_incidence(.domain(relations), I, meta)
}

### ** .relation_consensus_Condorcet

.relation_consensus_Condorcet <-
function(relations, weights, control)
{
    I <- .make_fit_relation_symdiff_M(relations, weights) > 0
    .make_relation_from_domain_and_incidence(.domain(relations), I)
}

### ** .relation_consensus_symdiff_E

.relation_consensus_symdiff_E <-
function(relations, weights, control)
{
    k <- control$k
    I <- if(!is.null(k))
        .relation_consensus_symdiff_E_k(relations, weights, k, control)
    else
        .relation_consensus_symdiff(relations, "E", weights, control)
}

### ** .relation_consensus_symdiff_L

.relation_consensus_symdiff_L <-
function(relations, weights, control)
    .relation_consensus_symdiff(relations, "L", weights, control)

### ** .relation_consensus_symdiff_O

.relation_consensus_symdiff_O <-
function(relations, weights, control)
    .relation_consensus_symdiff(relations, "O", weights, control)

### ** .relation_consensus_symdiff_P

.relation_consensus_symdiff_P <-
function(relations, weights, control)
{
    k <- control$k
    I <- if(!is.null(k))
        .relation_consensus_symdiff_P_k(relations, weights, k, control)
    else
        .relation_consensus_symdiff(relations, "P", weights, control)
}

### ** .relation_consensus_symdiff_T

.relation_consensus_symdiff_T <-
function(relations, weights, control)
    .relation_consensus_symdiff(relations, "T", weights, control)

### * Relations consensus workers

### ** .relation_consensus_symdiff

.relation_consensus_symdiff <-
function(relations, family, weights, control)
{
    if(!.is_ensemble_of_endorelations(relations))
        stop("Need an ensemble of endorelations.")

    M <- .make_fit_relation_symdiff_M(relations, weights)
    I <- fit_relation_symdiff(M, family, control)
    meta <- .relation_meta_db[[family]]
    
    .make_consensus_from_incidences(.domain(relations), I, meta)
}

### ** .relation_consensus_symdiff_E_k

.relation_consensus_symdiff_E_k <-
function(relations, weights, k, control)
{
    if(!.is_ensemble_of_endorelations(relations))
        stop("Need an ensemble of endorelations.")

    M <- .make_fit_relation_symdiff_M(relations, weights)
    I <- fit_relation_symdiff_E_k(M, k, control)
    meta <- .relation_meta_db[["E"]]
    
    .make_consensus_from_incidences(.domain(relations), I, meta)
}

### ** .relation_consensus_symdiff_P_k

.relation_consensus_symdiff_P_k <-
function(relations, weights, k, control)
{
    if(!.is_ensemble_of_endorelations(relations))
        stop("Need an ensemble of endorelations.")

    M <- .make_fit_relation_symdiff_M(relations, weights)
    I <- fit_relation_symdiff_P_k(M, k, control)
    meta <- .relation_meta_db[["P"]]
    
    .make_consensus_from_incidences(.domain(relations), I, meta)
}

### * Utilities

### ** .is_ensemble_of_endorelations

.is_ensemble_of_endorelations <-
function(x)
{
    ## Check whether we have an ensemble of endorelations (assuming that
    ## ensembles are known to have identical identical domains).
    relation_is_endorelation(x[[1L]])
}

### ** .make_fit_relation_symdiff_M

.make_fit_relation_symdiff_M <-
function(relations, weights)
{
    ## Compute the array
    ##   \sum_b w_b (2 incidence(b) - 1)
    ## used in fit_relation_symdiff() and also for the Condorcet
    ## consensus solution.
    
    w <- rep(weights, length = length(relations))
    incidences <- lapply(relations, relation_incidence)
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

### ** .non_transitivity

.non_transitivity <-
function(x)
{
    x <- relation_incidence(x)
    ## <NOTE>
    ## Here, we really want the number of triples for which the
    ## transitivity constraint is violated.
    ## </NOTE>
    n <- nrow(x)
    sum(sapply(seq_len(n),
               function(j) outer(x[, j], x[j, ], "+") - x) > 1)
}    

### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "### [*]+" ***
### End: ***
