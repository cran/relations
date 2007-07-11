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
function(relation, weights, control)
    .relation_consensus_Borda_like(relation, weights, control,
                                   function(x) colSums(x))

### ** .relation_consensus_Copeland

## Copeland method
.relation_consensus_Copeland <-
function(relation, weights, control)
    .relation_consensus_Borda_like(relation, weights, control,
                                   function(x) colSums(x) - rowSums(x))

### ** .relation_consensus_Borda_like

.relation_consensus_Borda_like <-
function(relations, weights, control, FUN)
{
    ## Several sanity checks could be done here.
    ## In particular, check whether all relations are in fact complete
    ## preferences (or whatever is really necessary).

    if(!.is_ensemble_of_endorelations(relations))
        stop("Need an ensemble of endorelations.")

    ## First, get the incidences.
    incidences <- lapply(relations, relation_incidence)
    ## From this, get scores.
    scores <- lapply(incidences, FUN)
    ## Multiply by the weights and compute the total scores.
    scores <- rowSums(mapply("*", scores, weights))

    ind <- seq_along(scores)
    out <- integer(length = length(ind))
    out[order(scores)] <- ind
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

### ** .relation_consensus_symdiff_C

.relation_consensus_symdiff_C <-
function(relations, weights, control)
    .relation_consensus_symdiff(relations, "C", weights, control)

### ** .relation_consensus_symdiff_A

.relation_consensus_symdiff_A <-
function(relations, weights, control)
    .relation_consensus_symdiff(relations, "A", weights, control)

### ** .relation_consensus_symdiff_S

.relation_consensus_symdiff_S <-
function(relations, weights, control)
    .relation_consensus_symdiff(relations, "S", weights, control)

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
set_relation_consensus_method("Condorcet",
                              .relation_consensus_Condorcet,
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
set_relation_consensus_method("SD/O",
                              .relation_consensus_symdiff_O,
                              dissimilarity = "symdiff",
                              exponent = 1)
set_relation_consensus_method("SD/P",
                              .relation_consensus_symdiff_P,
                              dissimilarity = "symdiff",
                              exponent = 1)
set_relation_consensus_method("SD/T",
                              .relation_consensus_symdiff_T,
                              dissimilarity = "symdiff",
                              exponent = 1)
set_relation_consensus_method("SD/C",
                              .relation_consensus_symdiff_C,
                              dissimilarity = "symdiff",
                              exponent = 1)
set_relation_consensus_method("SD/A",
                              .relation_consensus_symdiff_A,
                              dissimilarity = "symdiff",
                              exponent = 1)
set_relation_consensus_method("SD/S",
                              .relation_consensus_symdiff_S,
                              dissimilarity = "symdiff",
                              exponent = 1)

### * Relation consensus workers

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

    w <- rep(weights, length.out = length(relations))
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

### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "### [*]+" ***
### End: ***
