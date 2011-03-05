### methods for closure and reduction

closure.relation <-
function(x, operation = c("transitive", "reflexive"), ...)
{
    operation <- match.arg(operation)
    if (operation == "transitive")
        transitive_closure(x)
    else
        reflexive_closure(x)
}

reduction.relation <-
function(x, operation = c("transitive", "reflexive"), ...)
{
    operation <- match.arg(operation)
    if (operation == "transitive")
        transitive_reduction(x)
    else
        reflexive_reduction(x)
}

### * transitive_reduction

## !!The following code is erroneous:!!
##
## transitive_reduction <-
## function(x)
## {
##     if(!(is.relation(x) && relation_is_endorelation(x)))
##         stop("Argument 'x' must be an endorelation.")
##     if(!relation_is_crisp(x, na.rm = TRUE))
##         stop("Argument 'x' must be a crisp relation.")
##     I <- relation_incidence(x)
##     diag_hold <- diag(I)
##     diag(I) <- 0
##     is_transitive <- FALSE
##     for (i in seq_len(ncol(I))) {
##         tmp <- outer(I[,i], I[i,], .T.)
##         if (any(is.na(tmp))) is_transitive <- NA
##         I <- .T.(.S.(I, is.na(I)), .N.(tmp))
##     }
##     meta <- list(is_endorelation = TRUE,
##                  is_transitive = is_transitive)
##     diag(I) <- diag_hold
##     .make_relation_from_domain_and_incidence(.domain(x), I, meta)
## }

## <FIXME>
## The following code only works for acyclic relations!
## </FIXME>

transitive_reduction <-
function(x)
{
    if(!(is.relation(x) && relation_is_endorelation(x)))
        stop("Argument 'x' must be an endorelation.")
    if(!relation_is_crisp(x, na.rm = TRUE))
        stop("Argument 'x' must be a crisp relation.")

    R <- transitive_closure(x)

    if(!relation_is_antisymmetric(R))
        stop("Argument 'x' must be an acyclic relation.")

    x - x * R
}

### * transitive_closure

## Warshall's algorithm
transitive_closure <-
function(x)
{
    if(!(is.relation(x) && relation_is_endorelation(x)))
        stop("Argument 'x' must be an endorelation.")
    I <- relation_incidence(x)
    diag_hold <- diag(I)
    diag(I) <- 1
    is_transitive <- TRUE
    for (i in seq_len(ncol(I))) {
        tmp <- outer(I[,i], I[i,], .T.)
        if(any(is.na(tmp))) is_transitive <- NA
        I <- .S.(I, tmp)
    }
    meta <- list(is_endorelation = TRUE,
                 is_transitive = is_transitive)
    diag(I) <- diag_hold
    .make_relation_from_domain_and_incidence(.domain(x), I, meta)
}

### * reflexive_closure

reflexive_closure <-
function(x)
{
    if(!(is.relation(x) && relation_is_endorelation(x)))
        stop("Argument 'x' must be an endorelation.")
    I <- relation_incidence(x)
    if (isTRUE(all(diag(I) == 1))) return(x)
    diag(I) <- 1
    meta <- list(is_endorelation = TRUE,
                 is_reflexive = TRUE)
    .make_relation_from_domain_and_incidence(.domain(x), I, meta)
}

### * reflexive_reduction

reflexive_reduction <-
function(x)
{
    if(!(is.relation(x) && relation_is_endorelation(x)))
        stop("Argument 'x' must be an endorelation.")
    I <- relation_incidence(x)
    if (isTRUE(all(diag(I) == 0))) return(x)
    diag(I) <- 0
    meta <- list(is_endorelation = TRUE,
                 is_irreflexive = TRUE)
    .make_relation_from_domain_and_incidence(.domain(x), I, meta)
}

### * relation_trace

relation_trace <-
function(x, which)
{
    if(!(is.relation(x) && relation_is_endorelation(x)))
        stop("Argument 'x' must be an endorelation.")
    which <- match.arg(which, c("left", "right"))
    D <- .domain(x)
    x <- relation_incidence(x)
    n <- nrow(x)
    I <- matrix(1, nrow = n, ncol = n)
    if(which == "left") {
        for(k in seq_len(n))
            I <- pmin(I, outer(x[k, ], x[k, ], .I.))
        .make_relation_from_domain_and_incidence(D, I)
    } else {
        for(k in seq_len(n))
            I <- pmin(I, outer(x[, k], x[, k], .I.))
        .make_relation_from_domain_and_incidence(D, t(I))
    }
}

### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "### [*]+" ***
### End: ***
