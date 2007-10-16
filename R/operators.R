### * transitive_reduction

transitive_reduction <-
function(x)
{
    if(!(is.relation(x) && relation_is_endorelation(x)))
        stop("Argument 'x' must be an endorelation.")
    if(!relation_is_crisp(x))
        stop("Argument 'x' must be a crisp relation.")
    I <- relation_incidence(x)
    diag(I) <- 0
    has_children <- rowSums(I) > 0
    for (i in which(has_children)) {
        visited <- rep.int(FALSE, ncol(I))
        for (j in which(has_children & (I[i,] > 0)))
            for (k in which(I[j,] > 0))
                if(I[i,k] && !visited[k]) {
                    I[i,k] <- 0
                    visited[k] <- TRUE
                }
    }
    meta <- list(is_endorelation = TRUE,
                 is_transitive = TRUE)
    .make_relation_from_domain_and_incidence(.domain(x), I, meta)
}

### * transitive_closure

## Warshall's algorithm
transitive_closure <-
function(x)
{
    if(!(is.relation(x) && relation_is_endorelation(x)))
        stop("Argument 'x' must be an endorelation.")
    I <- relation_incidence(x)
    diag(I) <- 1
    for (i in seq_len(ncol(I)))
        I <- .S.(I, outer(I[,i], I[i,], ".T."))
    meta <- list(is_endorelation = TRUE,
                 is_transitive = TRUE)
    .make_relation_from_domain_and_incidence(.domain(x), I, meta)
}

### * reflexive_closure

reflexive_closure <-
function(x)
{
    if(!(is.relation(x) && relation_is_endorelation(x)))
        stop("Argument 'x' must be an endorelation.")
    I <- relation_incidence(x)
    if(identical(all(diag(I) == 1), TRUE))
        return(x)
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
    if(identical(all(diag(I) == 0), TRUE))
        return(x)
    diag(I) <- 0
    meta <- list(is_endorelation = TRUE,
                 is_irreflexive = TRUE)
    .make_relation_from_domain_and_incidence(.domain(x), I, meta)
}

### * relation_left_trace

relation_left_trace <-
function(x)
{
    if(!(is.relation(x) && relation_is_endorelation(x)))
        stop("Argument 'x' must be an endorelation.")
    D <- .domain(x)
    x <- relation_incidence(x)
    n <- nrow(x)
    I <- matrix(1, nrow = n, ncol = n)
    for(k in seq_len(n))
        I <- pmin(I, outer(x[k, ], x[k, ], .I.))
    .make_relation_from_domain_and_incidence(D, I)
}

### * relation_right_trace

relation_right_trace <-
function(x)
{
    if(!(is.relation(x) && relation_is_endorelation(x)))
        stop("Argument 'x' must be an endorelation.")
    D <- .domain(x)
    x <- relation_incidence(x)
    n <- nrow(x)
    I <- matrix(1, nrow = n, ncol = n)
    for(k in seq_len(n))
        I <- pmin(I, outer(x[, k], x[, k], .I.))
    .make_relation_from_domain_and_incidence(D, t(I))
}


### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "### [*]+" ***
### End: ***
