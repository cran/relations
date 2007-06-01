### Relational algebra-like operations.

### * relation_projection

relation_projection <-
function(x, margin = NULL)
{
    if (is.null(margin))
        return(x)
    D <- relation_domain(x)
    
    ind <- if (is.character(margin))
        match(margin, names(D))
    else if (is.numeric(margin))
        match(margin, seq_along(D))
    else NULL
    
    if (!length(ind) || any(is.na(ind)))
        stop("Invalid projection margin.")
    
    I <- apply(relation_incidence(x), ind, any)
    .make_relation_from_domain_and_incidence(.domain(x)[ind], I)
}

### * relation_selection (= subset)

## <FIXME>
## Do we really have to copy the code from subset.data.frame()?
## </FIXME>
relation_selection <-
function(x, subset)
{
    .stop_if_not_relation_has_unique_domain_names(x)
    df <- as.data.frame(x)
    if (missing(subset)) 
        r <- TRUE
    else {
        e <- substitute(subset)
        r <- eval(e, df, parent.frame())
        if (!is.logical(r)) 
            stop("'subset' must evaluate to logical")
        r <- r & !is.na(r)
    }
    relation_graph(x) <- df[r,]
    x
}

"%><%" <- 
relation_cartesian <-
function(x, y, ...)
{
    if (missing(y)) return(x)
    l <- list(...)
    if (length(l) > 0L)
        return(Recall(x, Recall(y, ...)))
    .make_relation_from_domain_and_incidence(c(relation_domain(x),
                                               relation_domain(y)),
                                             outer(relation_incidence(x),
                                                   relation_incidence(y))
                                             )
}

### * relation_union

## For union of relations, we also allow non-identical domain levels.

"%U%" <-
relation_union <-
function(x, y, ...)
{
    ## handle 1 and more than 2 arguments
    if (missing(y)) return(x)
    l <- list(...)
    if (length(l) > 0L)
        return(Recall(x, Recall(y, ...)))

    ## check arities
    ## (use relation_domain() instead of .domain() here,
    ## since we need tuples of sets for the combination.)
    D1 <- relation_domain(x)
    D2 <- relation_domain(y)
    if (length(x) != length(y))
        stop("Relation arity mismatch.")

    ## combine domains
    D12 <- mapply(c, D1, D2, SIMPLIFY = FALSE)

    ## combine graph components
    GC <- mapply(c,
                 .make_relation_graph_components(x),
                 .make_relation_graph_components(y),
                 SIMPLIFY = FALSE)
    
    ## And put things together
    .make_relation_from_domain_and_graph_components(D12, GC)
}

### * relation_complement

## For the complement, we also allow non-identical domain levels.

relation_complement <-
function(x, y)
{
    ## check arities
    D <- .domain(x)
    if (length(D) != length(.domain(y)))
        stop("Relation arity mismatch.")

    I <- relation_incidence(x)
    ## <FIXME>
    ## match() behaves strangely on list of factors with different levels
    ## so we have to compare character strings in this case
    G <- .transform_factors_into_characters(.make_relation_graph_components(y))
    D <- .transform_factors_into_characters(D)
    ## </FIXME>
    I[rbind(mapply(match, G, D))] <- 0
  
    ## And put things together
    .make_relation_from_domain_and_incidence(.domain(x), I)
}

### * relation_division

relation_division <-
function(x, y)
{
    .stop_if_not_relation_has_unique_domain_names(x)
    .stop_if_not_relation_has_unique_domain_names(y)

    if (length(relation_graph(y)) < 1L)
        stop("Division by empty relations not defined.")

    dx <- relation_domain_names(x)
    dy <- relation_domain_names(y)

    if (!all(dy %in% dx))
        stop("Divisor domain must be a subset of the dividend domain.")
  
    ## find attributes unique to x
    dxunique <- dx[!dx %in% dy]
    if (length(dxunique) < 1L)
        stop("Dividend needs at least one unique domain.")
    
    ## create projection of x to its unique attributes
    xunique <- relation_projection(x, dxunique)
    
    ## compute "maximum" set of tuples
    T <- relation_cartesian(xunique, y)
    
    ## remove actual set of tuples, and remove the projection
    ## of the remaining sets from the dividend
    relation_complement(xunique,
                        relation_projection(relation_complement(T, x),
                                            dxunique))
}

### relation_remainder

relation_remainder <-
function(x, y)
    relation_complement(x,
                        relation_cartesian(relation_division(x, y), y))

### * relation_join et al

"%|><|%" <-
relation_join <-
function(x, y, ...)
{
    .stop_if_not_relation_has_unique_domain_names(x)
    .stop_if_not_relation_has_unique_domain_names(y)
    tmp <- merge(X <- as.data.frame(x), Y <- as.data.frame(y), ...)
    nms <- unique(c(names(X), names(Y)))
    if(nrow(tmp) < 1L)
        stop("Join is empty!")
    as.relation(tmp[,nms])
}

"%><=%" <-
function(x, y, ...)
    relation_join(x, y, all.y = TRUE, ...)

"%=><%" <-
function(x, y, ...)
    relation_join(x, y, all.x = TRUE, ...)

"%=><=%" <-
function(x, y, ...)
    relation_join(x, y, all = TRUE, ...)

"%|><%" <-
relation_semijoin <-
function(x, y, ...)
    relation_projection(relation_join(x, y, ...),
                        relation_domain_names(x))

"%><|%" <-
function(x, y, ...)
    relation_semijoin(y, x, ...)

"%|>%" <-
relation_antijoin <-
function(x, y, ...)
    x - relation_semijoin(x, y, ...)

"%<|%" <-
function(x, y, ...)
    relation_antijoin(y, x, ...)

### * relation_symdiff

relation_symdiff <-
function(x, y)
    relation_union(relation_complement(x, y),
                   relation_complement(y, x))

### * relation_intersection

relation_intersection <-
function(x, y)
    relation_complement(relation_union(x, y), relation_symdiff(x, y))

### * .stop_if_not_relation_has_unique_domain_names

.stop_if_not_relation_has_unique_domain_names <-
function(x)
{
    nms <- relation_domain_names(x)
    if(is.null(nms) || (length(nms) < .arity(x)) || any(duplicated(nms)))
        stop("Relation(s) with unique domain names required.")
}

### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "### [*]+" ***
### End: ***
