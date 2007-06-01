### * relation

relation <-
function(domain = NULL, incidence = NULL, graph = NULL, charfun = NULL)
{
    if(sum(is.null(incidence), is.null(graph), is.null(charfun)) != 2L)
        stop("Need exactly one of 'incidence', 'graph', and 'charfun'.")

    if(!is.null(domain)) {
        ## Be nice first ...
        if(!is.list(domain) || is.set(domain))
            domain <- list(X = domain)
        ## ... and then check.
        if(!.is_valid_relation_domain(domain))
            stop("Invalid relation domain.")
    }

    if(!is.null(incidence)) {
        incidence <- as.array(incidence)
        if(!.is_valid_relation_incidence(incidence))
            stop("Invalid relation incidence.")
        size <- dim(incidence)
        if(!is.null(domain)) {
            ## Be nice first ...
            domain <- rep(domain, length.out = length(size))
            ## ... and then check.
            if(any(size != sapply(domain, length)))
                stop("Relation size mismatch between domain and incidence.")
        }
        return(.make_relation_from_domain_and_incidence(domain, incidence))
    }

    if(!is.null(graph)) {
        G <- .make_relation_graph_components(graph)
        ## Be nice and recycle domain (useful for endorelations).
        domain <- rep(domain, length.out = length(G))
        return(.make_relation_from_domain_and_graph_components(domain, G))
    }
    
    if(!is.null(charfun)) {
        if(is.null(domain))
            stop("Need domain along with characteristic function.")
        ## No recycling here, as we really do not know the arity of a
        ## function (nor is this a well-defined notion).
        I <- array(do.call(charfun, .cartesian_product(domain)),
                   dim = sapply(domain, length))
        return(.make_relation_from_domain_and_incidence(domain, I))
    }
}

### * is.relation

is.relation <-
function(x)
    inherits(x, "relation")

### * as.relation

as.relation <-
function(x)
    UseMethod("as.relation")

## Obviously.
as.relation.relation <- .identity

## Logical vectors are taken as unary relations (predicates).
as.relation.logical <-
function(x)
{
    D <- if(!is.null(nms <- names(x)))
        list(nms)
    else
        NULL
    I <- as.array(as.integer(x))
    meta <- list(is_endorelation = FALSE,
                 is_complete = all(is.finite(x)))
    .make_relation_from_domain_and_incidence(D, I, meta)
}

## Numeric vectors and ordered factors are taken as order relations.
as.relation.numeric <-
function(x)
{
    D <- if(!is.null(nms <- names(x)))
        list(nms, nms)
    else if(!any(duplicated(x)))
        rep.int(list(as.character(x)), 2L)
    else
        NULL
    I <- outer(x, x, "<=")
    meta <- if(any(is.na(x)))
        list(is_endorelation = TRUE,
             is_complete = NA,
             is_reflexive = NA,
             is_antisymmetric = NA,
             is_transitive = NA)
    else
        list(is_endorelation = TRUE,
             is_complete = TRUE,
             is_reflexive = TRUE,
             is_antisymmetric = !any(duplicated(x)),
             is_transitive = TRUE)
    .make_relation_from_domain_and_incidence(D, I, meta)
}
as.relation.integer <- as.relation.numeric
as.relation.ordered <- as.relation.numeric

## Unordered factors are taken as equivalence relations.
as.relation.factor <-
function(x)
{
    D <- if(!is.null(nms <- names(x)))
        list(nms, nms)
    else if(!any(duplicated(x)))
        rep.int(list(as.character(x)), 2L)
    else
        NULL
    I <- outer(x, x, "==")
    meta <- if(any(is.na(x)))
        list(is_endorelation = TRUE,
             is_complete = NA,
             is_reflexive = NA,
             is_symmetric = NA,
             is_transitive = NA)
    else
        list(is_endorelation = TRUE,
             is_complete = all(is.finite(x)),
             is_reflexive = TRUE,
             is_symmetric = TRUE,
             is_transitive = TRUE)
    .make_relation_from_domain_and_incidence(D, I, meta)
}

## Matrices and arrays are taken as incidences of relations, provided
## that they are feasible.
as.relation.matrix <-
function(x)
{
    if(!.is_valid_relation_incidence(x))
        stop("Invalid relation incidence.")
    meta <- list(is_endorelation =
                 .relation_is_endorelation_using_incidence(x))
    .make_relation_from_domain_and_incidence(dimnames(x), x, meta)
}
as.relation.array <- as.relation.matrix

## Data frames are taken as relation graph components.
as.relation.data.frame <-
function(x)
{
    ## Get the domain.
    D <- lapply(x, unique)
    names(D) <- names(x)
    ## Get the incidences.
    I <- .make_incidence_from_domain_and_graph_components(D, x)
    ## And put things together.
    .make_relation_from_domain_and_incidence(D, I)
}

## Package clue: cl_partition objects.
as.relation.cl_partition <-
function(x)
    as.relation(factor(clue::cl_class_ids(x)))

### * Relation methods

### ** as.data.frame.relation

as.data.frame.relation <-
function(x, row.names = NULL, ...)
{
    ## Get the "raw" graph components.
    out <- .make_relation_graph_components(x)
    names(out) <-
        .make_domain_names_from_relation_graph_components(out,
                                                          relation_is_endorelation(x))
    ## Flatten.
    out <- lapply(out, unlist, recursive = FALSE)
    ## And put into "some kind of" data frame.
    .make_data_frame_from_list(out, row.names)
}

### ** as.tuple.relation

as.tuple.relation <-
function(x)
{
    D <- as.tuple(relation_domain(x))
    G <- as.set(relation_graph(x))
    if(is.null(names(D)))
        names(D) <- names(G)
    names(G) <- NULL
    pair(Domain = D, Graph = G)
}

as.tuple.relation_domain <-
function(x)
    structure(x, class = "tuple")

### ** print.relation

print.relation <-
function(x, ...)
{
    a <- .arity(x)
    s <- paste(.size(x), collapse = " x ")
    if(a == 1L)
        writeLines(gettextf("A unary relation of size %s.", s))
    else if(a == 2L)
        writeLines(gettextf("A binary relation of size %s.", s))
    else
        writeLines(gettextf("A %d-ary relation of size %s.", a, s))
    invisible(x)
}

### * Group methods and related.

## Here is what we do.

## * Comparisons are obvious.
## * We use & and | for intersection and union.
## * We use min/max for meet and join (which of course is the same as
##   the above).
## * We use * for the composition and unary ! for the converse.
## * Finally, t() is used for the dual.

Summary.relation <-
function(..., na.rm = FALSE)
{
    ok <- switch(.Generic, max = , min = , range = TRUE, FALSE)
    if(!ok)
        stop(gettextf("Generic '%s' not defined for \"%s\" objects.",
                      .Generic, .Class))
    args <- list(...)
    x <- relation_ensemble(list = args)
    switch(.Generic,
           "min" = .relation_meet(x),
           "max" = .relation_join(x),
           "range" = {
               relation_ensemble(min = .relation_meet(x),
                                 max = .relation_join(x))
           })
}

Ops.relation <-
function(e1, e2)
{
    if(nargs() == 1L) {
        if(!(as.character(.Generic) %in% "!"))
            stop(gettextf("Unary '%s' not defined for \"%s\" objects.",
                          .Generic, .Class))
        return(.make_relation_from_domain_and_incidence(.domain(e1),
                                                        1 - .incidence(e1)))
    }

    ## In addition to comparisons, we support & | * + - %/% %%.
    if(!(as.character(.Generic)
         %in% c("<", "<=", ">", ">=", "==", "!=",
                "&", "|", "*", "+", "-", "%/%", "%%")))
        stop(gettextf("Generic '%s' not defined for \"%s\" objects.",
                      .Generic, .Class))

    switch(.Generic,
           "+" = return(relation_union(e1, e2)),
           "-" = return(relation_complement(e1, e2)),
           "%/%" = return(relation_division(e1, e2)),
           "%%" = return(relation_remainder(e1, e2))
           )
           
    D1 <- relation_domain(e1)
    D2 <- relation_domain(e2)
    I1 <- .incidence(e1)
    I2 <- .incidence(e2)

    ## Composition (*) is only defined for binary relations with
    ## matching 2nd/1st domain elements.
    if(as.character(.Generic) == "*") {
        if((length(D1) != 2L)
           || (length(D2) != 2L)
           || !set_is_equal(D1[[2L]], D2[[1L]]))
            stop("Composition of given relations not defined.")
        ## When composing indicidences, need the same *internal* order
        ## for D2[[1L]] as for D1[[2L]].
        I <- ((I1 %*% I2[match(D1[[2L]], D2[[1L]]), ]) > 0)
        ## The composition has domain (D1[[1L]], D2[[2L]]) (and
        ## appropriate names).  Information about auto-generation of
        ## domains is currently ignored.
        D1 <- .domain(e1)
        D2 <- .domain(e2)
        D <- list(D1[[1L]], D2[[2L]])
        if(!is.null(nms <- names(D1)))
            names(D)[1L] <- nms[1L]
        if(!is.null(nms <- names(D2)))
            names(D)[2L] <- nms[2L]
        return(.make_relation_from_domain_and_incidence(D, I))
    }

    ## In the remaining cases, the relations must have equal domains.
    if(!.domain_is_equal(D1, D2))
        stop("Relations need equal domains.")
    ## Use the same *internal* order for D2 as for D1.
    I2 <- .reorder_incidence(I2, .match_domain_components(D1, D2))
    ## And now do it.
    D <- .domain(e1)
    switch(.Generic,
           "<=" = all(I1 <= I2),
           "<"  = all(I1 <= I2) && any(I1 < I2),
           ">=" = all(I1 >= I2),
           ">"  = all(I1 >= I2) && any(I1 > I2),
           "==" = all(I1 == I2),
           "!=" = any(I1 != I2),
           "&" = .make_relation_from_domain_and_incidence(D, I1 & I2),
           "|" = .make_relation_from_domain_and_incidence(D, I1 | I2)
           )
}

t.relation <-
function(x)
{
    if(!relation_is_binary(x))
        stop("Can only compute duals of binary relations.")
    .make_relation_from_domain_and_incidence(rev(.domain(x)),
                                             t(.incidence(x)))
}

### * Relation representations

### ** .make_relation_by_domain_and_incidence

## (That's all for the time being ...)

.make_relation_by_domain_and_incidence <-
function(D, I)
{
    ## Canonicalize a "valid" incidence, but generate a valid domain if
    ## needed.  Note the difference: valid domains are ok as is, but
    ## currently e.g. "valid" incidences include vectors and logicals.
    I <- as.array(I)
    if(is.logical(I)) I <- I + 0L
    if(!.is_valid_relation_domain(D)) {
        D <- dimnames(I)
        if(!.is_valid_relation_domain(D)) {
            D <- lapply(dim(I), function(s) as.character(seq_len(s)))
            ## Just to make a point ...
            attr(D, "auto") <- TRUE
        }
    }
    ## No dimnames in incidence.
    dimnames(I) <- NULL
    size <- dim(I)

    structure(list(domain = D,
                   incidence = I,
                   .arity = length(size),
                   .size = size),
              class = "relation_by_domain_and_incidence")
}

### * Relation generators

### ** .make_relation_from_representation_and_meta

.make_relation_from_representation_and_meta <-
function(x, meta = NULL)
    .make_container(x, "relation", meta)

### ** .make_relation_from_domain_and_incidence

.make_relation_from_domain_and_incidence <-
function(D, I, meta = NULL)
{
    R <- .make_relation_by_domain_and_incidence(D, I)
    .make_relation_from_representation_and_meta(R, meta)
}

### ** .make_relation_from_domain_and_graph_components

.make_relation_from_domain_and_graph_components <-
function(D, G)
{
    ## (Assuming that G really has the graph *components* as obtained by
    ## .make_relation_graph_components().)
    values <- lapply(G, as.set)

    ## Get the domain.
    if(!is.null(D)) {
        if(length(D) != length(G))
            stop("Relation arity mismatch between domain and graph.")
        D <- lapply(D, as.set)
        ## Check containment.
        if(!all(mapply(function(s, t) all(s %in% t),
                       .transform_factors_into_characters(values),
                       .transform_factors_into_characters(D))))
            stop("Invalid graph with out-of-domain elements.")
    }
    else
        D <- values

    ## Get the incidences.
    I  <- .make_incidence_from_domain_and_graph_components(D, G)
    
    ## And put things together.
    .make_relation_from_domain_and_incidence(D, I)
}


### * Utilities

### ** .canonicalize_relation

.canonicalize_relation <-
function(R, D, pos = NULL)
{
    ## For a relation R with domain known to equal D in the sense that
    ## the respective domain elements are the same sets (as tested for
    ## by .domain_is_equal()), "canonicalize" R to have its domain
    ## elements use the same *internal* order as the elements of D.
    if(is.null(pos)) {
        pos <- .match_domain_components(lapply(D, as.set),
                                        relation_domain(R))
        ## If already canonical, do nothing.
        if(!any(sapply(pos, is.unsorted))) return(R)
    }
    ## Use the reference domain, and reorder incidences.
    I <- .reorder_incidence(.incidence(R), pos)
    meta <- relation_properties(R)
    .make_relation_from_domain_and_incidence(D, I, meta)
}

### ** .make_data_frame_from_list

.make_data_frame_from_list <-
function(x, row.names = NULL)
{
    if(length(x) > 0L) {
        len <- length(x[[1L]])
        row.names <- if(!is.null(row.names)) {
            ## Do some checking ...
            if(length(row.names) != len)
                stop("Incorrect length of given 'row.names'.")
            as.character(row.names)
        }
        else
            .set_row_names(len)
    }
    structure(x, class = "data.frame", row.names = row.names)
}

### * .make_domain_names_from_relation_graph_components

.make_domain_names_from_relation_graph_components <-
function(x, endorelation = FALSE) {
    ## In case the domains do not have names, use X_i.
    ## Needed for as.data.frame.relation();
    ## Not sure if we want this for relation_graph(), too.
    nms <- names(x)
    arity <- length(x)
    if (is.null(nms)) {
      nms <- if (endorelation) {
        rep("X", 2L)
        ## (Assuming that endorelations are always binary.)
      }
      else if (arity == 1L) "X"
      else sprintf("X_%d", seq_len(arity))
    }
    nms
}

### ** .make_incidence_from_domain_and_graph_components

.make_incidence_from_domain_and_graph_components <-
function(D, G, size = NULL)
{
    if(is.null(size)) size <- sapply(D, length)
    I <- array(0, size)
    I[rbind(mapply(match, G, D))] <- 1
    I
}


### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "### [*]+" ***
### End: ***

