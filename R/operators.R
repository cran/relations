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

transitive_reduction <-
function(x)
{
    if(!(is.relation(x) && relation_is_endorelation(x)))
        stop("Argument 'x' must be an endorelation.")
    if(!relation_is_crisp(x, na.rm = TRUE))
        stop("Argument 'x' must be a crisp relation.")

    diag_hold <- diag(relation_incidence(x))
    diag(relation_incidence(x)) <- 0

    R <- transitive_closure(x)

    if(!relation_is_antisymmetric(R)) {
        ## handle cyclic case:

        ## compute connected components and leaders
        I <- relation_incidence(x)
        scc <- .connected_components(I)
        leaders <- sapply(scc, min)

        ## compute transitive reduction from condensation
        M <- .condensation_incidences(I, scc, leaders)
        M <- .T.(M, .N.(M %*% .transitive_closure_incidences(M)))

        ## merge with component representation
        I <- .component_representation_incidences(I, scc)
        I[leaders, leaders] <- .S.(I[leaders, leaders], M)

        diag(I) <- diag_hold
        .make_relation_from_domain_and_incidence(.domain(x), I)
    } else {
        x <- x - x * R
        diag(relation_incidence(x)) <- diag_hold
        x
    }
}

### * transitive_closure

transitive_closure <-
function(x)
{
    if(!(is.relation(x) && relation_is_endorelation(x)))
        stop("Argument 'x' must be an endorelation.")
    I <- .transitive_closure_incidences(relation_incidence(x))
    .make_relation_from_domain_and_incidence(.domain(x), I, attr(I, "meta"))
}

## Warshall's algorithm
.transitive_closure_incidences <-
function(I)
{
    diag_hold <- diag(I)
    diag(I) <- 1
    is_transitive <- TRUE
    for (i in seq_len(ncol(I))) {
        tmp <- outer(I[,i], I[i,], .T.)
        if(any(is.na(tmp))) is_transitive <- NA
        I <- .S.(I, tmp)
    }
    diag(I) <- diag_hold

    structure(I,
              meta = list(is_endorelation = TRUE,
                          is_transitive = is_transitive)
              )
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

### * relation_connected_components

relation_connected_components <-
function(x, type = c("strongly", "weakly"))
{
    if(!relation_is_endorelation(x) && !isTRUE(relation_is_crisp(x)))
        stop("Argument 'x' must be a crisp endorelation without missings.")

    type <- match.arg(type)

    I <- relation_incidence(x)
    if (type == "weakly")
        I <- .S.(I, t(I)) ## symmetric completion

    scc <- .connected_components(I)
    vertices <- as.list(.domain(x)[[1L]])
    leaders <- vertices[sapply(scc, min)]
    names(scc) <- LABELS(leaders)

    structure(lapply(scc, function(i) as.set(vertices[i])),
              leaders = leaders,
              class = "relation_classes_of_objects"
              )
}

.connected_components <-
function(I)
{
    tarjan <- function(v) {
        lowlink[v] <<- indices[v] <<- index
        index <<- index + 1L
        stack <<- c(v, stack) # "push" vertice on stack

        for (w in which(I[v,] > 0)) # for all successors of v ...
            if (indices[w] < 0L) {
                ## successor has not been visited -> recurse on it
                tarjan(w)
                lowlink[v] <<- min(lowlink[v], lowlink[w])
            }
            else if (any(w == stack))
                ## successor is on stack, hence in current strongly connected component
                lowlink[v] <<- min(lowlink[v], indices[w])

        ## if v is a root node ("Leader"), pop the stack and generate strongly connected component
        if (lowlink[v] == indices[v]) {
            to <- which(stack == v)[1L]
            scc <<- c(scc, list(stack[seq_len(to)]))
            stack <<- stack[-seq_len(to)]
        }
    }

    N <- ncol(I)

    indices <- lowlink <- rep.int(-1L, N)
    index <- 0L
    stack <- c()
    scc <- list()

    for(i in seq_len(N))
        if (indices[i] < 0L)
            tarjan(i)

    lapply(scc, sort)
}

### . relation_condensation

relation_condensation <-
function(x)
{
    if(!relation_is_endorelation(x) && !isTRUE(relation_is_crisp(x)))
        stop("Argument 'x' must be a crisp endorelation without missings.")

    I <- relation_incidence(x)
    scc <- .connected_components(I)
    leaders <- sapply(scc, min)

    M <- .condensation_incidences(I, scc, leaders)
    D <- rep.int(list(as.list(.domain(x)[[1L]])[leaders]), 2L)

    .make_relation_from_domain_and_incidence(D, M)
}

.condensation_incidences <-
function(I, scc, leaders)
{
    N <- length(scc)
    M <- matrix(0, nrow = N, ncol = N)

    s <- seq_len(N)
    for (i in s)
        for (j in s[-i])
            M[i, j] <- any(I[scc[[i]], scc[[j]]] > 0)

    M
}

### * relation_component_representation

relation_component_representation <-
function(x)
{
    if(!relation_is_endorelation(x) && !isTRUE(relation_is_crisp(x)))
        stop("Argument 'x' must be a crisp endorelation without missings.")

    I <- relation_incidence(x)
    scc <- .connected_components(I)

    `relation_incidence<-`(x, .component_representation_incidences(I, scc))
}

.component_representation_incidences<-
function(I, scc)
{
   M <- `[<-`(I, 0)
   for (i in scc)
       if (length(i) > 1L)
           for (j in seq_along(i))
               M[i[j], c(i, i[1L])[j + 1L]] <- 1
   M
}

### symmetric and asymmetric parts as methods for sym() and asy()
### generics.

### * sym

sym <-
function(x)
    UseMethod("sym")

sym.relation <-
function(x)
{
    if(!relation_is_endorelation(x))
        stop("Argument 'x' must be an endorelation.")

    ## For crisp relations, the symmetric part is
    ##   R \cap t(R)
    ## which is symmetric if there are no missings.
    ## The natural extension to fuzzy relations is
    ##   T( R(x, y), R(y, x) )
    ## which is clearly symmetric in x and y.  Apparently, this
    ## definition is at least given in
    ##   Ovchinnikov (2000)
    ##     <http://link.springer.com/chapter/10.1007%2F978-1-4615-4429-6_5>
    ## so we also use it here.

    I <- .incidence(x)

    meta <- if(any(is.na(I))) NULL else list(is_symmetric = TRUE)

    .make_relation_from_domain_and_incidence(.domain(x),
                                             .T.(I, t(I)),
                                             meta)
}

### * asy

asy <-
function(x)
    UseMethod("asy")

asy.relation <-
function(x)
{
    ## For crisp relations, the asymmetric part is
    ##   R \cap not(t(R))
    ## which is asymmetric (and hence irreflexive) if there are no
    ## missings.
    ## A possible extension to fuzzy relations is
    ##   T( R(x, y), N(R(y, x)) )
    ## but apparently this does not necessarily yield asymmetry in the
    ## sense of satisfying T( S(x, y), S(y, x) ) = 0, i.e.
    ##   T ( T(R(x, y), N(R(y, x))), T(R(y, x), N(R(x, y))) ) = 0.
    ## (It would be interesting to investigate this for common choices
    ## for T and N, and check the literature once more.)
    ## So for now we only do asymmetric parts of *crisp* endorelations
    ## which requires no missings (as otherwise we cannot know whether
    ## we have a crisp relation).

    if(!relation_is_endorelation(x) && !isTRUE(relation_is_crisp(x)))
        stop("Argument 'x' must be a crisp endorelation without missings.")

    I <- .incidence(x)

    .make_relation_from_domain_and_incidence(.domain(x),
                                             pmin(I, 1 - t(I)),
                                             list(is_asymmetric = TRUE,
                                                  is_irreflexive = TRUE))
}

### * codual

## This used to be 'dual', but then this is not use consistently in the
## literature, with the (fuzzy) preference modeling community typically
## using 'dual' for the complement of the transpose, e.g.,
##   Ovchinnikov (2000)
##     <http://link.springer.com/chapter/10.1007%2F978-1-4615-4429-6_5>
## or chapter 2 in Fodor & Roubens (1994)
##     <http://www.springer.com/us/book/9780792331162>
## but Fishburn and the utility theory community using 'dual' as synonym
## for transpose/inverse.
## To avoid confusion, we now use 'codual' for the complement of the
## transpose, with
##   Clark (1990)
##     <doi:10.1016/0165-4896(90)90065-F>
##     <http://www.sciencedirect.com/science/article/pii/016548969090065F>
## as references.

codual <-
function(x, ...)
    UseMethod("codual")

codual.relation <-
function(x, ...)
{
    if(!relation_is_binary(x))
        stop("Argument 'x' must be a binary relation.")
    I <- .incidence(x)
    meta <- if(relation_is_endorelation(x)) {
        ## Predicates for the codual relation of an endorelation R can
        ## be inferred from those of R, see e.g. Fodor & Roubens, "Fuzzy
        ## Preference Modelling and Multicriteria Decision Support",
        ## Table 2.2, page 41.
        ## <http://www.springer.com/us/book/9780792331162>
        ## For valued relations, only the correspondencies
        ##   reflexive <-> irreflexive, symmetric <-> symmetric
        ## are always true: the others require a deMorgan triple of
        ## fuzzy connectives N/T/S.
        db <- c(is_reflexive = "is_irreflexive",
                is_irreflexive = "is_reflexive",
                is_symmetric = "is_symmetric")
        if(fuzzy_logic_predicates()$is_de_Morgan_triple) {
            db <-
                c(db,
                  is_antisymmetric = "is_complete",
                  is_complete = "is_antisymmetric",
                  is_asymmetric = "is_strongly_complete",
                  is_strongly_complete = "is_asymmetric",
                  is_transitive = "is_negatively_transitive",
                  is_negatively_transitive = "is_transitive",
                  is_Ferrers = "is_Ferrers",
                  is_semitransitive = "is_semitransitive"
                  )
        }
        predicates <-
            names(Filter(function(e) identical(e, TRUE),
                         relation_properties(x)[names(db)]))
        c(list(is_endorelation = TRUE),
          .structure(as.list(rep.int(TRUE, length(predicates))),
                     names = db[predicates]))
    } else NULL
    .make_relation_from_domain_and_incidence(.domain(x), .N.(t(I)), meta)
}

### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "### [*]+" ***
### End: ***
