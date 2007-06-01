### Relation graph.

## <NOTE>
## The sparser the relation incidences, the less efficient/reasonable
## the implementation below ...
## We could also try improving printing ...
## </NOTE>

### * relation_graph

relation_graph <-
function(x) {
    if(!is.relation(x))
        stop("Argument 'x' must be a relation.")
    out <- .make_relation_graph_components(x)
    nms <- names(out) 
    ## Note that we no longer use auto-generated names for the relation
    ## graph.
    names(out) <- NULL
    .make_relation_graph(.make_set_of_tuples_from_relation_graph_components(out),
                         domain_names = nms)
}

.make_relation_graph <-
function(x, domain_names = NULL)
    structure(x,
              domain_names = domain_names,
              class = c("relation_graph", class(x)))

### * relation_graph<-

"relation_graph<-" <-
function(x, value)
{
    if(!is.relation(x))
        stop("Argument 'x' must be a relation.")
    .make_relation_from_domain_and_graph_components(.domain(x),
                                                    .make_relation_graph_components(value))
}

### * relation_graph methods

print.relation_graph <-
function(x, ...)
{
    writeLines("Relation graph:\n")
    if (length(x) < 1L)
       writeLines("The empty set.")
    else {
       domain_names<- attr(x, "domain_names")
       writeLines(sprintf("A set with %s%s:",
                          .ntuple(x[[1L]], plural = length(x) > 1),
                          if (is.null(domain_names)) "" else 
                          paste(" ", format(as.tuple(domain_names)), sep = "")
                          )
                  )
       for (i in x) print(i)
    }
    invisible(x)
}

as.list.relation_graph <-
function(x, ...)
    unclass(x)

### * .make_relation_graph_components

.make_relation_graph_components <-
function(x)
    UseMethod(".make_relation_graph_components")

.make_relation_graph_components.default <-
function(x)
    stop("Not implemented.")

.make_relation_graph_components.list <-
.make_relation_graph_components.set <-
function(x)
{
    l <- length(x[[1L]])
    if (!all(sapply(x, length) == l))
        stop("All elements need to be of same length!")
    if (l < 2L)
        as.list(x)
    else {
        ## <NOTE>
        ## We could try making this more efficient by building a big
        ## list and slicing into components ...
        lapply(seq_len(l), function(i) sapply(x, function(j) j[[i]]))
        ## </NOTE>
    }
}

.make_relation_graph_components.relation <-
function(x)
{
    I <- relation_incidence(x)
    D <- relation_domain(x)
    ind <- as.logical(c(I))
    ## Set missings to zero.
    ind[is.na(ind)] <- FALSE
    ## <NOTE>
    ## We would like to use
    ##   out <- do.call("expand.grid", D)[ind, ]
    ## but this turns character vectors into factors with no way of
    ## turning this off.
    ## (Yes this gains storage efficiency but ...)
    ## Seems we need our own simplified expand.grid() ...
    ## Of course, we could try to generate only the replications needed
    ## (as we already have ind).
    structure(lapply(.cartesian_product(D), "[", ind),
              names = names(D))
    ## </NOTE>
}

.make_relation_graph_components.data.frame <-
function(x)
    as.list(x)

.make_relation_graph_components.matrix <-
function(x)
    .split_into_components(x)

### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "### [*]+" ***
### End: ***
