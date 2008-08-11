### * plot.relation

plot.relation <-
function(x,
         attrs = list(graph = list(rankdir = "BT"),
                      edge = list(arrowsize = "0"),
                      node = list(shape = "rectangle", fixedsize = FALSE)),
         limit = 6L, labels = NULL,
         ...)
{
    if(!relation_is_endorelation(x))
        stop("Plot only available for endorelations.")
    if(!relation_is_crisp(x))
        stop("Plot only available for crisp relations.")
    if(!(relation_is_transitive(x)
         && (relation_is_antisymmetric(x) || relation_is_complete(x))))
        stop("Plot only available for antisymmetric or complete transitive relations.")
    if(!require("Rgraphviz"))
        stop("Need Rgraphviz package (obtainable from bioconductor.org))!")

    ## if x is a preference, use dual of complement instead
    if (!relation_is_antisymmetric(x))
      x <- t(!x)

    ## compute transitive reduction to avoid cluttered graph,
    ## and extract incidence
    I <- unclass(relation_incidence(transitive_reduction(x), limit = limit))

    ## transform to graphViz-compatible incidence
    if (is.null(labels))
        labels <- labels(I)
    dimnames(I) <- lapply(labels, .make_unique_labels)

    plot(as(I, "graphNEL"), attrs = attrs, ...)
}

.make_unique_labels <-
function(x)
{
    I <- which(duplicated(x))
    x[I] <- paste(x[I], ".", seq_along(I), sep="")
    x
}


### * plot.relation_ensemble

plot.relation_ensemble <-
function(x, attrs = list(list(graph = list(rankdir = "BT"),
                              edge = list(arrowsize = "0"),
                              node = list(shape = "rectangle",
                                          fixedsize = FALSE))),
         ..., layout = NULL)
{
    ## Why not?
    ## Of course, we maybe should only have one thing ...
    if(!is.relation_ensemble(x))
        stop("Wrong class.")
    if(!all(sapply(x,
                   function(e)
                   (relation_is_crisp(e)
                    && relation_is_endorelation(e)
                    && relation_is_transitive(e)
                    && (relation_is_antisymmetric(e)
                        || relation_is_complete(e))))))
        stop("Plotting only available for ensembles of antisymmetric or complete transitive crisp relations.")
    ## Make things a bit more efficient.
    x <- unclass(x)
    if(!require("Rgraphviz"))
        stop("Need Rgraphviz package (obtainable from bioconductor.org))!")
    ## Number of elements.
    n <- length(x)
    ## Layout.
    byrow <- TRUE
    if(is.null(layout)) {
        nc <- ceiling(sqrt(n))
        nr <- ceiling(n / nc)
    }
    else {
        layout <- c(as.list(layout), byrow)[seq_len(3)]
        if(is.null(names(layout)))
            names(layout) <- c("nr", "nc", "byrow")
        nr <- layout[["nr"]]
        nc <- layout[["nc"]]
        byrow <- layout[["byrow"]]
    }
    op <- if(byrow)
        par(mfrow = c(nr, nc))
    else
        par(mfcol = c(nr, nc))
    on.exit(par(op))
    attrs <- rep(attrs, length.out = length(x))
    for(i in seq_along(x)) {
        ## if x is a preference, use dual of complement instead
        if (!relation_is_antisymmetric(x[[i]]))
          x[[i]] <- t(!x[[i]])

        ## compute transitive reduction to avoid cluttered graph
        ## and extract incidence
        I <- unclass(relation_incidence(transitive_reduction(x[[i]])))
        dimnames(I) <- lapply(labels(I), .make_unique_labels)
        plot(as(I, "graphNEL"), attrs = attrs[[i]], ...)
    }
}

### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "### [*]+" ***
### End: ***
