### * plot.relation

plot.relation <-
function(x, attrs = list(edge = list(arrowsize = "0")), ...)
{
    if(!(is.relation(x) && relation_is_endorelation(x)))
        stop("Plot only available for endorelations.")
    if(!(relation_is_transitive(x) && relation_is_antisymmetric(x)))
        stop("Plot only available for antisymmetric and transitive relations.")
    if(!require("Rgraphviz"))
        stop("Need Rgraphviz package (obtainable from bioconductor.org))!")
      
    x <- unclass(relation_incidence(transitive_reduction(x)))
    dimnames(x) <- labels(x)
    plot(as(t(x), "graphNEL"), attrs = attrs, ...)
}

### * plot.relation_ensemble

plot.relation_ensemble <-
function(x, attrs = list(edge = list(arrowsize = "0")), ...,
         layout = NULL)
{
    ## Why not?
    ## Of course, we maybe should only have one thing ...
    if(!is.relation_ensemble(x))
        stop("Wrong class.")
    if(!all(sapply(x,
                   function(e)
                   (relation_is_endorelation(e)
                    && relation_is_transitive(e)
                    && relation_is_antisymmetric(e)))))
        stop("Plotting only available for ensembles of antisymmetric and transitive relations.")
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
    for(i in seq_along(x)) {
        I <- unclass(relation_incidence(transitive_reduction(x[[i]])))
        dimnames(I) <- labels(I)
        plot(as(t(I), "graphNEL"), attrs = attrs, ...)
    }
}    

### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "### [*]+" ***
### End: ***
