### * plot.relation

plot.relation <-
function(x,
         attrs = list(graph = list(rankdir = "BT"),
                      edge = list(arrowsize = "0"),
                      node = list(shape = "rectangle", fixedsize = FALSE)),
         limit = 6L, labels = NULL, main = NULL,
         ...)
{
    if(!relation_is_endorelation(x))
        stop("Plotting only available for endorelations.")
    if(!relation_is_crisp(x))
        stop("Plotting only available for crisp relations.")

    if(!require("Rgraphviz"))
        stop("Plotting requires package 'Rgraphviz'.")

    ## possibly, transform relation to obtain a poset,
    ## and change main title if not specified.
    if(relation_is_linear_order(x)) {
        if (is.null(main))
            main <- "Linear Order"
    } else if(relation_is_partial_order(x)) {
        if (is.null(main))
            main <- "Partial Order"
    } else if(relation_is_weak_order(x)) {
        ## If x is a preference, use dual instead.
        x <- dual(x)
        if (is.null(main))
            main <- "Weak Order"
    } else { ## extract asymmetric part
        x <- x & dual(x)
        if (is.null(main))
            main <- "Strict Preference Part"
    }

    ## Compute transitive reduction to avoid cluttered graph, and
    ## extract incidence.
    I <- unclass(relation_incidence(transitive_reduction(x), limit = limit))

    ## Perform reflexive reduction.
    diag(I) <- 0

    ## Transform to graphViz-compatible incidence.
    if (is.null(labels))
        labels <- labels(I)
    dimnames(I) <- lapply(labels, .make_unique_labels)

    plot(as(I, "graphNEL"), attrs = attrs, main = main, ...)
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
         ..., layout = NULL, main = NULL)
{
    ## Why not?
    ## Of course, we maybe should only have one thing ...
    if(!is.relation_ensemble(x))
        stop("Wrong class.")
    if(!all(sapply(x,
                   function(e)
                   (relation_is_crisp(e) && relation_is_endorelation(e)))))
        stop("Plotting only available for ensembles of crisp endorelations.")

    ## Make things a bit more efficient.
    x <- unclass(x)
    if(!require("Rgraphviz"))
        stop("Plotting requires package 'Rgraphviz'.")

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
        ## possibly, transform relation to obtain a poset,
        ## and change main title if not specified.
        if (relation_is_linear_order(x[[i]])) {
            if (is.null(main[i]) || is.na(main[i]))
                main[i] <- "Linear Order"
        } else if (relation_is_partial_order(x[[i]])) {
            if (is.null(main[i]) || is.na(main[i]))
                main[i] <- "Partial Order"
        } else if(relation_is_weak_order(x[[i]])) {
            ## If x is a preference, use dual instead.
            x[[i]] <- dual(x[[i]])
            if (is.null(main[i]) || is.na(main[i]))
                main[i] <- "Weak Order"
        } else { ## extract asymmetric part
            x[[i]] <- x[[i]] & dual(x[[i]])
            if (is.null(main[i]) || is.na(main[i]))
                main[i] <- "Strict Preference Part"
        }

        ## Compute transitive reduction to avoid cluttered graph and
        ## extract incidence.
        I <- unclass(relation_incidence(transitive_reduction(x[[i]])))

        ## Perform reflexive reduction.
        diag(I) <- 0

        dimnames(I) <- lapply(labels(I), .make_unique_labels)
        plot(as(I, "graphNEL"), attrs = attrs[[i]], main = main[i], ...)
    }
}

####### plot incidence matrices

plot.relation_incidence <-
function(x, oma = c(3, 3, 1, 1), ...)
{
    o <- rev(order(colSums(x, na.rm = TRUE)))
    parhold <- par(no.readonly = TRUE)
    on.exit(par(parhold))
    par(las = 2L, oma = oma)
    x <- x[o, o]
    image.default(1:dim(x)[2L], 1:dim(x)[1L], t(x)[, dim(x)[1L]:1L],
                  axes = FALSE, col = gray(c(1, 0.5)),
                  xlab = "", ylab = "", ...)
    axis(1L, at = 1:dim(x)[1L], labels = rev(labels(x)[[1L]]))
    axis(2L, at = 1:dim(x)[2L], labels = rev(labels(x)[[2L]]))
}

### plot rankings

plot.ranking <-
function(x, ...)
    plot(as.relation(x), ...)


### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "### [*]+" ***
### End: ***

