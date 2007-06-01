### Relation metrics (dissimilarities).

## We use the CLUE approach: have a front-end
##   relation_dissimilarity <-
##       function(x, y = NULL, method = "symdiff", ...)
## which returns dissimilarities or cross-dissimilarities, eventually
## with classes and subscript methods as in CLUE, and a mechanism for 
## registering dissimilarity methods.

relation_dissimilarity <-
function(x, y = NULL, method = "symdiff", ...)
{
    x <- as.relation_ensemble(x)

    ## Be nice.
    if(is.character(y) || is.function(y)) {
        method <- y
        y <- NULL
    }

    known_methods <-
        list("symdiff" = c(".relation_dissimilarity_symdiff",
             "symmetric difference distance"))
    if(is.character(method)) {
        ## Hopefully of length one, add some tests eventually ...
        if(is.na(ind <- pmatch(method, names(known_methods))))
            stop(gettextf("Method '%s' is not a valid dissimilarity method.",
                          method))
        method <- get(known_methods[[ind]][1L])
        method_name <- known_methods[[ind]][2L]
    }
    else if(is.function(method))
        method_name <- "user-defined method"
    else
        stop("Invalid argument 'method'.")

    if(!is.null(y)) {
        y <- as.relation_ensemble(y)
        D <- relation_domain(x)
        if(!.domain_is_equal(relation_domain(y), D))
            stop("All relations must have the same domain.")
        y <- .canonicalize_relation_ensemble(y, D)
        ## Build a cross-proximity object of cross-dissimilarities.
        ## <FIXME>
        ## Not yet: if we don't want to require clue, all we can return
        ## is a matrix ...
        d <- matrix(0, length(x), length(y))
        for(j in seq_along(y))
            d[, j] <- sapply(x, method, y[[j]])
        dimnames(d) <- list(names(x), names(y))
        return(d)
        ## </FIXME>
    }

    ## Otherwise, build a proximity object of dissimilarities.
    ## <FIXME>
    ## Not yet: if we don't want to require clue, all we can return
    ## is a dist object ...
    n <- length(x)
    d <- vector("list", length = n - 1)
    ind <- seq_len(n)
    while(length(ind) > 1L) {
        j <- ind[1L]
        ind <- ind[-1L]
        d[[j]] <- sapply(x[ind], method, x[[j]])
    }
    ## Grr ... see clue:::.dist_from_vector().
    structure(unlist(d), Size = n, Labels = names(x),
              class = "dist")
}

.relation_dissimilarity_symdiff <-
function(x, y)
    sum(abs(relation_incidence(x) - relation_incidence(y)))
