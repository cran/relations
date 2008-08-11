relation_class_ids <-
function(x)
{
    if(!is.relation(x))
        stop("Argument 'x' must be a relation.")
    if(relation_is_weak_order(x)) {
        ## Get the class ids of the corresponding indifference relation.
        ## One possibility:
        ##   I <- relation_incidence(x)
        ##   get_class_ids_from_incidence(I & t(I))
        ## But this is faster and also gives the class ids in the
        ## natural order of the indifference classes.
        s <- relation_scores(x, "ranks", decreasing = FALSE)
        ids <- match(s, sort(unique(s)))
        names(ids) <- names(s)
        ids
    }
    else if(relation_is_equivalence(x))
        get_class_ids_from_incidence(relation_incidence(x))
    else
        stop("Can only determine class ids for equivalences and weak orders.")
}

relation_classes <-
function(x)
{
    ids <- relation_class_ids(x)
    out <- split(seq_along(ids), ids)
    class(out) <- c("relation_classes_of_objects")
    attr(out, "labels") <- names(ids)
    out
}

print.relation_classes_of_objects <-
function(x, ...)
{
    labels <- attr(x, "labels")
    y <- lapply(x, function(i) paste(labels[i], collapse = ", "))
    writeLines(formatDL(rev(names(x)), sprintf("{%s}", rev(unlist(y))),
                        style = "list", ...))
    invisible(x)
}

get_class_ids_from_incidence <-
function(x)
{
    ## Ugly ...
    y <- integer(nrow(x))
    c <- 1L
    pos <- seq_along(y)
    while(length(pos)) {
        ind <- x[pos[1L], pos] == 1
        y[pos[ind]] <- c
        pos <- pos[!ind]
        c <- c + 1L
    }
    names(y) <- rownames(x)
    y
}

## Let R be an endorelation.
## An element x is minimal if there is no "smaller" one, i.e.:
##   There is no y != x with y R x
## An element x is a first element if it is "not larger" than any
## other element, i.e.:
##   x R y for all y != x
## An element x is maximal if there is no "larger" one, i.e.:
##   There is no y != x with x R y
## An element x is a last element if it is "not smaller" than any
## other element, i.e.:
##   y R x for all y != x

## Note that sets cannot directly be indexed positionally.
## (Well, as of 2008-08-08 there is sets:::.set_subset() ...)

relation_minimal_elements <-
function(x, na.rm = TRUE)
{
    if(!relation_is_endorelation(x))
        stop("Argument 'x' must be an endorelation.")
    X <- as.list(relation_domain(x)[[1L]])
    I <- .incidence(x)
    diag(I) <- 0
    as.set(X[colSums(I, na.rm = na.rm) == 0])
}

relation_first_elements <-
function(x, na.rm = FALSE)
{
    if(!relation_is_endorelation(x))
        stop("Argument 'x' must be an endorelation.")
    X <- as.list(relation_domain(x)[[1L]])
    I <- .incidence(x)
    diag(I) <- 1
    ind <- rowSums(I != 1, na.rm = na.rm) == 0
    as.set(X[!is.na(ind) & ind])
}

relation_last_elements <-
function(x, na.rm = FALSE)
{
    if(!relation_is_endorelation(x))
        stop("Argument 'x' must be an endorelation.")
    X <- as.list(relation_domain(x)[[1L]])
    I <- .incidence(x)
    diag(I) <- 1
    ind <- colSums(I != 1, na.rm = na.rm) == 0
    as.set(X[!is.na(ind) & ind])
}

relation_maximal_elements <-
function(x, na.rm = TRUE)
{
    if(!relation_is_endorelation(x))
        stop("Argument 'x' must be an endorelation.")
    X <- as.list(relation_domain(x)[[1L]])
    I <- .incidence(x)
    diag(I) <- 0
    as.set(X[rowSums(I, na.rm = na.rm) == 0])
}

## Let R be an endorelation.
## An element y is a *successor* of an element x if y != x, x R y
## and for all z != y with x R z we have y R z.
## <FIXME>
## One should verify that this is really the "official" definition.
## </FIXME>

relation_successors <-
function(x, e = NULL)
{
    if(!relation_is_endorelation(x))
        stop("Argument 'x' must be an endorelation.")

    X <- as.list(relation_domain(x)[[1L]])
    ## Need to find the positions of e in X.
    if(is.null(e)) {
        pos <- seq_along(X)
    } else {
        pos <- sets:::.exact_match(e, X)
        if(any(is.na(pos)))
            stop("Elements of 'e' must be contained in the domain components.")
    }
    ## Argh, terminology is really a nuisance.  If we had gone with the
    ## domain/codomain (domain/range) terminology commonly employed for
    ## endorelations, we could say that X is the domain of the relation,
    ## but then what is the tuple (X, X) called?  (And what is used in
    ## the general k-ary case?)  [Wikipedia says that X_1, ..., X_k are
    ## the domains of the relation.]
    I <- relation_incidence(x)
    class(I) <- "matrix"
    diag(I) <- 0
    out <- lapply(pos,
                  function(p) {
                      candidates <- which(I[p, ] == 1)
                      ## Need to find those candidates y for which there
                      ## is no z != y with z R y.
                      ind <- candidates[colSums(I[candidates,
                                                  candidates,
                                                  drop = FALSE]) == 0]
                      as.set(X[ind])
                  })
    names(out) <- rownames(I)[pos]
    out
}

## Should similarly have something for finding precursors.
