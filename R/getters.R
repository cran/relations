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
        s <- relation_scores(x)
        ids <- match(s, sort(unique(s)))
        names(ids) <- names(s)
        ids
    }
    else if(relation_is_equivalence(x))
        get_class_ids_from_incidence(relation_incidence(x))
    else
        stop("Can only determine class ids for equivalences and weak orders.")
}

get_class_ids_from_incidence <-
function(x)
{
    ## Ugly ...
    y <- integer(nrow(x))
    c <- 1
    pos <- seq_along(y)
    while(length(pos)) {
        ind <- x[pos[1L], pos] == 1
        y[pos[ind]] <- c
        pos <- pos[!ind]
        c <- c + 1
    }
    names(y) <- rownames(x)
    y
}

## <FIXME>
## Are these still needed?
## If so, why not use relation_scores()?

get_order_from_incidence <-
function(x)
    labels(x)[[1L]][order(rowSums(x), decreasing = TRUE)]

get_preference_from_incidence <-
function(x)
    labels(x)[[1L]][order(rowSums(x))]

## </FIXME>
