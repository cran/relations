## * constructor

relation_multiset <-
function(..., list = NULL, times = 1L)
{
    mset <- unique(c(list(...), list))
    if(!length(mset))
        stop("Cannot create empty multisets.")
    .check_domains(mset)
    .make_relation_multiset_from_gset(
        gset(mset, rep(times, length.out = length(mset))))
}

.make_relation_multiset_from_gset <-
function(x)
    structure(x, class = unique(c("relation_multiset", class(x))))

.check_domains <-
function(x)
{
    if(length(as.list(x)) > 1L) {
        domains <- lapply(x, .domain)
        ## Check whether all domains are the same.
        if(!all(sapply(domains, .domain_is_equal, domains[[1L]]))) {
            ## Only need to test domains[-1L] of course, but it is not
            ## clear whether doing so really saves time.
            stop("All relations must have the same domain.")
        }
    }
}

## * is. method

is.relation_multiset <-
function(x)
    inherits(x, "relation_multiset")

### * as.foo methods

as.relation_multiset <-
function(x)
    UseMethod("as.relation_multiset")

as.relation_multiset.default <-
function(x)
    relation_multiset(x)

as.relation_multiset.relation_ensemble <-
function(x)
    .make_relation_multiset_from_gset(as.gset(x))

as.relation_multiset.relation_multiset <- identity

as.relation_ensemble.relation_multiset <-
function(x)
    relation_ensemble(list = rep.int(unclass(x), gset_memberships(x)))

### * print method

print.relation_multiset <-
function(x, ...)
{
    len <- length(as.list(x))
    if(len)
        writeLines(sprintf(ngettext(len,
                                    "A multiset of %d relation of size %s.",
                                    "A multiset of %d relations of size %s."),
                           len,
                           paste(.size(x), collapse = " x ")))
    else
        writeLines(gettext("An empty relation multiset."))
    invisible(x)
}

### c method

c.relation_multiset <-
function(...)
{
    ret <- do.call(gset_union, lapply(list(...), as.relation_multiset))
    .check_domains(ret)
    .make_relation_multiset_from_gset(ret)
}

### Summary methods

Summary.relation_multiset <-
function(..., na.rm = FALSE)
{
    ok <- switch(.Generic, max = , min = , range = TRUE, FALSE)
    if(!ok)
        stop(gettextf("Generic '%s' not defined for \"%s\" objects.",
                      .Generic, .Class))

    ## Combine the given multisets, and convert to ensemble
    ## by discarding multiplicity
#    x <- relation_ensemble(list = as.list(do.call(c, list(...))))
    x <- do.call(relation_ensemble, c(...))

    switch(.Generic,
           "min" = .relation_meet(x),
           "max" = .relation_join(x),
           "range" = {
               relation_ensemble(min = .relation_meet(x),
                                 max = .relation_join(x))
           })
}

