## Relation ensembls, cf. the cluster ensembles in CLUE.

## Note that ensembles are not implemented as containers.
## Hence, properties (if at all) are stored as attributes with "hidden"
## (leading dot) names.
##
## <NOTE>
## It is possible to conveniently iterate through containers via
## lapply() [CLUE style] provided there is an as.list() method for the
## container class.  E.g.,
##   as.list.FOO <- function(x, ...) x$ZZZ
##   FOO <- function(...)
##     structure(list(ZZZ = list(...)), class = "FOO")
##   lapply(FOO(1, 2), function(x) x)
## gives list(1, 2).
## We could thus implement ensembles as containers rather than list
## "structures".  But there is little benefit (if at all) in this as
## long as we only have a trivial "list of relations" representation.
## As all relations in an ensemble must have the same domain, we could
## represent ensembles by a single domain and a list of (sparse or
## dense) incidences (again with unclear benefits apart from space
## savings).  If we added alternative representations, realizing
## relation ensembles via containers would be much more attractive.
## </NOTE>
##
## Note also that (unlike CLUE), we do not allow the creation of empty
## relation ensembles, so that we always have the domain information.
## It is possible, however, to obtain empty ensembles by subscripting,
## see below.

### * relation_ensemble

relation_ensemble <-
function(..., list = NULL)
{
    relations <- lapply(c(list(...), list), as.relation)
    if(!length(relations))
        stop("Empty relation ensembles cannot be created.")
    domains <- lapply(relations, relation_domain)
    if(!all(sapply(domains, .domain_is_equal, domains[[1L]])))
        stop("All relations must have the same domain.")
    relations <-
        lapply(relations, .canonicalize_relation, domains[[1L]])
     .make_relation_ensemble_from_list_and_meta(relations)
}

.make_relation_ensemble_from_list_and_meta <-
function(x, meta = NULL)
{
    if(is.null(meta))
        meta <- .make_ensemble_meta_from_element(x[[1L]])
    structure(x, .Meta = meta, class = c("relation_ensemble", "tuple"))
}

.make_ensemble_meta_from_element <-
function(x)
    list(arity = .arity(x), domain = .domain(x), size = .size(x))

### * is.relation_ensemble

is.relation_ensemble <-
function(x)
    inherits(x, "relation_ensemble")

### * as.relation_ensemble

as.relation_ensemble <-
function(x)
    UseMethod("as.relation_ensemble")
as.relation_ensemble.data.frame <-
function(x)
{
    nms <- rownames(x)
    relations <- lapply(x,
                        function(u) {
                            names(u) <- nms
                            as.relation(u)
                        })
    .make_relation_ensemble_from_list_and_meta(relations)
}
as.relation_ensemble.default <-
function(x)
    relation_ensemble(x)
as.relation_ensemble.relation_ensemble <- identity

### * Methods

"[.relation_ensemble" <-
function(x, i)
{
    ## Make subscripting empty ensembles a noop.
    if(length(x) == 0L) return(x)
    ## What should we do when subscripting out of bounds?
    ## For lists, this gives NULL elements.  We could convert them to
    ## null relations (all incidences zero).  For now, we remove them.
    y <- NextMethod("[")
    .make_relation_ensemble_from_list_and_meta(y[!sapply(y, is.null)],
                                               attr(x, ".Meta"))
}

## <NOTE>
## We could have a [[ method which ensures that the ensemble metadata
## (arity, domain and size) get tucked on.  But should we?
# </NOTE>

c.relation_ensemble <-
function(..., recursive = FALSE)
{
    relations <- unlist(lapply(list(...), as.relation_ensemble),
                        recursive = FALSE)
    relation_ensemble(list = relations)
}

t.relation_ensemble <-
function(x)
    .make_relation_ensemble_from_list_and_meta(lapply(x, t),
                                               attr(x, ".Meta"))

rep.relation_ensemble <-
function(x, times, ...)
    .make_relation_ensemble_from_list_and_meta(NextMethod("rep"),
                                               attr(x, ".Meta"))

unique.relation_ensemble <-
function(x, incomparables = FALSE, ...)
    .make_relation_ensemble_from_list_and_meta(NextMethod("unique"),
                                               attr(x, ".Meta"))

print.relation_ensemble <-
function(x, ...)
{
    len <- length(x)
    if(len > 0L)
        writeLines(sprintf(ngettext(len,
                                    "An ensemble of %d relation of size %s.",
                                    "An ensemble of %d relations of size %s."),
                           len,
                           paste(.size(x), collapse = " x ")))
    else
        writeLines(gettext("An empty relation ensemble."))
    invisible(x)
}

Ops.relation_ensemble <-
function(e1, e2)
{
    ## Relation ensembles are of the tuple type, so that comparisons etc
    ## are performed elementwise.  It might be useful to recycle ...
    if(missing(e2)) {
        ## Could have unary '+', '-' and '!'.
        ## Only complement makes sense.  Could make the other a no-op.
        if(!(as.character(.Generic) %in% "!"))
            stop(gettextf("Unary '%s' not defined for \"%s\" objects.",
                          .Generic, .Class))
        ## Could also do relation_ensemble(list = NextMethod()).
        return(.make_relation_ensemble_from_list_and_meta(lapply(e1,
                                                                 .Generic),
                                                          attr(e1,
                                                               ".Meta")))
    }
    ## Otherwise, keep things simple.  Dispatch to tuple method:
    out <- NextMethod()
    ## And massage the return type:
    if(as.character(.Generic) %in% c("<", "<=", ">", ">=", "==", "!="))
        as.logical(out)
    else
        relation_ensemble(list = out)
    ## (Could make the last more efficient, of course.)
}

Summary.relation_ensemble <-
function(..., na.rm = FALSE)
{
    ok <- switch(.Generic, max = , min = , range = TRUE, FALSE)
    if(!ok)
        stop(gettextf("Generic '%s' not defined for \"%s\" objects.",
                      .Generic, .Class))
    args <- list(...)
    ## Combine the given relation ensembles.
    x <- do.call(c, args)
    switch(.Generic,
           "min" = .relation_meet(x),
           "max" = .relation_join(x),
           "range" = {
               relation_ensemble(min = .relation_meet(x),
                                 max = .relation_join(x))
           })
}

as.list.relation_ensemble <-
function(x, ...)
{
    attributes(x) <- NULL
    x
}

### * Utilities

### ** Meet and Join.

.relation_meet <-
function(x)
{
    if(!is.relation_ensemble(x))
        stop("Argument 'x' must be a relation ensemble.")
    I <- do.call(pmin, lapply(x, relation_incidence))
    .make_relation_from_domain_and_incidence(.domain(x), I)
}

.relation_join <-
function(x)
{
    if(!is.relation_ensemble(x))
        stop("Argument 'x' must be a relation ensemble.")
    I <- do.call(pmax, lapply(x, relation_incidence))
    .make_relation_from_domain_and_incidence(.domain(x), I)
}

### ** .canonicalize_relation_ensemble

.canonicalize_relation_ensemble <-
function(x, D)
{
    ## Canonicalize a relation ensemble to have its domain elements use
    ## the same *internal* order as the elements of D (assuming that the
    ## domain of the ensemble is known to equal D in the sense that the
    ## respective elements are the same sets).
    ##
    ## Note that relation ensembles are always canonicalized.

    if(length(x) == 0L) return(x)
    pos <- .match_domain_components(lapply(D, as.set),
                                    relation_domain(x))
    if(!any(sapply(pos, is.unsorted))) return(x)
    x <- lapply(x, .canonicalize_relation, D, pos)
    size <- lapply(D, length)
    meta <- list(arity = length(size), domain = D, size = size)
    .make_relation_ensemble_from_list_and_meta(x, meta)
}

### ** .is_ensemble_of_endorelations

.is_ensemble_of_endorelations <-
function(x)
{
    ## Check whether we have an ensemble of endorelations (assuming that
    ## ensembles are known to have identical identical domains).
    relation_is_endorelation(x[[1L]])
}

### ** .is_ensemble_of_crisp_relations

.is_ensemble_of_crisp_relations <-
function(x)
    all(sapply(x, relation_is_crisp))


### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "### [*]+" ***
### End: ***

