############
### Sets ###
############

### Basic stuff (constructors, print/summary methods)
set <-
function(...)
    .make_set_from_list(.set_unique(list(...)))

print.set <-
function(x, ...)
{
    writeLines(strwrap(format(x), exdent = 1L))
    invisible(x)
}

summary.set <-
function(object, ...)
{
    len <- length(object)
    out <- if (len == 0L)
        gettext("The empty set.")
    else if (len == 1L)
        gettext("A set with 1 element.")
    else
        gettextf("A set with %d elements.", len)
    structure(out, class = "summary.set")
}

print.summary.set <-
function(x, ...)
{
    writeLines(x)
    invisible(x)
}

### converters

as.set <-
function(x)
    UseMethod("as.set")

as.set.default <-
function(x)
    stop("Not implemented.")

as.set.set <- .identity

as.set.tuple <-
function(x)
    do.call(set, unclass(x))

as.set.numeric <-
as.set.factor <-
as.set.character <-
as.set.integer <-
as.set.ordered <-
as.set.logical <-
function(x)
    .set_unique(as.list(x))

as.set.list <-
function(x)
    .make_set_from_list(x)

as.set.data.frame <-
function(x)
    .make_set_from_list(lapply(split(x, seq_len(nrow(x))), as.tuple))

as.set.relation_graph <-
function(x)
    structure(x, class = "set")

format.set <-
function(x, ...) {
    .format_set_or_tuple(x, "{", "}")
}
  
as.relation.set <-
function(x)
    relation(graph = x)

as.list.set <-
function(x, ...)
    unclass(x)

## predicates

is.set <-
function(x)
    inherits(x, "set")

set_is_empty <-
function(x)
{
    if(is.set(x))
        length(x) < 1L
    else
        sapply(x, length) < 1L
}

set_is_subset <-
function(a, b)
{
    .help <- function(a,b) set_is_empty(set_complement(b, a))
    if(is.set(a))
        .help(a, b)
    else
        Vectorize(.help)(a, b)
}

set_is_proper_subset <-
function(a, b)
{
    set_is_subset(a, b) &
    if(is.set(a))
        length(a) != length(b)
    else
        sapply(a, length) != sapply(b, length)
}

set_is_equal <-
function(a, b)
{
    .help <- function(a, b)
        ((length(a) == length(b))
         && (length(set_intersection(a,b)) == length(a)))
    if(is.set(a))
        .help(a,b)
    else
        Vectorize(.help)(a,b)
}

"%e%" <-
set_is_element <-
function(e, b)
{
    if(set_is_empty(b))
        return(FALSE)
    if(is.set(e)) {
        if(!any(sapply(b, is.set)))
            stop("Set to look into does not contain any set.")
        e <- list(e)
    }
    if(is.tuple(e))
        e <- list(e)
    e %in% b
}

### methods

c.set <-
set_union <-
function(...)
    .set_unique(do.call(c, lapply(list(...), unclass)))

set_intersection <-
function(...)
{
    len <- length(l <- list(...))
    if(len < 1L)
        set()
    else if(len < 2L)
        l[[1L]]
    else if(len < 3L)
        .make_set_from_list(l[[2L]][unique(na.omit(match(l[[1L]],
                                                         l[[2L]])))])
    else
        do.call(Recall, c(l[1L], list(do.call(Recall, l[-1L]))))
}

set_complement <-
function(a, b)
{
    ind <- unique(na.omit(match(a, b)))
    .make_set_from_list(if(length(ind)) b[-ind] else b)
}

"%D%" <-
set_symdiff <-
function(...)
    set_complement(set_intersection(...), set_union(...))

set_power <-
function(x)
    set_union(set(set()),
              as.set(unlist(lapply(seq_along(x),
                                   function(i)
                                   apply(combn(x, i), 2L, as.set)
                                   ),
                            recursive = FALSE)
                     )
              )

set_cartesian <-
function(...)
{
    if(nargs() < 2L)
        return(..1)
    l <- list(...)
    if(!all(len <- sapply(l, length)))
        return(set())
    .make_set_of_tuples_from_list_of_lists(.cartesian_product(l))
}

set_combn <-
function(x, m)
{
    if(m == 0)
        set()
    else
        do.call(set, apply(combn(x, m), 2L, as.set))
}

set_outer <-
function(X, Y, FUN = "*", ..., SIMPLIFY = TRUE)
{
    ## convenience
    nx <- deparse(substitute(X))
    if(missing(Y)) {
        Y <- X
        ny <- nx
    } else if(is.function(Y) || is.character(Y)) {
        FUN <- Y
        Y <- X
        ny <- nx
    } else ny <- deparse(substitute(Y))

    FUN <- match.fun(FUN)
  
    ## loop
    xrep <- rep(unclass(X), times = (ylen <- length(Y)))
    yrep <- rep(unclass(Y), each = (xlen <- length(X)))
    ret <- mapply(FUN, xrep, yrep, MoreArgs = list(...), SIMPLIFY = FALSE)
    
    ## simplify if sensible
    if(SIMPLIFY && all(sapply(ret, is.atomic)))
        ret <- unlist(ret, recursive = FALSE)
    
    ## make matrix
    dim(ret) <- c(xlen, ylen)
    dimnames(ret) <- list(LABELS(X), LABELS(Y))
    ret
}

### operators

Ops.set <-
function(e1, e2)
{
    if(nargs() == 1L)
        stop("Unary operators not defined for \"set\" objects.")
  
    if(!(as.character(.Generic)
         %in% c("<", "<=", ">", ">=", "==", "!=",
                "&", "|", "*", "+", "-", "^")))
        stop(gettextf("Generic '%s' not defined for \"%s\" objects.",
                      .Generic, .Class))

    if(as.character(.Generic) == "^") {
        if(is.set(e1) &&
            ((trunc(e2) != e2) || (e2 < 1L)))
            stop("Cartesian product only defined for positive integers.")
        if(is.set(e2) && (e1 != 2L))
            stop("Operator not defined.")
    }
    
    switch(.Generic,
           "+"  =,
           "|"  = set_union(e1, e2),
           "-"  = set_complement(e2, e1),
           "&"  = set_intersection(e1, e2),
           "*"  = set_cartesian(e1, e2),
           "<"  = set_is_proper_subset(e1, e2),
           "<=" = set_is_subset(e1, e2),
           ">"  = set_is_proper_subset(e2, e1),
           ">=" = set_is_subset(e2, e1),
           "==" = set_is_equal(e1, e2),
           "!=" = !set_is_equal(e1, e2),
           "^"  = {
               if(is.set(e2))
                   set_power(e2)
               else
                   do.call(set_cartesian, rep(list(e1), e2))}
           )
         
}

rep.set <-
function(x, ...)
    x                                   # for the time being ...

### internal stuff

.set_unique <-
function(x)
    .make_set_from_list(x[!duplicated(x)])

.make_set_from_list <-
function(x)
    structure(x, class = "set")

.format_set_or_tuple <-
function(x, left, right)
{
    nms <- names(x)
    names(x) <- NULL
    SEP <- rep("", length(x))
    if (!is.null(nms))
      SEP[nms != ""] <- " = "
    paste(left,
          if (length(x) > 0)
              paste(nms, SEP, LABELS(unclass(x)),
                    sep = "", collapse = ", "),
          right,
          sep = "")
  }
