### Relation incidences.

## * relation_incidence

relation_incidence <-
function(x)
{
    if(!is.relation(x))
        stop("Argument 'x' must be a relation.")
    I <- .incidence(.get_representation(x))
    structure(as.array(I),
              dimnames = lapply(relation_domain(x), LABELS),
              class = "relation_incidence")
}

.incidence <-
function(x)
    UseMethod(".incidence")
.incidence.relation <-
function(x)
    .get_property_from_object_or_representation(x, "incidence", .incidence)
.incidence.relation_by_domain_and_incidence <-
function(x)
    x$incidence

print.relation_incidence <-
function(x, ...)
{
    writeLines("Incidences:")
    print(array(as.vector(x),
                dim = dim(x),
                dimnames = dimnames(x)),
          ...)
    invisible(x)
}

### * relation_incidence<-

"relation_incidence<-" <-
function(x, value)
{
    if(!is.relation(x))
        stop("Argument 'x' must be a relation.")
    size <- .size(x)
    I <- as.array(value)
    if(length(size) != length(dim(I)))
        stop("Relation arity mismatch between 'x' and 'value'.")
    if(any(size != dim(I)))
        stop("Relation size mismatch between 'x' and 'value'.")
    .make_relation_from_domain_and_incidence(.domain(x), I)
}

### * .is_valid_relation_incidence

.is_valid_relation_incidence <-
function(x)
{
    if(length(x) == 0L) return(FALSE)
    x <- as.array(x)
    if(any(dim(x) == 0L)) return(FALSE)
    (is.logical(x)
     || (is.numeric(x)
         && all(x >= 0, na.rm = TRUE)
         && all(x <= 1, na.rm = TRUE)))
}

### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "### [*]+" ***
### End: ***
