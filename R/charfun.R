## Characteristic functions.

### * relation_charfun

relation_charfun <-
function(x, components = FALSE)
{
    if(!is.relation(x))
        stop("Argument 'x' must be a relation.")
    local({I <- relation_incidence(x)
           D <- relation_domain(x)
           a <- .arity(x)
           if (!components)
               function(...) {
                   args <- list(...)
                   if (a == 2) {
                       ## recycle for binary relations
                       maxlen <- max(sapply(args, length))
                       args <- lapply(args, rep, length.out = maxlen)
                   }
                   if(length(args) != a)
                       stop("Wrong number of arguments.")
                   t <- .split_into_components(do.call("cbind", args))
                   as.logical(I[rbind(mapply(match, t, D))])
               }
           else
               function(t) {
                   t <- .split_into_components(t)
                   as.logical(I[rbind(mapply(match, t, D))])
               }
       })
}

### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "### [*]+" ***
### End: ***
