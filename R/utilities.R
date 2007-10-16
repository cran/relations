### * .domain_is_equal

.domain_is_equal <-
function(D1, D2)
    all(mapply(set_is_equal, lapply(D1, as.set), lapply(D2, as.set)))

### * .is_subsettable
### FIXME: not needed anymore?
#.is_subsettable <-
#function(x)
#    tryCatch({x[[1L]]; TRUE}, error = function(e) FALSE)

### * .is_valid_relation_domain

## <NOTE>
## We deliberately only check on lists and not on tuple class because we
## want to be "nice" and allow "rawer" specifications for convenience.
## </NOTE>

.is_valid_relation_domain <-
function(x)
{
    (is.list(x)
     && (length(x) != 0L)
     && all(sapply(x, length) > 0L))
}

### * .make_set_of_tuples_from_relation_graph_components

.make_set_of_tuples_from_relation_graph_components <-
sets:::.make_set_of_tuples_from_list_of_lists

### * .match_domain_components

.match_domain_components <-
function(x, y)
    mapply(match, x, y, SIMPLIFY = FALSE)

### * .negate

## <FIXME>
## Remove eventually once we have Negate() or equivalent in a released
## version of R.
.negate <-
function(f)
    function(...) ! match.fun(f)(...)
## </FIXME>

### * .offdiag

.offdiag <-
function(x)
    row(x) != col(x)

### * .relation_meta_db

.relation_meta_db <-
    list(E =
         list(is_endorelation = TRUE,
              is_reflexive = TRUE,
              is_symmetric = TRUE,
              is_transitive = TRUE),
         L =
         list(is_endorelation = TRUE,
              is_complete = TRUE,
              is_antisymmetric = TRUE,
              is_transitive = TRUE),
         O =
         list(is_endorelation = TRUE,
              is_antisymmetric = TRUE,
              is_transitive = TRUE),
         P =
         list(is_endorelation = TRUE,
              is_complete = TRUE,
              is_reflexive = TRUE,
              is_transitive = TRUE),
         T =
         list(is_endorelation = TRUE,
              is_complete = TRUE,
              is_irreflexive = TRUE,
              is_antisymmetric = TRUE,
              is_asymmetric = TRUE),
         C =
         list(is_endorelation = TRUE,
              is_complete = TRUE),
         A =
         list(is_endorelation = TRUE,
              is_antisymmetric = TRUE),
         S =
         list(is_endorelation = TRUE,
              is_symmetric = TRUE),
         M =
         list(is_endorelation = TRUE,
              is_complete = TRUE,
              is_reflexive = TRUE)
         )

### * .reorder_incidence

.reorder_incidence <-
function(I, pos)
    do.call("[", c(list(I), pos, list(drop = FALSE)))

### * .split_into_components

.split_into_components <-
function(x)
{
    ## If x is a matrix: make list of columns; otherwise, leave alone
    ## which is ok for at least lists and data frames.
    if(is.matrix(x))
        split(x, col(x))
    else
        x
}

### * .transform_factors_into_characters

## <FIXME>
## This seems necessary because match() incorrectly (?) handles
## comparisons of lists of characters and lists of factors
.transform_factors_into_characters <-
function(x)
    lapply(x,
           lapply,
           function(i) if(is.factor(i)) as.character(i) else i)
## </FIXME>

### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "### [*]+" ***
### End: ***
