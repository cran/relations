## Predicates.

## Note that we strongly prefer to use is_FOO rather than is.foo: we use
## this for distinguishing between class predicates (is.foo) and others.
## Alternatively, we could use a single predicate function
##   relation_test(x, predicate)
## The relation_is_${predicate}() approach has the advantage that some
## of these functions could be made generic eventually: provided we
## allow for relations without explicit coercion, we can simplify some
## of the tests (non-ordered factors given equivalence relations, etc).

### * Arity predicates

relation_is_binary <-
function(x)
{
    if(!is.relation(x))
        stop("Argument 'x' must be a relation.")
    .arity(x) == 2
}
relation_is_ternary <-
function(x)
{
    if(!is.relation(x))
        stop("Argument 'x' must be a relation.")
    .arity(x) == 3
}
relation_is_quaternary <-
function(x)
{
    if(!is.relation(x))
        stop("Argument 'x' must be a relation.")
    .arity(x) == 4
}

### * Predicates for general binary relations

## Note that typically relations will have arity metadata, so binarity
## can be checked without computing incidences.

relation_is_left_total <-
function(x)
{
    if(!is.relation(x))
        stop("Argument 'x' must be a relation.")
    if(!relation_is_binary(x)) return(FALSE)
    all(rowSums(relation_incidence(x)) >= 1)
}

relation_is_right_total <-
function(x)
{
    if(!is.relation(x))
        stop("Argument 'x' must be a relation.")
    if(!relation_is_binary(x)) return(FALSE)
    all(colSums(relation_incidence(x)) >= 1)
}

relation_is_surjective <- relation_is_right_total

relation_is_functional <- 
function(x)
{
    if(!is.relation(x))
        stop("Argument 'x' must be a relation.")
    if(!relation_is_binary(x)) return(FALSE)
    all(rowSums(relation_incidence(x)) <= 1)
}

relation_is_injective <-
function(x)
{
    if(!is.relation(x))
        stop("Argument 'x' must be a relation.")
    if(!relation_is_binary(x)) return(FALSE)
    all(colSums(relation_incidence(x)) <= 1)
}

relation_is_bijective <-
function(x)
{
    if(!is.relation(x))
        stop("Argument 'x' must be a relation.")
    if(!relation_is_binary(x)) return(FALSE)
    (all(rowSums(relation_incidence(x)) == 1)
     && all(colSums(relation_incidence(x)) == 1))
}

### * Endorelations and predicates of such

relation_is_endorelation <-
function(x)
{
    if(!is.relation(x))
        stop("Argument 'x' must be a relation.")
    if(.has_property(x, "is_endorelation"))
        return(.get_property(x, "is_endorelation"))
    .relation_is_endorelation_using_incidence(relation_incidence(x))
}

.relation_is_endorelation_using_incidence <-
function(x)
{
    ## For internal purposes only: used to avoid the need for possibly
    ## computing incidences at least twice in some of the code below
    ## (not that much of an issue as long as incidences are the only
    ## possible representation).

    ## Need some heuristic to determine whether we have an endorelation
    ## or not.  Idea: assume yes if nrow = ncol and either there are no
    ## dimnames (argh) or rownames and colnames are identical (better);
    ## otherwise, assume no.
    (is.matrix(x)
     && (nrow(x) == ncol(x))
     && identical(rownames(x), colnames(x)))
}
    
## <NOTE>
## Sometimes "total" is used synonymously to "complete".
## http://en.wikipedia.org/wiki/Binary_relation has two different usages
## of "total" ...
## </NOTE>
relation_is_complete <-
function(x)
{
    if(!is.relation(x))
        stop("Argument 'x' must be a relation.")
    if(.has_property(x, "is_complete"))
        return(.get_property(x, "is_complete"))
    x <- relation_incidence(x)
    if(!.relation_is_endorelation_using_incidence(x)) return(FALSE)
    all((x + t(x))[row(x) != col(x)] >= 1)
}

relation_is_reflexive <-
function(x)
{
    if(!is.relation(x))
        stop("Argument 'x' must be a relation.")
    if(.has_property(x, "is_reflexive"))
        return(.get_property(x, "is_reflexive"))
    x <- relation_incidence(x)
    if(!.relation_is_endorelation_using_incidence(x)) return(FALSE)
    all(diag(x) == 1)
}

relation_is_irreflexive <-
function(x)
{
    if(!is.relation(x))
        stop("Argument 'x' must be a relation.")
    x <- relation_incidence(x)
    if(!.relation_is_endorelation_using_incidence(x)) return(FALSE)
    all(diag(x) == 0)
}

relation_is_coreflexive <-
function(x)
{
    if(!is.relation(x))
        stop("Argument 'x' must be a relation.")
    x <- relation_incidence(x)
    if(!.relation_is_endorelation_using_incidence(x)) return(FALSE)
    all(x[row(x) != col(x)] == 0)
}

relation_is_symmetric <-
function(x)
{
    if(!is.relation(x))
        stop("Argument 'x' must be a relation.")
    if(.has_property(x, "is_symmetric"))
        return(.get_property(x, "is_symmetric"))
    x <- relation_incidence(x)
    if(!.relation_is_endorelation_using_incidence(x)) return(FALSE)
    all(x == t(x))
}

relation_is_asymmetric <-
function(x)
{
    if(!is.relation(x))
        stop("Argument 'x' must be a relation.")
    x <- relation_incidence(x)
    if(!.relation_is_endorelation_using_incidence(x)) return(FALSE)
    all((x & t(x)) == 0)
}

relation_is_antisymmetric <-
function(x)
{
    if(!is.relation(x))
        stop("Argument 'x' must be a relation.")
    if(.has_property(x, "is_antisymmetric"))
        return(.get_property(x, "is_antisymmetric"))
    x <- relation_incidence(x)
    if(!.relation_is_endorelation_using_incidence(x)) return(FALSE)
    all((x & t(x))[row(x) != col(x)] == 0)
}

relation_is_transitive <-
function(x)
{
    if(!is.relation(x))
        stop("Argument 'x' must be a relation.")
    if(.has_property(x, "is_transitive"))
        return(.get_property(x, "is_transitive"))
    x <- relation_incidence(x)
    if(!.relation_is_endorelation_using_incidence(x)) return(FALSE)
    all(((x %*% x) > 0) <= x)
}

## <FIXME>
## Add predicates for the following:
## Trichotomous: exactly one of xRy, yRx, or x=y holds.
## Euclidean: xRy & xRz => yRz.
## </FIXME>

## And now combine:
## Of course, these could be made more efficient by doing all
## computations [on incidences] just once ...

relation_is_equivalence <-
function(x)
    (relation_is_endorelation(x)
     && relation_is_reflexive(x)
     && relation_is_symmetric(x)
     && relation_is_transitive(x))

relation_is_weak_order <- relation_is_preference <-
function(x)
    (relation_is_endorelation(x)
     && relation_is_complete(x)
     && relation_is_transitive(x))

relation_is_preorder <- relation_is_quasiorder <-
function(x)
    (relation_is_endorelation(x)
     && relation_is_reflexive(x)
     && relation_is_transitive(x))

relation_is_partial_order <-
function(x)
    (relation_is_endorelation(x)
     && relation_is_reflexive(x)
     && relation_is_antisymmetric(x)
     && relation_is_transitive(x))

relation_is_linear_order <-
function(x)
    relation_is_partial_order(x) && relation_is_complete(x)

relation_is_tournament <-
function(x)    
    (relation_is_endorelation(x)
     && relation_is_complete(x)
     && relation_is_antisymmetric(x))

### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "### [*]+" ***
### End: ***
