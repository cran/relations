## Trivial duality infrastructure.

dual <-
function(x, ...)
    UseMethod("dual")

dual.relation <-
function(x, ...)
{
    if(!(is.relation(x) && relation_is_endorelation(x)))
        stop("Argument 'x' must be an endorelation.")
    I <- .incidence(x)
    ## Predicates for the dual relation of R can be inferred from those
    ## of R, see e.g. Fodor & Roubens, Table 2.2, page 41.
    db <- c(is_reflexive = "is_irreflexive",
            is_irreflexive = "is_reflexive",
            is_symmetric = "is_symmetric",
            is_antisymmetric = "is_complete",
            is_complete = "is_antisymmetric",
            is_asymmetric = "is_strongly_complete",
            is_strongly_complete = "is_asymmetric",
            is_transitive = "is_negatively_transitive",
            is_negatively_transitive = "is_transitive",
            is_Ferrers = "is_Ferrers",
            is_semitransitive = "is_semitransitive")
    predicates <-
        names(Filter(function(e) identical(e, TRUE),
                     relation_properties(x)[names(db)]))
    meta <- c(list(is_endorelation = TRUE),
              structure(as.list(rep.int(TRUE, length(predicates))),
                        names = db[predicates]))
    .make_relation_from_domain_and_incidence(.domain(x), .N.(t(I)), meta)
}
