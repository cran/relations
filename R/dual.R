## Trivial duality infrastructure.

dual <-
function(x, ...)
    UseMethod("dual")

## Ideally there would be a good reference which defines the dual as the
## negative of the inverse, as e.g. Fishburn uses "dual" differently.
## Possible references are
##   Ovchinnikov (2000)
##     <http://link.springer.com/chapter/10.1007%2F978-1-4615-4429-6_5>
## or chapter 2 in Fodor & Roubens (1994)
##     <http://www.springer.com/us/book/9780792331162>
## but there should really be better/earlier ones ...

dual.relation <-
function(x, ...)
{
    if(!relation_is_binary(x))
        stop("Argument 'x' must be a binary relation.")
    I <- .incidence(x)
    meta <- if(relation_is_endorelation(x)) {
        ## Predicates for the dual relation of an endorelation R can be
        ## inferred from those of R, see e.g. Fodor & Roubens, "Fuzzy
        ## Preference Modelling and Multicriteria Decision Support",
        ## Table 2.2, page 41.
        ## <http://www.springer.com/us/book/9780792331162>
        ## For valued relations, only the correspondencies
        ##   reflexive <-> irreflexive, symmetric <-> symmetric
        ## are always true: the others require a deMorgan triple of
        ## fuzzy connectives N/T/S.
        db <- c(is_reflexive = "is_irreflexive",
                is_irreflexive = "is_reflexive",
                is_symmetric = "is_symmetric")
        if(fuzzy_logic_predicates()$is_de_Morgan_triple) {
            db <-
                c(db,
                  is_antisymmetric = "is_complete",
                  is_complete = "is_antisymmetric",
                  is_asymmetric = "is_strongly_complete",
                  is_strongly_complete = "is_asymmetric",
                  is_transitive = "is_negatively_transitive",
                  is_negatively_transitive = "is_transitive",
                  is_Ferrers = "is_Ferrers",
                  is_semitransitive = "is_semitransitive"
                  )
        }
        predicates <-
            names(Filter(function(e) identical(e, TRUE),
                         relation_properties(x)[names(db)]))
        c(list(is_endorelation = TRUE),
          .structure(as.list(rep.int(TRUE, length(predicates))),
                     names = db[predicates]))
    } else NULL
    .make_relation_from_domain_and_incidence(.domain(x), .N.(t(I)), meta)
}
