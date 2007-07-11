relation_violations <-
function(x, family = c("T", "transitive",
                         "S", "symmetric",
                         "A", "antisymmetric",
                         "C", "complete",
                         "R", "reflexive"))
{
    if (!relation_is_endorelation(x))
        stop("Relation violations only defined for endorelations.")

    family <- match.arg(family)
    I <- .incidence(x)

    switch(family,
           T =, transitive = .non_transitivity(I),
           S =, symmetric = .non_symmetry(I),
           A =, antisymmetric = .non_antisymmetry(I),
           C =, complete = .non_completeness(I),
           R =, reflexive = .non_reflexivity(I)
           )
}

.non_transitivity <-
function(I)
    sum(sapply(seq_len(nrow(I)),
               function(j) outer(I[, j], I[j, ], "+") - I) > 1)

.non_symmetry <-
function(I)
    sum(I & I != t(I))

.non_antisymmetry <-
function(I)
{
    diag(I) <- 0
    sum(I & t(I)) / 2
}

.non_completeness <-
function(I)
{
    I <- 1 - I
    D <- diag(I)
    diag(I) <- 0
    sum(I & t(I)) / 2 + sum(D)
}

.non_reflexivity <-
function(I)
     sum(1 - diag(I))


