relation_violations <-
function(x,
         property =
         c("complete", "match",
           "reflexive", "irreflexive", "coreflexive",
           "symmetric", "antisymmetric", "asymmetric",
           "transitive", "negatively_transitive", "Ferrers",
           "semitransitive",
           "trichotomous",
           "Euclidean"),
         tuples = FALSE)
{
    if (!relation_is_endorelation(x))
        stop("Relation violations only defined for endorelations.")

    property <- match.arg(property)
    I <- .incidence(x)

    if(!tuples) {
        do.call(sprintf(".amount_by_which_relation_is_not_%s", property),
                list(I))
    } else {
        ## First get a matrix of indices of violating tuples.
        ind <- do.call(sprintf(".tuples_for_which_relation_is_not_%s",
                               property),
                       list(I))
        if(!nrow(ind)) return(set())
        ## And construct a set of violating tuples from this.
        D <- rep.int(list(.get_elements_in_homorelation(x)), ncol(ind))
        as.set(apply(ind, 1L, function(i) as.tuple(Map(`[[`, D, i))))
    }
}

.amount_by_which_relation_is_not_complete <-
function(I)
{
    I <- .N.(I)
    diag(I) <- 0
    sum(.T.(I, t(I))) / 2
}

.tuples_for_which_relation_is_not_complete <-
function(I)
{
    I <- .N.(I)
    diag(I) <- 0
    ind <- which(.T.(I, t(I)) > 0, arr.ind = TRUE)
    ind[ind[, 1L] <= ind[, 2L], , drop = FALSE]
}

.amount_by_which_relation_is_not_match <-
function(I)
{
    I <- .N.(I)
    D <- diag(I)
    diag(I) <- 0
    sum(.T.(I, t(I))) / 2 + sum(D)
}

.tuples_for_which_relation_is_not_match <-
function(I)
{
    I <- .N.(I)
    ind <- which(diag(I) > 0)
    diag(I) <- 0
    ind <- rbind(cbind(ind, ind),
                 which(.T.(I, t(I)) > 0, arr.ind = TRUE))
    ind[ind[, 1L] <= ind[, 2L], , drop = FALSE]
}

.amount_by_which_relation_is_not_reflexive <-
function(I)
    sum(1 - diag(I))

.tuples_for_which_relation_is_not_reflexive <-
function(I)
    matrix(which(diag(I) < 1), ncol = 1L)

.amount_by_which_relation_is_not_irreflexive <-
function(I)
    sum(diag(I))

.tuples_for_which_relation_is_not_irreflexive <-
function(I)
    matrix(which(diag(I) > 0), ncol = 1L)
    
.amount_by_which_relation_is_not_coreflexive <-
function(I)
    sum(I[row(I) != col(I)])

.tuples_for_which_relation_is_not_coreflexive <-
function(I)
{
    ind <- which(I > 0, arr.ind = TRUE)
    ind[ind[, 1L] != ind[, 2L], , drop = FALSE]
}

.amount_by_which_relation_is_not_symmetric <-
function(I)
    sum(abs(I - t(I))) / 2

.tuples_for_which_relation_is_not_symmetric <-
function(I)
{
    ind <- which(I != t(I), arr.ind = TRUE)
    ind[ind[, 1L] <= ind[, 2L], , drop = FALSE]
}

.amount_by_which_relation_is_not_asymmetric <-
function(I)
{
    D <- diag(I)
    diag(I) <- 0
    sum(.T.(I, t(I))) / 2 + sum(D)
}

.tuples_for_which_relation_is_not_asymmetric <-
function(I)
{
    ind <- which(.T.(I, t(I)) > 0, arr.ind = TRUE)
    ind[ind[, 1L] <= ind[, 2L], , drop = FALSE]
}

.amount_by_which_relation_is_not_antisymmetric <-
function(I)
{
    diag(I) <- 0
    sum(.T.(I, t(I))) / 2
}

.tuples_for_which_relation_is_not_antisymmetric <-
function(I)
{
    ind <- which(.T.(I, t(I)) > 0, arr.ind = TRUE)
    ind[ind[, 1L] < ind[, 2L], , drop = FALSE]
}    

.amount_by_which_relation_is_not_transitive <-
function(I)
    sum(sapply(seq_len(nrow(I)),
               function(j) pmax(outer(I[, j], I[j, ], .T.) - I, 0)))

.tuples_for_which_relation_is_not_transitive <-
function(I)
{
    ind <- lapply(seq_len(nrow(I)),
                  function(j) {
                      pos <- which(outer(I[, j], I[j, ], .T.) > I,
                                   arr.ind = TRUE)
                      cbind(pos, rep.int(j, nrow(pos)))
                  })
    do.call(rbind, ind)[, c(1L, 3L, 2L), drop = FALSE]
}

.amount_by_which_relation_is_not_negatively_transitive <-
function(I)
    sum(sapply(seq_len(nrow(I)),
               function(j) pmax(I - outer(I[, j], I[j, ], .S.), 0)))

.tuples_for_which_relation_is_not_negatively_transitive <-
function(I)
{
    ind <- lapply(seq_len(nrow(I)),
                  function(j) {
                      pos <- which(outer(I[, j], I[j, ], .S.) < I,
                                   arr.ind = TRUE)
                      cbind(pos, rep.int(j, nrow(pos)))
                  })
    do.call(rbind, ind)[, c(1L, 3L, 2L), drop = FALSE]
}

.amount_by_which_relation_is_not_Ferrers <-
function(I)
{
    out <- 0
    for(j in seq_len(nrow(I)))
        for(l in seq_len(nrow(I)))
            out <- out + sum(pmax(outer(I[, j], I[, l], .T.) -
                                  outer(I[, l], I[, j], .S.),
                                  0))
    out
}

.tuples_for_which_relation_is_not_Ferrers <-
function(I)
{
    n <- nrow(I)
    ind <- Map(function(j, l) {
                   pos <- which(outer(I[, j], I[, l], .T.) >
                                outer(I[, l], I[, j], .S.),
                                arr.ind = TRUE)
                   cbind(pos,
                         rep.int(j, nrow(pos)),
                         rep.int(l, nrow(pos)))
               },
               rep.int(seq_len(n), n),
               rep(seq_len(n), each = n))
    do.call(rbind, ind)[, c(1L, 3L, 2L, 4L), drop = FALSE]
}

.amount_by_which_relation_is_not_semitransitive <-
function(I)
{
    out <- 0
    for(k in seq_len(nrow(I)))
        for(l in seq_len(nrow(I)))
            out <- out + sum(pmax(outer(I[, k], I[k, ], .T.) -
                                  outer(I[, l], I[l, ], .S.),
                                  0))
    out
}

.tuples_for_which_relation_is_not_semitransitive <-
function(I)
{
    n <- nrow(I)
    ind <- Map(function(k, l) {
                   pos <- which(outer(I[, k], I[k, ], .T.) >
                                outer(I[, l], I[l, ], .S.),
                                arr.ind = TRUE)
                   cbind(pos,
                         rep.int(k, nrow(pos)),
                         rep.int(l, nrow(pos)))
               },
               rep.int(seq_len(n), n),
               rep(seq_len(n), each = n))
    do.call(rbind, ind)
}   

.amount_by_which_relation_is_not_trichotomous <-
function(I)
    sum(diag(I)) + sum(1 - abs(I - t(I))[row(I) != col(I)]) / 2

.tuples_for_which_relation_is_not_trichotomous <-
function(I)
{
    ind <- which(abs(I - t(I)) < 1, arr.ind = TRUE)
    ind <- ind[ind[, 1L] < ind[, 2L], , drop = FALSE]
    pos <- which(diag(I) > 0)
    rbind(ind, cbind(pos, pos))
}

.amount_by_which_relation_is_not_Euclidean <-
function(I)
    sum(sapply(seq_len(nrow(I)),
               function(i) pmax(outer(I[i, ], I[i, ], .T.) - I, 0)))

.tuples_for_which_relation_is_not_Euclidean <-
function(I)
{
    ind <- lapply(seq_len(nrow(I)),
                  function(i) {
                      pos <- which(outer(I[i, ], I[i, ], .T.) > I,
                                   arr.ind = TRUE)
                      cbind(rep.int(i, nrow(pos)), pos)
                  })
    do.call(rbind, ind)
}
