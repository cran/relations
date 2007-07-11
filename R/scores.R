### Relation scores.

relation_scores <-
function(x,
         method = c("Barthelemy/Monjardet", "Wei", "Borda",
                    "Kendall", "differential", "rankB"),
         normalize = FALSE)
{
    method <- match.arg(method)
    if(!relation_is_endorelation(x))
        stop("Relation scores are only available for endorelations.")
    labs <- LABELS(relation_domain(x)[[1]])
    x <- relation_incidence(x)

    ret <- switch(method,
                  ## Use formula in Barthelemy & Monjardet, p. 258.
                  ## See also
                  ## http://mathworld.wolfram.com/ScoreSequence.html.
                  "Barthelemy/Monjardet" =
                  (colSums(x * (1 - t(x))) + colSums(x) - 1) / 2,
                  Borda =, Kendall = colSums(x),
                  differential = colSums(x) - rowSums(x),
                  rankB = (colSums(x) - rowSums(x) + ncol(x) + 1) / 2,
                  ## FIXME: Cook & Kress use "preference matrices", so to
                  ## be consistent, take complement?
                  Wei = abs(Re(eigen(1 - x)$vectors[,1]))
                  )
    names(ret) <- labs
    if(normalize)
        ret / sum(ret)
    else
        ret
}
