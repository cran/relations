### Relation scores.

relation_scores <-
function(x,
         method = c("ranks", "Barthelemy/Monjardet", "Borda",
                    "Kendall", "Wei", "differential"),
         normalize = FALSE, ...)
{
    method <- match.arg(method)
    if(!relation_is_endorelation(x))
        stop("Relation scores are only available for endorelations.")
    labs <- LABELS(relation_domain(x)[[1L]])

    ## <NOTE>
    ## When adding .relation_score_FOO(x, ...) method functions, we
    ## might want to pass x as the relation itself rather than its
    ## incidence.
    x <- relation_incidence(x)
    ## </NOTE>

    ret <- switch(method,
                  ranks = .relation_scores_ranks(x, ...),
                  "Barthelemy/Monjardet" = {
                      ## Use formula in Barthelemy & Monjardet, p. 258.
                      ## See also
                      ## http://mathworld.wolfram.com/ScoreSequence.html.
                      (colSums(x * (1 - t(x))) + colSums(x) - 1) / 2
                  },
                  Borda =, Kendall = colSums(x),
                  differential = colSums(x) - rowSums(x),
                  Wei = {
                      ## <FIXME>
                      ## Cook & Kress use "preference matrices", so to
                      ## be consistent, take complement?
                      abs(Re(eigen(1 - x)$vectors[, 1L]))
                      ## </FIXME>
                  })
    
    names(ret) <- labs
    if(normalize)
        ret / sum(ret)
    else
        ret
}

.relation_scores_ranks <-
function(x, decreasing = TRUE)
{
    n <- ncol(x)
    if(decreasing)
        (n + 1 + rowSums(x) - colSums(x)) / 2
    else
        (n + 1 + colSums(x) - rowSums(x)) / 2
}
    
