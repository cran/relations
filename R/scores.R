### Relation scores.

relation_scores <-
function(x)
{
    if(!relation_is_endorelation(x)) return(NA)
    ## Or whatever ...

    ## Use formula in Barthelemy & Monjardet, p. 258
    ## See also http://mathworld.wolfram.com/ScoreSequence.html.
    x <- relation_incidence(x)
    (colSums(x * (1 - t(x))) + colSums(x) - 1) / 2
}
