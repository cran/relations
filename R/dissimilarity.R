### Relation metrics (dissimilarities).

## We use the CLUE approach: have a front-end
##   relation_dissimilarity <-
##       function(x, y = NULL, method = "symdiff", ...)
## which returns dissimilarities or cross-dissimilarities, eventually
## with classes and subscript methods as in CLUE.

## If dissimilarity computations can be factored as f(tx(x), ty(y)) with
## non-trivial transformations tx() and/or ty() it is more time
## efficient to transform the elements of x and/or y once instead of
## every time a dissimilarity between elements of x and y is computed.
## Hence, for built-in dissimilarity methods we provide a mechanism for
## pre-transforming.
## (The same could be done for user-defined methods by providing and
## exposing a relation_dissimilarity_method() generator.)

relation_dissimilarity <-
function(x, y = NULL, method = "symdiff", ...)
{
    x <- as.relation_ensemble(x)

    ## Be nice.
    if(is.character(y) || is.function(y)) {
        method <- y
        y <- NULL
    }

    if(is.character(method)) {
        ## Hopefully of length one, add some tests eventually ...
        entry <- get_relation_dissimilarity_method(method)
        ## Check additional arguments if necessary.
        if(length(nms <- names(list(...)))) {
            ind <- pmatch(nms, entry$additional)
            if(any(is.na(ind))) {
                stop(gettextf("Invalid additional argument '%s'.",
                              nms[which(ind)[1L]]),
                     domain = NA)
            }
        }
        description <- entry$description
        method <- entry$method
        xtrafo <- entry$xtrafo
        ytrafo <- entry$ytrafo
    }
    else if(is.function(method)) {
        description <- "user-defined method"
        xtrafo <- ytrafo <- NULL
    }
    else
        stop("Invalid 'method' argument.")

    ## Canonicalize method and trafos, expanding generators.
    ## Note that generators are currently called with *all* args in
    ## ... and hence should accept all possible additional args.

    fmethod <- if(inherits(method, "Generator"))
        method(...)
    else if(length(formals(method)) > 2L)
        function(x, y) method(x, y, ...)
    else
        method

    fxtrafo <- if(is.null(xtrafo))
        NULL
    else if(inherits(xtrafo, "Generator"))
        xtrafo(...)
    else if(length(formals(xtrafo)) > 1L)
        function(x) xtrafo(x, ...)
    else
        xtrafo

    fytrafo <- if(is.null(ytrafo))
        NULL
    else if(inherits(ytrafo, "Generator"))
        ytrafo(...)
    else if(length(formals(ytrafo)) > 1L)
        function(y) ytrafo(y, ...)
    else
        ytrafo

    ## Build a cross-proximity object of cross-dissimilarities if y is
    ## given, and a proximity object of dissimilarities otherwise.
    ## <NOTE>
    ## If we do not require clue (or provide a similar infrastructure
    ## for (cross-)proximity objects, all we can return is a matrix or a
    ## dist object.
    ## </NOTE>

    if(!is.null(y)) {
        y <- as.relation_ensemble(y)
        D <- relation_domain(x)
        if(!.domain_is_equal(relation_domain(y), D))
            stop("All relations must have the same domain.")
        if(!is.null(fxtrafo)) x <- lapply(x, fxtrafo)
        if(!is.null(fytrafo)) y <- lapply(y, fytrafo)
        d <- matrix(0, length(x), length(y))
        for(j in seq_along(y))
            d[, j] <- sapply(x, fmethod, y[[j]])
        dimnames(d) <- list(names(x), names(y))
        attr(x, "description") <- description
        return(d)
    }

    y <- x
    if(!is.null(fxtrafo)) x <- lapply(x, fxtrafo)
    if(!is.null(fytrafo)) y <- lapply(y, fytrafo)
    n <- length(x)
    d <- vector("list", length = n - 1L)
    ind <- seq_len(n)
    while(length(ind) > 1L) {
        j <- ind[1L]
        ind <- ind[-1L]
        d[[j]] <- sapply(x[ind], fmethod, y[[j]])
    }
    ## Grr ... no generator for dist objects.
    ## See clue:::.dist_from_vector().
    .structure(unlist(d), Size = n, Labels = names(x),
               class = "dist", description = description)
}

## * relation_dissimilarity infrastructure

relation_dissimilarity_methods_db <- new.env()

set_relation_dissimilarity_method <-
function(name,
         description, method, xtrafo = NULL, ytrafo = NULL,
         additional = character(), ...)
{
    entry <- c(list(method = method,
                    description = description,
                    xtrafo = xtrafo,
                    ytrafo = ytrafo,
                    additional = additional),
               list(...))
    ## This could check that generators take all additional arguments.
    class(entry) <- "relation_dissimilarity_method"
    relation_dissimilarity_methods_db[[name]] <- entry
}

get_relation_dissimilarity_method <-
function(name)
{
    keys <- names(relation_dissimilarity_methods_db)
    if(is.na(ind <- pmatch(name, keys)))
        stop(gettextf("Method '%s' is not a valid dissimilarity method.",
                      name),
             domain = NA)
    relation_dissimilarity_methods_db[[keys[ind]]]
}

Generator <-
function(f)
{
    class(f) <- "Generator"
    f
}

## * relation_dissimilarity methods

## Note that the actual registration happens in zzz.R to make sure that
## everything needed is already available.

## ** relation_dissimilarity method 'symdiff'

.relation_dissimilarity_symdiff_params <- c("na.rm")

.relation_dissimilarity_symdiff_xtrafo <-
function(x)
{
    if(!identical(relation_is_crisp(x), TRUE))
        stop("Not implemented.")        # David said so ...
    relation_incidence(x)
}

.relation_dissimilarity_symdiff_method <-
Generator(function(na.rm = FALSE) {
              function(x, y)
                  sum(abs(x - y), na.rm = na.rm)
          })

## set_relation_dissimilarity_method("symdiff",
##                                   "symmetric difference distance",
##                                   .relation_dissimilarity_symdiff_method,
##                                   .relation_dissimilarity_symdiff_xtrafo,
##                                   .relation_dissimilarity_symdiff_xtrafo,
##                                   .relation_dissimilarity_symdiff_params)

## Old style without trafos:

## .relation_dissimilarity_symdiff <-
## function(x, y, na.rm = FALSE)
## {
##     if(!identical(relation_is_crisp(x), TRUE) ||
##        !identical(relation_is_crisp(y), TRUE))
##         stop("Not implemented.")        # David said so ...
##     .incidence_dissimilarity_symdiff(relation_incidence(x),
##                                      relation_incidence(y),
##                                      na.rm = na.rm)
## }

## (Still used in consensus code.)
.incidence_dissimilarity_symdiff <-
function(x, y, na.rm = FALSE)
    sum(abs(x - y), na.rm = na.rm)

## ** relation_dissimilarity method 'SD'

## set_relation_dissimilarity_method("SD",
##                                   "symmetric difference distance",
##                                   .relation_dissimilarity_symdiff_method,
##                                   .relation_dissimilarity_symdiff_xtrafo,
##                                   .relation_dissimilarity_symdiff_xtrafo,
##                                   .relation_dissimilarity_symdiff_params)

## ** relation_dissimilarity method 'CKS'

## Wade D. Cook and Moshe Kress and Lawrence M. Seiford
## Information and preference in partial orders: a bimatrix
## representation.
## Psychometrika 51/2. 197-207.
## Unique paired comparison metric for partial rankings (Definition 2.1)
## under the assumption that indifference is the centroid between strict
## preferences and incomparability.
## <NOTE>
## Originally only defined between (partial) rankings, but applicable
## more generally.
## </NOTE>

.relation_dissimilarity_CKS_xtrafo <-
function(x)
{
    x <- relation_incidence(x)
    list(P = pmin(t(x), 1 - x),
         I = pmax(x, t(x)))
}

.relation_dissimilarity_CKS_method <-
function(x, y)
{
    sum(abs(x$P - y$P)) + sum(abs(x$I - y$I)) / 2
}

## set_relation_dissimilarity_method("CKS",
##                                   "Cook-Kress-Seiford distance",
##                                   .relation_dissimilarity_CKS_method,
##                                   .relation_dissimilarity_CKS_xtrafo,
##                                   .relation_dissimilarity_CKS_xtrafo)

## Old style without trafos:

## .relation_dissimilarity_CKS <-
## function(x, y)
## {
##     .incidence_dissimilarity_CKS(relation_incidence(x),
##                                  relation_incidence(y))
## }

## (Still used in consensus code.)
.incidence_dissimilarity_CKS <-
function(x, y)
{
    P_x <- pmin(t(x), 1 - x)
    I_x <- pmax(x, t(x))
    P_y <- pmin(t(y), 1 - y)
    I_y <- pmax(y, t(y))
    sum(abs(P_x - P_y)) + sum(abs(I_x - I_y)) / 2
}

## ** relation_dissimilarity method 'CS'

## Wade D. Cook and Lawrence M. Seiford
## Priority Ranking and Consensus Formation
## Management Science 24/16, 1721--1732
## Ordinal ranking in the sense of Kendall CHECK
## Equivalent to: complete and transitive?
## I.e., weak order aka preference ...
## <NOTE>
## Originally only defined between (complete) rankings, but applicable
## more generally.
## </NOTE>

.relation_dissimilarity_CS_xtrafo <-
function(x)
    .incidence_scores_ranks(relation_incidence(x))

.relation_dissimilarity_CS_method <-
function(x, y)
    sum(abs(x - y))

## set_relation_dissimilarity_method("CS",
##                                   "Cook-Seiford distance",
##                                   .relation_dissimilarity_CS_method,
##                                   .relation_dissimilarity_CS_xtrafo,
##                                   .relation_dissimilarity_CS_xtrafo)

## Old style without trafos:

## .relation_dissimilarity_CS <-
## function(x, y)
## {
##     .incidence_dissimilarity_CS(relation_incidence(x),
##                                 relation_incidence(y))
## }

## (Still used in consensus code.)
.incidence_dissimilarity_CS <-
function(x, y)
    sum(abs(.incidence_scores_ranks(x) - .incidence_scores_ranks(y)))

## ** relation_dissimilarity method 'score'

## Score-based distance Delta(score(x), score(y)).
## If score is NULL, the default relation_score() is used.
## If it is a character string, relation_scores(x, score) is used.
## Otherwise, it must be a function (the score function itself).
## If Delta is a number p (\ge 1), the p-norm is used as Delta().
## Otherwise, it must be a fucntion.
## The generalized Cook-Seiford dissimilarity is a special case of a
## score-based distance (corresponding to the defaults).

.relation_dissimilarity_score_params <- c("score", "Delta")

.relation_dissimilarity_score_xtrafo <-
Generator(function(score = NULL, Delta = 1) {
              if(is.null(score))
                  relation_scores
              else if(is.character(score) && (length(score) == 1L))
                  function(x) relation_scores(x, score)
              else if(is.function(score))
                  score
              else
                  stop("Invalid 'score' argument.")
          })

.relation_dissimilarity_score_method <-
Generator(function(score = NULL, Delta = 1) {
              if(is.numeric(p <- Delta) &&
                 (length(p) == 1L) && (p >= 1))
                  function(x, y) sum(abs(x - y) ^ p) ^ (1 / p)
              else if(is.function(Delta))
                  Delta
              else
                  stop("Invalid 'Delta' argument.")
          })

## set_relation_dissimilarity_method("score",
##                                   "score-based distance",
##                                   .relation_dissimilarity_score_method,
##                                   .relation_dissimilarity_score_xtrafo,
##                                   .relation_dissimilarity_score_xtrafo,
##                                   .relation_dissimilarity_score_params)

## Old style without trafos:

## .relation_dissimilarity_score <
## function(x, y, score = NULL, Delta = 1)
## {
##     if(is.null(score)) {
##         s_x <- relation_scores(x)
##         s_y <- relation_scores(y)
##     } else if(is.character(score) && (length(score) == 1L)) {
##         s_x <- relation_scores(x, score)
##         s_y <- relation_scores(y, score)
##     } else if(is.function(score)) {
##         s_x <- score(x)
##         s_y <- score(y)
##     }
##     else
##         stop("Invalid 'score' argument.")
##
##     if(is.numeric(p <- Delta) && (length(p) == 1L) && (p >= 1)) {
##         Delta <- function(u, v) sum(abs(u - v) ^ p) ^ (1 / p)
##     } else if(!is.function(Delta))
##         stop("Invalid 'Delta' argument.")
##
##     Delta(s_x, s_y)
## }

## ** relation_dissimilarity method 'manhattan'

.relation_dissimilarity_manhattan_params <- c("na.rm")

.relation_dissimilarity_manhattan_method <-
Generator(function(na.rm = FALSE) {
              function(x, y)
                  sum(abs(x - y), na.rm = na.rm)
          })

## set_relation_dissimilarity_method("manhattan",
##                                   "Manhattan distance",
##                                   .relation_dissimilarity_manhattan_method,
##                                   relation_incidence,
##                                   relation_incidence,
##                                   .relation_dissimilarity_manhattan_params)

## Old style without trafos:

## .relation_dissimilarity_manhattan <-
## function(x, y, na.rm = FALSE)
##     .incidence_dissimilarity_manhattan(relation_incidence(x),
##                                        relation_incidence(y),
##                                        na.rm = na.rm)

## (Still used in consensus code.)
.incidence_dissimilarity_manhattan <-
function(x, y, na.rm = FALSE)
    sum(abs(x - y), na.rm = na.rm)

## ** relation_dissimilarity method 'euclidean'

.relation_dissimilarity_euclidean_params <- c("na.rm")

.relation_dissimilarity_euclidean_method <-
Generator(function(na.rm = FALSE) {
              function(x, y)
                  sqrt(sum((x - y) ^ 2, na.rm = na.rm))
          })

## set_relation_dissimilarity_method("euclidean",
##                                   "Euclidean distance",
##                                   .relation_dissimilarity_euclidean_method,
##                                   relation_incidence,
##                                   relation_incidence,
##                                   .relation_dissimilarity_euclidean_params)

## Old style without trafos:

## .relation_dissimilarity_euclidean <-
## function(x, y, na.rm = FALSE)
##     .incidence_dissimilarity_euclidean(relation_incidence(x),
##                                        relation_incidence(y),
##                                        na.rm = na.rm)

## (Still used in consensus code.)
.incidence_dissimilarity_euclidean <-
function(x, y, na.rm = FALSE)
    sqrt(sum((x - y) ^ 2, na.rm = na.rm))

## ** relation_dissimilarity method 'Jaccard'

.relation_dissimilarity_Jaccard_params <- c("na.rm")

.relation_dissimilarity_Jaccard_method <-
Generator(function(na.rm = FALSE) {
              function(x, y) {
                  if(identical(all(x == 0, na.rm = na.rm), TRUE) &&
                     identical(all(y == 0, na.rm = na.rm), TRUE))
                      return(0)
                  1 - sum(.T.(x, y), na.rm = na.rm) /
                      sum(.S.(x, y), na.rm = na.rm)
              }
          })


## set_relation_dissimilarity_method("Jaccard",
##                                   "Jaccard distance",
##                                   .relation_dissimilarity_Jaccard_method,
##                                   relation_incidence,
##                                   relation_incidence,
##                                   .relation_dissimilarity_Jaccard_params)

## Old style without trafos:

## ## One could also use
## ##   1 - gset_similarity(relation_graph(x), relation_graph(y))
## ## but proceeding directly should be more efficient.
##
## .relation_dissimilarity_Jaccard <-
## function(x, y, na.rm = FALSE)
##     .incidence_dissimilarity_Jaccard(relation_incidence(x),
##                                      relation_incidence(y),
##                                      na.rm = na.rm)
##
## .incidence_dissimilarity_Jaccard <-
## function(x, y, na.rm = FALSE)
## {
##     if(identical(all(x == 0, na.rm = na.rm), TRUE) &&
##        identical(all(y == 0, na.rm = na.rm), TRUE))
##         return(0)
##     1 - sum(.T.(x, y), na.rm = na.rm) /
##         sum(.S.(x, y), na.rm = na.rm)
## }

## ** relation_dissimilarity method 'PC'

.relation_dissimilarity_PC_params <-
    c("delta", "gamma", "family")

.relation_dissimilarity_PC_xtrafo <-
Generator(function(delta = "symdiff", gamma = NULL, family = NULL) {
              D <- .relation_dissimilarity_PC_Delta(delta)
              M <- .relation_dissimilarity_PC_M(D)
              function(x) {
                  IP <- .relation_dissimilarity_PC_IP(x, family)
                  .relation_dissimilarity_PC_ABC(IP, delta, gamma,
                                                 D, M)
              }
          })

.relation_dissimilarity_PC_ytrafo <-
Generator(function(delta = "symdiff", gamma = NULL, family = NULL) {
              function(y)
                  .relation_dissimilarity_PC_IP(y, family)
          })

.relation_dissimilarity_PC_method <-
function(x, y) {
    sum(x$C) + sum(x$B * y$I) + sum(x$A * y$P)
}

## set_relation_dissimilarity_method("PC",
##                                   "paired comparison distance",
##                                   .relation_dissimilarity_PC_method,
##                                   .relation_dissimilarity_PC_xtrafo,
##                                   .relation_dissimilarity_PC_ytrafo,
##                                   .relation_dissimilarity_PC_params)

## Old style without trafos:

## .relation_dissimilarity_PC <-
## function(x, y, delta = "symdiff", gamma = NULL)
##     .incidence_dissimilarity_PC(relation_incidence(x),
##                                 relation_incidence(y),
##                                 delta,
##                                 gamma)

## (Still used in consensus code.)
.incidence_dissimilarity_PC <-
function(IP, Y, delta, gamma)
{
    ABC <- .relation_dissimilarity_PC_ABC(IP, delta, gamma)
    sum(ABC$C) + sum(ABC$B * Y) + sum(ABC$A * Y * t(Y))
}

.relation_dissimilarity_PC_ABC <-
function(IP, delta, gamma, D = NULL, M = NULL)
{
    X <- IP$I
    P <- IP$P

    if(is.null(D))
        D <- .relation_dissimilarity_PC_Delta(delta)
    if(is.null(M))
        M <- .relation_dissimilarity_PC_M(D)

    T <- t(X)
    if(is.null(P))
        P <- X * T

    nrx <- nrow(X)
    dnx <- dimnames(X)

    x <- diag(X)

    M1 <-             M[1L, 2L] * X + M[1L, 3L] * T + M[1L, 4L] * P
    M2 <- M[2L, 1L] + M[2L, 2L] * X + M[2L, 3L] * T + M[2L, 4L] * P
    ## Actually compute t(M3): more efficient to transpose here than
    ## when subassigning below.
    M3 <- M[3L, 1L] + M[3L, 2L] * T + M[3L, 3L] * X + M[3L, 4L] * P
    M4 <- M[4L, 1L] + M[4L, 2L] * X + M[4L, 3L] * T + M[4L, 4L] * P

    upper <- (row(X) < col(X))
    lower <- (row(X) > col(X))

    d <- D[4L, 1L]

    if(is.null(gamma)) {
        ## Avoid unnecessary multiplications by \gamma_{ij} = 1.
        C <- diag(d * x / 2, nrx, nrx)
        C[upper] <- M1[upper]
        B <- diag(d * (1/2 - x), nrx, nrx)
        B[upper] <- M2[upper]
        B[lower] <- M3[lower]
        A <- matrix(0, nrx, nrx)
        A[upper] <- M4[upper]
    } else {
        if(!identical(dim(gamma), c(nrx, nrx)))
            stop("Invalid PC pair weights.")
        gii <- diag(gamma)
        gij <- gamma[upper]
        gji <- t(gamma)[lower]
        C <- diag(gii * d * x, nrx, nrx)
        C[upper] <- gij * M1[upper]
        B <- diag(gii * d * (1 - 2 * x), nrx, nrx)
        B[upper] <- gij * M2[upper]
        B[lower] <- gji * M3[lower]
        A <- matrix(0, nrx, nrx)
        A[upper] <- M4[upper]
    }

    dimnames(C) <- dnx
    dimnames(B) <- dnx
    dimnames(A) <- dnx

    list(C = C, B = B, A = A)
}

.relation_dissimilarity_PC_T <-
    matrix(c( 1,  0,  0, 0,
             -1,  1,  0, 0,
             -1,  0,  1, 0,
              1, -1, -1, 1),
           4L, 4L)

.relation_dissimilarity_PC_M <-
function(D)
{
    T <- .relation_dissimilarity_PC_T
    if(is.character(D))
        D <- .relation_dissimilarity_PC_Delta(D)
    crossprod(T, D %*% T)
}

.relation_dissimilarity_PC_Delta <-
function(x)
{
    if(identical(dim(x), c(4L, 4L)))
        x
    else {
        x <- .relation_dissimilarity_PC_delta(x)
        D <- matrix(0, 4L, 4L)
        i <- c(2L, 3L, 4L, 3L, 4L, 4L)
        j <- c(1L, 1L, 1L, 2L, 2L, 3L)
        D[cbind(c(i, j), c(j, i))] <- c(x, x)
        D
    }
}

.relation_dissimilarity_PC_delta <-
function(x)
{
    known <-
        list(symdiff =
                 c(1, 1, 2, 2, 1, 1),
             SD =
                 c(1, 1, 2, 2, 1, 1),
             CKS =
                 c(2, 2, 1, 2, 1, 1),
             KS =
                 c(1, 1, 0, 2, 1, 1),
             EM =
                 c(1, 1, 1, 2, 1, 1),
             JMB =
                 c(4/3, 4/3, 4/3, 5/3, 1, 1),
             discrete =
                 rep.int(1, 6L)
             )
    if(is.character(x) && !is.na(pos <- pmatch(x, names(known))))
        known[[pos]]
    else if(is.numeric(x) && (length(x) == 6L))
        x
    else
        stop("Invalid PC discrepancies.")
}

## averaged incidences for W with PC dist
.impute_P <-
function(I, n)
{
    o <- do.call(order, split(I, col(I)))
    ro <- `[<-`(seq_along(o), o, seq_along(o))
    I <- I[o, o, drop = FALSE]

    mn <- ncol(I)
    m <- mn - n

    A <- I[1:m, 1:m]
    c  <- length(table(cumsum(!duplicated(A))))

    A <- A * t(A)

    f <- .nsol_W(c, n - 1) / .nsol_W(c, n)
    B <- matrix(f, m, n)
    D <- matrix(f, n, n)
    diag(D) <- 1

    rbind(cbind(A, B), cbind(t(B), D))[ro, ro, drop = FALSE]
}


.relation_dissimilarity_PC_IP <-
function(x, family = NULL)
{
    X <- relation_incidence(x)
    ## If X has no missings, return a list with X and P = X * t(X) [the
    ## incidences of the symmetric part of x].
    ## Otherwise, if family is suitable, return a list of the average X
    ## and P over all imputations of x within the family (if possible).

    if(is.character(family) && anyNA(X) && (length(family) == 1L)) {
        if(family == "L") {
            ## FIXME: make testing for possible imputations in L more
            ## efficient.
            if(!relation_is_linear_order(x, na.rm = TRUE))
                stop("Cannot impute within family 'L'.")
            X <- .incidence(relation_impute(x, "average/L"))
            P <- diag(1, nrow(X), ncol(X))
        } else if(family == "W") {
            ## FIXME: make testing for possible imputations in W more
            ## efficient.
            if(!relation_is_weak_order(x, na.rm = TRUE))
                stop("Cannot impute within family 'W'.")
            P <- .impute_P(X, length(.missing_objects(X)))
            X <- .incidence(relation_impute(x, "average/W"))
        } else
            stop(gettextf("Invalid family '%s'.", family),
                 domain = NA)
    } else {
        P <- X * t(X)
    }

    list(I = X, P = P)
}
