################################################################################
### Calculate rowSums for stratified columns (and rows) of a matrix (an array)
###
### Copyright (C) 2015 Sebastian Meyer
###
### This file is part of the R package "hhh4contacts",
### free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at http://www.r-project.org/Licenses/.
################################################################################

rowSumsBy <- function (x, by, ...) UseMethod("rowSumsBy")

rowSumsBy.matrix <- function (x, by, ...)
{
    by <- as.factor(by)
    stopifnot(length(by) == ncol(x))
    vapply(X = levels(by),
           FUN = function (g) rowSums(x[, by == g, drop = FALSE]),
           FUN.VALUE = numeric(nrow(x)), USE.NAMES = TRUE)
}

rowSumsBy.array <- function (x, by,
                             which3 = c("total", "own", "neighbours", "each"),
                             ...)
{
    by <- as.factor(by)
    which3 <- match.arg(which3)
    dimx <- dim(x)
    stopifnot(length(dimx) == 3L, length(by) == dimx[2L])
    if (which3 != "total")
        stopifnot(length(by) == dimx[3L])

    rowSums1 <- switch(which3,
        "total"      = function (g) rowSums(x[, by == g,        , drop = FALSE]),
        "own"        = function (g) rowSums(x[, by == g, by == g, drop = FALSE]),
        "neighbours" = function (g) rowSums(x[, by == g, by != g, drop = FALSE]),
        ## "each"       = function (g) apply(x[, by == g, , drop = FALSE], 3, rowSums)
        "each"       = function (g) vapply(X = levels(by),
                FUN = function (g2) rowSums(x[, by == g, by == g2, drop = FALSE]),
                FUN.VALUE = numeric(dimx[1L]), USE.NAMES = TRUE)
    )
    res <- vapply(X = levels(by),
                  FUN = rowSums1,
                  FUN.VALUE = if (which3 == "each") {
                      matrix(0, dimx[1L], nlevels(by))  # (t, j)
                  } else {
                      numeric(dimx[1L])
                  }, USE.NAMES = TRUE)
    if (which3 == "each") {
        ## dimension order of res is (t, j, i), change to (t, i, j) as 'x'
        names(dimnames(res)) <- names(dimnames(x))[c(1, 3, 2)]
        aperm(res, perm = c(1, 3, 2))
    } else {
        res
    }
}
