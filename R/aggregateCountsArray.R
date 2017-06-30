################################################################################
### Aggregate an array in one dimension
###
### Copyright (C) 2014-2017 Sebastian Meyer
###
### This file is part of the R package "hhh4contacts",
### free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at http://www.r-project.org/Licenses/.
################################################################################

##' Aggregate an Array of Counts wrt One Dimension (Stratum)
##'
##' @param counts an (integer) array of counts with \code{dimnames},
##' e.g., \code{\link{counts}} or \code{\link{pop2011}}.
##' @param dim the dimension index of the stratum defining groups.
##' @param grouping,... how the groups should be built.
##' @param sort logical indicating if the resulting array should be ordered
##' by the grouping levels in the \code{dim} dimension.
##' @return an array with similar dimensions as the input \code{counts}
##' except for the \code{dim} dimension, which will be smaller due to the
##' aggregation as specified by the \code{grouping} argument.
##' @author Sebastian Meyer
##' @examples
##' ## works for matrices
##' aggregateCountsArray(pop2011, dim = 2, grouping = c(2,1,3,2,4))
##' aggregateCountsArray(pop2011, dim = 1, grouping = list(
##'     "a" = c("chwi","span","zehl"),
##'     "b" = c("neuk","scho")
##' ))
##' ## and of course for arrays
##' str(aggregateCountsArray(counts, dim = 3, grouping = c(1, 3, 4)))
##' @keywords manip
##' @export
aggregateCountsArray <- function (counts, dim, grouping, ..., sort = TRUE)
{
    stopifnot(is.array(counts), !is.null(odimnames <- dimnames(counts)),
              dim > 0, dim <= length(odim <- dim(counts)))
    if (missing(grouping)) {
        grouping <- list(...)
    }
    grouping <- grouping2list(grouping = grouping, groups = odimnames[[dim]])

    ## strata to keep
    strata2keep <- setdiff(odimnames[[dim]], unlist(grouping))
    counts_keep <- subset.array(counts, dim = dim, i = strata2keep)

    ## strata to aggregate
    counts_groups <- vapply(
        X = grouping,
        FUN = function (groups) {
            res <- if (dim == length(odim)) { # if stratum is last dimension
                rowSums(subset.array(counts, dim=dim, i=groups), dims = dim - 1L)
            } else { # have to use (slow) apply
                apply(subset.array(counts, dim=dim, i=groups), seq_along(odim)[-dim], sum)
            }
            storage.mode(res) <- storage.mode(counts) # keep "integer"
            res
        },
        FUN.VALUE = subset.array(counts, dim=dim, i=1, drop=TRUE),
        USE.NAMES = TRUE
    )
    if (dim != length(odim)) { ## fix order of dimensions
        counts_groups <- aperm(counts_groups,
                               perm = order(c(seq_along(odim)[-dim], dim)))
    }

    ## put them together
    res <- c(counts_keep, counts_groups)
    stopifnot(all.equal(sum(counts), sum(res)))
    dim(res) <- local({odim[dim] <- length(strata2keep) + length(grouping); odim})
    dimnames(res) <- local({odimnames[[dim]] <- c(strata2keep, names(grouping)); odimnames})

    ## sort
    if (sort) {
        subset.array(res, dim = dim, i = order(dimnames(res)[[dim]]))
    } else {
        res
    }
}


##' Subset an Array in one Dimension
##' @param x an array.
##' @param dim an integer specifying the dimension to subset.
##' @param i subset index, see \code{\link{[}}.
##' @param drop logical indicating if singular dimensions should be dropped,
##' see \code{\link{[}}.
##' @seealso the more general function \code{extract.array} in package
##' \pkg{R.utils}
##' @examples
##' hhh4contacts:::subset.array(counts, 1, 4:7)
##' @keywords internal
subset.array <- function (x, dim, i, drop = FALSE)
{
    ndim <- length(dim(x))
    stopifnot(dim > 0, dim <= ndim)
    args <- rep.int(alist(i=), ndim)
    args[[dim]] <- i
    do.call("[", c(alist(x), args, drop = drop))
}
