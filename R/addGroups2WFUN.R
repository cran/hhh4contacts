################################################################################
### Modify a parametric weight function to model group-dependent parameters
###
### Copyright (C) 2015-2017 Sebastian Meyer
###
### This file is part of the R package "hhh4contacts",
### free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at http://www.r-project.org/Licenses/.
################################################################################

##' Group-Dependent Parametric Weights
##'
##' This function takes a specification of parametric weights and returns
##' a modified version with group-dependent parameters.
##' Only single-parameter functions are currently supported.
##'
##' @param WFUN a list specification of parametric weights, e.g.,
##' as returned by the constructor functions \code{\link{W_powerlaw}} and
##' \code{\link{W_np}}.
##' @param groups a vector of length \code{nUnits} determining to which
##' group each unit belongs to. The supplied vector is converted to a
##' factor using \code{\link{as.factor}}.
##' @param initial (named) vector of initial parameters.
##' @return a list specifying group-dependent parametric weights for \code{hhh4}.
##' @author Sebastian Meyer
##' @importFrom surveillance clapply
##' @export
##' @examples
##' data("measlesWeserEms")
##' WPLgroups <- addGroups2WFUN(
##'   W_powerlaw(maxlag = 5, normalize = FALSE, log = FALSE),
##'   groups = factor(sample(2, ncol(measlesWeserEms), replace = TRUE)))

addGroups2WFUN <- function (WFUN, groups,
                            initial = rep.int(WFUN$initial, nlevels(groups)))
{
    ## parameters of the basic weight function
    npar <- length(WFUN$initial)
    stopifnot(npar == 1)
    parnames <- names(WFUN$initial)
    if (is.null(parnames)) parnames <- seq_len(npar)
    ## number of groups of units
    groups <- as.factor(groups)
    ngroups <- nlevels(groups)
    if (ngroups == 1L) {
        message("NOTE: 'groups' has only one level ...")
    }
    stopifnot(length(initial) == npar * ngroups)

    ## convert "neweights" parameter vector to a list (by group)
    coeflist <- function (pars_groups) {
        split.default(pars_groups,
                      gl(n = ngroups, k = npar, labels = levels(groups)))
    }

    ## evaluate weights for each parameter group (on the complete nbmat)
    evalWFUNbyGroup <- function (pars_groups, nbmat, data,
                                 which = c("w", "dw", "d2w")) {
        lapply(X = coeflist(pars_groups),
               FUN = function (pars) WFUN[[which]](pars, nbmat, data))
    }

    ## set rows to 0 if that units group is not the group of the list element
    ## returns a list of length ngroups of matrices (for "w") or of lists of
    ## length(WFUN$dw()), and length(WFUN$d2w()), respectively
    matchGroup <- function (wlist, groups) {
        mapply(
            FUN = function (weights, group) {
                ## weights: single matrix for "w", but a list for "dw"/"d2w"
                setTo0 <- function (W) {
                    W[groups != group,] <- 0  # for units not in this group
                    W
                }
                clapply(X = weights, FUN = setTo0)
            },
            weights = wlist, group = levels(groups),
            SIMPLIFY = FALSE, USE.NAMES = TRUE
            )
    }

    w <- function (pars_groups, nbmat, data) {
        wlist <- evalWFUNbyGroup(pars_groups, nbmat, data, which = "w")
        wlist0 <- matchGroup(wlist, groups)
        warray0 <- simplify2array(unlist(wlist0, recursive = FALSE, use.names = TRUE),
                                  higher = TRUE)
        res <- rowSums(warray0, dims = 2L)
        ## alternative: for each row j pick the results according to j's group
        stopifnot(all.equal(res,
            t(mapply(function (group, j) wlist[[group]][j,],
                     group = groups, j = seq_along(groups),
                     SIMPLIFY = TRUE, USE.NAMES = TRUE)),
                            check.attributes = FALSE))
        res
    }

    dw <- function (pars_groups, nbmat, data) {
        dwlist <- evalWFUNbyGroup(pars_groups, nbmat, data, which = "dw")
        dwlist0 <- matchGroup(dwlist, groups)
        unlist(dwlist0, recursive = FALSE, use.names = TRUE)
    }

    d2w <- function (pars_groups, nbmat, data) {
        d2wlist <- evalWFUNbyGroup(pars_groups, nbmat, data, which = "d2w")
        d2wlist0 <- matchGroup(d2wlist, groups)
        ## d2w has block diagonal structure
        isDiagEntry <- diag(ngroups) == 1
        isDiagEntry <- isDiagEntry[lower.tri(isDiagEntry, diag = TRUE)]
        ##<- upper triangle by row = lower triangle by column
        zeroEntry <- rep.int(list(matrix(0, nrow(nbmat), ncol(nbmat))), npar)
        res <- rep.int(list(zeroEntry), length(isDiagEntry))
        ## fill diagonal entries of upper triangle list
        res[isDiagEntry] <- d2wlist0
        unlist(res, recursive = FALSE, use.names = TRUE)
    }

    if (is.null(names(initial))) {
        names(initial) <- if (npar == 1) {
            levels(groups)
        } else {
            paste(rep.int(parnames, ngroups),
                  rep(levels(groups), each = npar),
                  sep = ".")
        }
    }
    list(w = w, dw = dw, d2w = d2w, initial = initial)
}
