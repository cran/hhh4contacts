################################################################################
### Compute the Dawid-Sebastiani Score on aggregated predictions
###
### Copyright (C) 2015 Sebastian Meyer
###
### This file is part of the R package "hhh4contacts",
### free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at http://www.r-project.org/Licenses/.
################################################################################

##' Compute the DSS on Aggregated Predictions and Observations
##'
##' The expectation and variance of aggregated predictions is just a sum if the
##' predictions are (conditionally) independent. This function computes the DSS
##' for a matrix of observations and a matrix of predictions
##' where the columns are to be summed according to a given factor.
##' @param observed a numeric matrix of observed counts.
##' @param pred a numeric matrix of predicted counts.
##' @param psi a numeric vector or matrix of overdispersion parameters such that
##' \code{pred * (1 + pred/exp(psi))} is the prediction's variance.
##' Alternatively, \code{psi = NULL} indicated Poisson predictions.
##' @param groups a factor variable of length \code{ncol(observed)} indicating
##' which columns should be aggregated.
##' @return a matrix of DSS values
##' @export

dssAggregate <- function (observed, pred, psi, groups)
{
    x <- rowSumsBy(observed, groups)
    EV <- aggregatePredictions(pred, psi, groups)
    DSS(meanP = EV[["mean"]], varP = EV[["var"]], x = x)
}

aggregatePredictions <- function (pred, psi, groups)
{
    ## mean and variance by group
    gmean <- rowSumsBy(pred, groups)
    gvar <- if (is.null(psi)) {
        gmean
    } else {
        ## psi is -log(overdispersion)
        size <- exp(c(psi))
        ## size <- rep_len(size, length(pred))
        ## attributes(size) <- attributes(pred)
        rowSumsBy(pred * (1 + pred/size), groups)
    }
    list(mean = gmean, var = gvar)
}

DSS <- function (meanP, varP, x)
{
    (x-meanP)^2 / varP + log(varP)
}
