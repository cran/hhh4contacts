################################################################################
### Plot fitted values aggregated by a stratum variable, group-specific season
###
### Copyright (C) 2015-2017 Sebastian Meyer
###
### This file is part of the R package "hhh4contacts",
### free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at http://www.r-project.org/Licenses/.
################################################################################

##' Plot Mean Components of a \code{hhh4} Fit by Group
##'
##' Fitted mean components for age-structured, areal time series
##' \code{\link[surveillance]{hhh4}} models can be aggregated over districts or age groups.
##' @param x an object of class \code{"hhh4"}.
##' @param groups a factor of grouping the units in the model, i.e., it must be
##' of length \code{x$nUnit}. There will be one plot for each factor level.
##' @param total a logical indicating if the group-wise mean components should
##' be subsequently summed up over all \code{groups} for an overall plot.
##' @param decompose,... see \code{\link[surveillance]{plotHHH4_fitted}}.
##' @return see \code{\link[surveillance]{plotHHH4_fitted}}.
##' @importFrom surveillance plotHHH4_fitted
##' @export

plotHHH4_fitted_groups <- function (x, groups, total = FALSE,
                                    decompose = NULL, ...)
{
    groups <- as.factor(groups)
    if (isTRUE(decompose)) decompose <- levels(groups)

    ## aggregate the mean components within each group
    meanHHH_groups <- meanHHH_groups(x = x, groups = groups,
                                     decomposed = !is.null(decompose))

    if (total) # sum the mean components over all 'groups'
        meanHHH_groups <- if (is.null(decompose)) {
            lapply(X = meanHHH_groups, FUN = rowSumsBy.matrix,
                   by = rep.int("Overall", nlevels(groups)))
        } else {
            array(data = apply(X = meanHHH_groups, MARGIN = c(1L, 3L), FUN = sum),
                  dim = dim(meanHHH_groups)^c(1,0,1),
                  dimnames = "[[<-"(dimnames(meanHHH_groups), 2L, "Overall"))
        }

    ## hack 'x' for plotHHH4_fitted
    x$stsObj@observed <- rowSumsBy.matrix(
        x = x$stsObj@observed,
        by = if (total) rep.int("Overall", x$nUnit) else groups)
    x$nUnit <- ncol(x$stsObj@observed)
    x$control$ar$inModel <- x$control$ar$inModel || x$control$ne$inModel
    plotHHH4_fitted(x, meanHHH = meanHHH_groups, decompose = decompose, ...)
}

##' Plot Mean Components of a \code{hhh4} Fit by District Averaged Over Time
##'
##' This is a wrapper for \code{\link[surveillance]{plotHHH4_maps}} with prior aggregation
##' over different (age) groups.
##' @param x an object of class \code{"hhh4"}.
##' @param map an object inheriting from \code{"SpatialPolygons"}.
##' @param districts a factor of length \code{x$nUnit} with as many levels as
##'     there are districts and names according to \code{row.names(map)}.
##' @param ... arguments passed to \code{\link[surveillance]{plotHHH4_maps}}.
##' @return see \code{\link[surveillance]{plotHHH4_maps}}
##' @importFrom surveillance plotHHH4_maps
##' @export
plotHHH4_maps_groups <- function (x, map, districts, ...)
{
    meanHHH_districts <- meanHHH_groups(x, districts)
    plotHHH4_maps(x, ..., map = map, meanHHH = meanHHH_districts)
}

##' @importFrom surveillance decompose.hhh4 meanHHH
meanHHH_groups <- function (x, groups, decomposed = FALSE)
{
    meanHHH_decomposed <- decompose.hhh4(x)
    res <- if (decomposed) {
        endemic <- rowSumsBy.matrix(meanHHH_decomposed[,,1L], groups)
        epidemic <- rowSumsBy.array(meanHHH_decomposed[,,-1L,drop=FALSE], groups, which3 = "each")
        array(c(endemic, epidemic),
              dim = dim(epidemic) + c(0,0,1),
              dimnames = c(dimnames(epidemic)[1:2],
                           list("j" = c("endemic", dimnames(epidemic)[[3]]))))
    } else {
        list(mean = rowSumsBy.array(meanHHH_decomposed, groups, which3 = "total"),
             endemic = rowSumsBy.matrix(meanHHH_decomposed[,,1L], groups),
             epi.own = rowSumsBy.array(meanHHH_decomposed[,,-1L,drop=FALSE], groups, which3 = "own"),
             epi.neighbours = rowSumsBy.array(meanHHH_decomposed[,,-1L,drop=FALSE], groups, which3 = "neighbours"))
    }

    ## consistency checks (FIXME: to be removed in the future)
    if (decomposed) {
        stopifnot(
            all.equal(res[,,1L], drop(endemic), check.attributes = FALSE),
            identical(res[,,-1L], drop(epidemic)),
            all.equal(rowSumsBy.matrix(
                meanHHH(x$coefficients, surveillance:::terms.hhh4(x), total.only=TRUE), groups),
                rowSums(res, dims = 2), check.attributes = FALSE))
    } else {
        stopifnot(
            all.equal(res[["mean"]], with(res, endemic + epi.own + epi.neighbours)),
            all.equal(res[["mean"]], rowSumsBy.matrix(
                meanHHH(x$coefficients, surveillance:::terms.hhh4(x), total.only=TRUE), groups)))
    }

    res
}


##' Plot Seasonality of a \code{hhh4} Fit by Group
##'
##' A plot method for models with group-specific seasonality terms that are not
##' handled correctly by \code{\link[surveillance]{plotHHH4_season}}.
##' @param x an object of class \code{"hhh4"}.
##' @param component character string indicating from which component
##'     seasonality terms should be extracted.
##' @param seasonStart an integer defining the \code{\link[surveillance]{epochInYear}} that
##'     starts a new season (by default the first).
##' @param conf.level,conf.B a confidence level for the pointwise confidence
##'     intervals around the group-specific seasonal effects. The confidence
##'     intervals are based on quantiles of \code{conf.B} samples from the
##'     asymptotic multivariate normal distribution of the maximum likelihood
##'     estimate. Alternatively, if \code{conf.level = NA}, the individual
##'     samples are drawn instead of the confidence lines.
##'     Set \code{conf.level = NULL} to disable confidence intervals.
##' @param col a vector of group-specific colors, recycled as necessary and
##'     passed to \code{\link{matplot}}.
##' @param xlab,ylab,... arguments passed to \code{\link{matplot}}.
##' @param refline.args a list of arguments for \code{\link{abline}} to change
##'     the style of the horizontal reference line at 1.
##'     This line is omitted if \code{refline.args} is not a list.
##' @param yearline.args a list of arguments for \code{\link{abline}} to change
##'     the style of the line marking the end of the year at
##'     \code{x$stsObj@freq} if \code{seasonStart} is not 1.
##'     This line is omitted if \code{yearline.args} is not a list.
##' @param legend.args a list of arguments for \code{\link{legend}} modifying
##'     the internal defaults. If \code{legend.args} is not a list, the legend
##'     is omitted.
##' @return a matrix of the plotted point estimates of the multiplicative
##'     seasonal effect by group.
##' @importFrom graphics matplot matlines
##' @importFrom grDevices rgb col2rgb adjustcolor
##' @importFrom stats quantile
##' @importFrom utils modifyList
##' @export
plotHHH4_season_groups <- function (x, component = "end", seasonStart = 1,
    conf.level = 0.95, conf.B = 999,
    col = 1:6, xlab = "time", ylab = "multiplicative effect", ...,
    refline.args = list(), yearline.args = list(), legend.args = list())
{
    freq <- x$stsObj@freq
    seasonIndex <- c(seasonStart:freq, seq_len(seasonStart-1))
    seasonx <- seq.int(seasonStart, length.out = freq)
    seasonEstimate <- exp(
        getSeason_groups(x = x, component = component)[seasonIndex, , drop=FALSE]
    )

    colFit <- if (is.null(conf.level)) {
        ## establish an empty plot window
        matplot(x = seasonx, y = seasonEstimate, type = "n",
                xlab = xlab, ylab = ylab, ...)
        col
    } else {
        ## sample from the asymptotic normal of the MLE
        coefSamples <- MASS::mvrnorm(conf.B, mu = x$coefficients, Sigma = x$cov)
        logSeasonByGroup <- apply(X = rbind(x$coefficients, coefSamples),
                                  MARGIN = 1, FUN = getSeason_groups,
                                  x = x, component = component)
        dim(logSeasonByGroup) <- c(freq, nrow(logSeasonByGroup)/freq, 1 + conf.B)

        ## exp and cycle season
        season2plot <- exp(logSeasonByGroup[seasonIndex, , , drop = FALSE])

        ## plot simulation-based CI
        col <- rep_len(col, dim(season2plot)[2L])
        if (is.na(conf.level)) { # plot individual samples
            matplot(x = seasonx, y = matrix(season2plot, nrow = freq),
                    type = "l", lty = 1, col = adjustcolor(col, alpha.f = 0.1),
                    xlab = xlab, ylab = ylab, ...)
            darken <- function (col, f = 0.5)
                apply(col2rgb(col)/255*f, 2L, function (x) rgb(x[1L], x[2L], x[3L]))
            darken(col)
        } else { # plot conf.level confidence polygon
            alpha <- 1 - conf.level
            seasonCI <- apply(X = season2plot, MARGIN = 1:2, FUN = quantile,
                              probs = c(alpha/2, 1-alpha/2), names = FALSE)
            matplot(x = seasonx, y = matrix(aperm(seasonCI, perm = c(2,3,1)), nrow = freq),
                    type = "l", lty = 2, col = col, xlab = xlab, ylab = ylab, ...)
            col
        }
    }

    ## add reference line
    if (is.list(refline.args))
        do.call("abline", modifyList(list(h = 1, col = "gray"), refline.args))

    ## add year delimiter
    if (seasonStart != 1 && is.list(yearline.args))
        do.call("abline", modifyList(list(v = freq+.5, lty = 3), yearline.args))

    ## add point estimate
    matlines(x = seasonx, y = seasonEstimate,
             type = "l", lty = 1, lwd = 3, col = colFit, ...)

    ## add legend
    if (is.list(legend.args)) {
        legend.args <- modifyList(
            list(x = "topright", legend = colnames(seasonEstimate),
                 col = col, lwd = 3, bty = "n"),
            legend.args)
        do.call("legend", legend.args)
    }

    invisible(seasonEstimate)
}

getSeason_groups <- function (x, component = "end", coef = x$coefficients)
{
    startseason <- surveillance:::getSeasonStart(x)
    freq <- x$stsObj@freq
    t <- startseason - 1 + seq_len(freq)

    ## extract coefficients of sine and cosine terms
    coefSinCos <- coef[grep(paste0("^", component, "\\.(sin|cos)\\("), names(coef))]

    ## split group indicator from sin/cos term
    termsAndGroups <- strsplit(
        x = sub(paste0("^", component, "\\."), "", names(coefSinCos)),
        split = ".", fixed = TRUE)
    terms <- vapply(X = termsAndGroups, FUN = "[", 1L, FUN.VALUE = "")
    groups <- vapply(X = termsAndGroups, FUN = "[", 2L, FUN.VALUE = "")
    groups[is.na(groups)] <- "joint" # season not group-specific

    ## evaluate sin/cos terms using the time vector 't' defined above
    seasonTerms <- mapply(function (coef, expr) coef * eval(parse(text = expr)),
                          coef = coefSinCos, expr = terms, SIMPLIFY = TRUE)

    ## sum sin/cos terms by group
    rowSumsBy.matrix(seasonTerms, groups)
}
