################################################################################
### Highlight epochs in a time series plot of an "sts" object
###
### Copyright (C) 2015 Sebastian Meyer
###
### This file is part of the R package "hhh4contacts",
### free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at http://www.r-project.org/Licenses/.
################################################################################

##' Hook functions for \code{stsplot_time1}
##'
##' Hook functions can be passed to \code{\link[surveillance]{stsplot_time1}},
##' which are evaluated after all the plotting has been done,
##' and with the hook function environment set to the evaluation environment
##' of \code{stsplot_time1} such that local variables can be accessed.
##' They are not intended to be called directly.
##' @name stsplothook
##' @export
##' @param christmas logical indicating if Christmas should be highlighted.
##' @param epochInYear integer vector of epochs to highlight.
##' @param col,lwd graphical parameters for the highlighting lines.
##' @author Sebastian Meyer
##' @examples
##' plot(noroBE("agegroups"), hookFunc = stsplothook_highlight(epochInYear=51))

stsplothook_highlight <- function (christmas = FALSE, epochInYear = NULL,
                                   col = 2, lwd = 2)
{
    indicators <- list(
        isChristmas = quote(isChristmas <- if (identical(x@freq, 52)) {
            epochInYear(x) == strftime(paste0(year(x), "-12-24"), "%V")
            ## with %W, strftime("2013-12-24", "%W") is 51 ?!?
        } else FALSE),
        isInKW = substitute(isInKW <-
            epochInYear(x) %in% KW, list(KW = epochInYear)
        )
    )[c(isTRUE(christmas), !is.null(epochInYear))]

    if (length(indicators) == 0L) {
        return(as.function(alist(invisible()), envir = .GlobalEnv))
    }

    as.function(list(as.call(c(
        as.name("{"),
        indicators,
        parse(text = paste0("toMark <- ",
                  paste0(names(indicators), collapse = " | ")))[[1L]],
        substitute(lines(seq_along(observed)[toMark], observed[toMark],
                         type = "h", lend = 1, lwd = lwd, col = col))
    ))), envir = parent.frame())
}
