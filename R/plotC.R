################################################################################
### Produce a levelplot or a filled.contour plot of a contact matrix
###
### Copyright (C) 2014-2017 Sebastian Meyer
###
### This file is part of the R package "hhh4contacts",
### free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at http://www.r-project.org/Licenses/.
################################################################################

##' Generate an Image of a Contact Matrix
##'
##' @param C a square numeric matrix.
##' @param grouping numeric vector of sizes of aggregated groups, e.g.,
##' \code{grouping = c(1,3)}, to draw separation lines after the first and
##' the forth subgroup. This is ignored if \code{contour = TRUE}.
##' @param xlab,ylab axis labels.
##' @param at numeric vector of break points of the color levels, or a single
##' integer specifying the number of \code{cuts} (which defaults to 15 as in
##' \code{\link[lattice]{levelplot}}).
##' @param col.regions vector of color levels.
##' @param ... further arguments passed to \code{\link[lattice]{levelplot}} or
##' \code{\link{filled.contour}} (if \code{contour = TRUE}).
##' @param contour logical indicating if a \code{\link{filled.contour}} should
##' be drawn instead of a \code{\link[lattice]{levelplot}} (the default).
##' @examples
##' ## contour plot
##' plotC(contactmatrix_POLYMOD, contour = TRUE)
##'
##' ## level plots illustrating aggregation of age groups
##' if (require("gridExtra")) {
##'     grid.arrange(plotC(contactmatrix_POLYMOD, grouping = c(1,2,2,4,4,2)),
##'                  plotC(contactmatrix(grouping = c(1,2,2,4,4,2))),
##'                  nrow = 1)
##' }
##'
##' @importFrom graphics filled.contour
##' @importFrom grDevices heat.colors
##' @importFrom stats update
##' @export
plotC <- function (C, grouping = NULL,
                   xlab = "age group of contact",
                   ylab = "age group of participant",
                   at = 15, col.regions = rev(heat.colors(length(at)-1)), ...,
                   contour = FALSE)
{
    if (length(at) == 1)
        at <- seq(0, max(C), length.out = at + 2)
    if (contour) {
        agegrid <- seq(0, 74, length.out = nrow(C))
        filled.contour(x = agegrid, y = agegrid, z = t(C),
                       levels = at, col = col.regions, ...)
    } else {
        tr <- lattice::levelplot(t(C), xlab = xlab, ylab = ylab,
                                 at = at, col.regions = col.regions, ...)
        if (is.vector(grouping, mode = "numeric")) {
            update(tr, panel = function (...) {
                cl <- sys.call()
                cl[[1]] <- match.fun(tr$panel)
                eval.parent(cl)
                at <- cumsum(grouping) + 0.5
                lattice::panel.refline(h = at, lwd = 2)
                lattice::panel.refline(v = at, lwd = 2, lty = 2)
            })
        } else {
            tr
        }
    }
}
