% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotC.R
\name{plotC}
\alias{plotC}
\title{Generate an Image of a Contact Matrix}
\usage{
plotC(
  C,
  grouping = NULL,
  xlab = "age group of contact",
  ylab = "age group of participant",
  at = 15,
  col.regions = rev(heat.colors(length(at) - 1)),
  ...,
  contour = FALSE
)
}
\arguments{
\item{C}{a square numeric matrix.}

\item{grouping}{numeric vector of sizes of aggregated groups, e.g.,
\code{grouping = c(1,3)}, to draw separation lines after the first and
the forth subgroup. This is ignored if \code{contour = TRUE}.}

\item{xlab, ylab}{axis labels.}

\item{at}{numeric vector of break points of the color levels, or a single
integer specifying the number of \code{cuts} (which defaults to 15 as in
\code{\link[lattice]{levelplot}}).}

\item{col.regions}{vector of color levels.}

\item{...}{further arguments passed to \code{\link[lattice]{levelplot}} or
\code{\link{filled.contour}} (if \code{contour = TRUE}).}

\item{contour}{logical indicating if a \code{\link{filled.contour}} should
be drawn instead of a \code{\link[lattice]{levelplot}} (the default).}
}
\description{
Generate an Image of a Contact Matrix
}
\examples{
## contour plot
plotC(contactmatrix_POLYMOD, contour = TRUE)

## level plots illustrating aggregation of age groups
if (require("gridExtra")) {
    grid.arrange(plotC(contactmatrix_POLYMOD, grouping = c(1,2,2,4,4,2)),
                 plotC(contactmatrix(grouping = c(1,2,2,4,4,2))),
                 nrow = 1)
}

}
