% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stsplothook.R
\name{stsplothook}
\alias{stsplothook}
\alias{stsplothook_highlight}
\title{Hook functions for \code{stsplot_time1}}
\usage{
stsplothook_highlight(christmas = FALSE, epochInYear = NULL, col = 2, lwd = 2)
}
\arguments{
\item{christmas}{logical indicating if Christmas should be highlighted.}

\item{epochInYear}{integer vector of epochs to highlight.}

\item{col, lwd}{graphical parameters for the highlighting lines.}
}
\description{
Hook functions can be passed to \code{\link[surveillance]{stsplot_time1}},
which are evaluated after all the plotting has been done,
and with the hook function environment set to the evaluation environment
of \code{stsplot_time1} such that local variables can be accessed.
They are not intended to be called directly.
}
\examples{
plot(noroBE("agegroups"), hookFunc = stsplothook_highlight(epochInYear=51))
}
\author{
Sebastian Meyer
}
