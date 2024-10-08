% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fitC.R
\name{fitC}
\alias{fitC}
\title{Estimate the Power of the Contact Matrix in a \code{"hhh4"} Model}
\usage{
fitC(object, C, normalize = TRUE, truncate = TRUE, optim.args = list(), ...)
}
\arguments{
\item{object}{a model fit of class \code{"\link[surveillance]{hhh4}"}.}

\item{C}{the contact matrix to use.}

\item{normalize, truncate}{see \code{\link{powerC}}.}

\item{optim.args}{a list to modify the default optimization parameters.}

\item{...}{additional arguments for each run of \code{\link[surveillance]{update.hhh4}}.}
}
\value{
an object of class \code{"fitC"}, which is an \code{"\link[surveillance]{hhh4}"}
    object with an additional element \code{logpower}.
}
\description{
The profile log-likelihood of the log(power) parameter of the contact matrix
(see \code{\link{powerC}}) is maximized using \code{\link{optim}}.
The \code{\link[surveillance]{hhh4}} fit for the optimal power value is returned with an
additional element \code{logpower} which holds information on the result of
the optimization.
}
\author{
Sebastian Meyer
}
