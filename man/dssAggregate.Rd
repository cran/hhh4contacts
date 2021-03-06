% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dssAggregate.R
\name{dssAggregate}
\alias{dssAggregate}
\title{Compute the DSS on Aggregated Predictions and Observations}
\usage{
dssAggregate(observed, pred, psi, groups)
}
\arguments{
\item{observed}{a numeric matrix of observed counts.}

\item{pred}{a numeric matrix of predicted counts.}

\item{psi}{a numeric vector or matrix of overdispersion parameters such that
\code{pred * (1 + pred/exp(psi))} is the prediction's variance.
Alternatively, \code{psi = NULL} indicated Poisson predictions.}

\item{groups}{a factor variable of length \code{ncol(observed)} indicating
which columns should be aggregated.}
}
\value{
a matrix of DSS values
}
\description{
The expectation and variance of aggregated predictions is just a sum if the
predictions are (conditionally) independent. This function computes the DSS
for a matrix of observations and a matrix of predictions
where the columns are to be summed according to a given factor.
}
