% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/contactmatrix.R
\name{aggregateC}
\alias{aggregateC}
\title{Aggregate a Contact Matrix}
\usage{
aggregateC(C, grouping, ..., weights = NULL)
}
\arguments{
\item{C}{a square numeric contact matrix such as
\code{\link{contactmatrix_POLYMOD}}.}

\item{grouping}{specification of how to aggregate groups (a named list or
an integer vector of group sizes) using \code{\link{aggregateC}} with the
\code{"agedistri"} attribute of the contact matrix as weights.
If \code{NULL}, the original 5-year intervals are returned.
The default setting produces the six age groups of Meyer and Held (2017).}

\item{...}{specification of how to aggregate groups (alternative to using a
named list as the \code{grouping} argument).}

\item{weights}{a named numeric vector containing the weights for the rows of
\code{C}, typically the age distribution of the participants. The names
are matched against \code{rownames(C)}. A value of \code{NULL} is
interpreted as uniform weights.}
}
\description{
The (age) groups of a contact matrix can be joined together by the
\code{grouping} argument, which first sums over contact groups (columns) and
then averages over the corresponding participant groups (rows), optionally
using weights such as the age distribution of the study participants.
}
\author{
Sebastian Meyer
}
\keyword{manip}
