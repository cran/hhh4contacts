################################################################################
### Generic function "stratum" to extract labels of the grouping variables
###
### Copyright (C) 2015 Sebastian Meyer
###
### This file is part of the R package "hhh4contacts",
### free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at http://www.r-project.org/Licenses/.
################################################################################

##' Extract Strata
##'
##' Methods to extract strata information from an object.
##' Here we only define a method for class \code{"sts"}.
##' @param x an object of class \code{"sts"}.
##' @param ... further arguments passed to methods.
##' @return a character vector of strata names of length \code{ncol(x)}.
##' @export
setGeneric("stratum", function (x, ...) standardGeneric("stratum"))

##' @param which an integer (strata dimension) or \code{NULL} (to get the plain
##' \code{colnames}, the default).
##' @describeIn stratum
##' Extract the names of the units, i.e., the \code{colnames},
##' from a multivariate \code{"sts"} object.
##' If the units result from the interaction of multiple strata
##' separated by dots, e.g., \code{"region.group"},
##' the function can also extract the names corresponding to a specific
##' strata dimension, e.g., \code{which = 2} to get the group names.
##' @importClassesFrom surveillance sts
##' @examples
##' noroBEall <- noroBE(by = "all", flatten = TRUE)
##' stratum(noroBEall)  # just colnames(noroBEall)
##' stratum(noroBEall, which = 2)  # the age groups
setMethod("stratum", "sts", function (x, which = NULL, ...)
    {
        units <- colnames(x)
        if (is.null(which))
            return(units)

        unlist(lapply(X = strsplit(units, split = ".", fixed = TRUE),
                      FUN = "[[", which),
               recursive = FALSE, use.names = FALSE)
    })
