################################################################################
### Convert a vector specification of a "grouping" to its list representation
###
### Copyright (C) 2014-2015 Sebastian Meyer
###
### This file is part of the R package "hhh4contacts",
### free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at http://www.r-project.org/Licenses/.
################################################################################

## Handle two ways of specifying 'grouping'
grouping2list <- function (grouping, groups)
{
    if (is.list(grouping)) {
        stopifnot(!is.null(names(grouping)),
                  unlist(grouping) %in% groups,
                  !anyDuplicated(unlist(grouping)))
    } else {
        stopifnot(is.vector(grouping, mode = "numeric"), !is.na(grouping))
        withNames <- !is.null(names(grouping))
        remainder <- length(groups) - sum(grouping)
        if (remainder > 0) {
            grouping <- c(grouping, remainder)
            if (withNames) {
                names(grouping)[length(grouping)] <- "other"
            }
        } else if (remainder < 0) {
            stop("sum(grouping) > length(groups)")
        }
        lg <- rep.int(if (withNames) names(grouping) else seq_along(grouping),
                      grouping)
        grouping <- split.default(x = groups, f = lg)[grouping > 1]
        if (!withNames) { # automatic names (only suitable for age groups)
            names(grouping) <- vapply(
                X = grouping,
                FUN = function (groupnames) {
                    L <- length(groupnames)
                    if (L == 1) {
                        groupnames
                    } else {
                        paste0(substr(groupnames[1L], 1, 2),
                               if (grepl("\\+$", groupnames[L])) "+" else
                               substr(groupnames[L], 3, 5))
                    }
                },
                FUN.VALUE = "", USE.NAMES = FALSE)
        }
    }

    return(grouping)
}
