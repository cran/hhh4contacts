################################################################################
### Social contact matrices based on the German POLYMOD sample
###
### Copyright (C) 2014-2017 Sebastian Meyer
###
### This file is part of the R package "hhh4contacts",
### free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at http://www.r-project.org/Licenses/.
################################################################################

##' POLYMOD Contact Matrices for Germany
##'
##' The function \code{contactmatrix} retrieves various social contact matrices
##' for Germany from the POLYMOD survey (Mossong et al., 2008). Such a matrix
##' contains the average numbers of reported contacts by participant age group.
##' The original age groups (5-year intervals) can be joined together by the
##' \code{grouping} argument, which first sums over contact groups (columns) and
##' then averages over the corresponding participant groups (rows) using the
##' corresponding age distribution as weights.
##'
##' @param which character string indicating which contact matrix to return.
##' \code{"mossong"} uses the average numbers of reported contacts as published
##' in Table S5 of Mossong et al. (2008), available as
##' \code{contactmatrix_mossong} or \code{contactmatrix_mossong_physical}.
##' \code{"corrected"} (from \code{contactmatrix_POLYMOD} or
##' \code{contactmatrix_POLYMOD_physical}) fixes an error in these numbers
##' related to the age group 70+ (see the Examples) and is the default.
##' If \code{which="reciprocal"} (corresponding to \code{contactmatrix_wallinga}
##' or \code{contactmatrix_wallinga_physical} as used by Meyer and Held, 2017),
##' the returned social contact matrix fulfils reciprocity of contacts with
##' respect to the age distribution of Berlin, \code{\link{pop2011}},
##' via the method of Wallinga et al. (2006).
##' @param type a character string to select the type of contacts to use:
##' either \code{"all"} contacts, i.e., count both physical and pure
##' conversational contacts, or only \code{"physical"} contacts.
##' @param grouping specification of how to aggregate groups (a named list or
##' an integer vector of group sizes) using \code{\link{aggregateC}} with the
##' \code{"agedistri"} attribute of the contact matrix as weights.
##' If \code{NULL}, the original 5-year intervals are returned.
##' The default setting produces the six age groups of Meyer and Held (2017).
##' @param normalize a logical indicating whether to normalize the matrix
##' such that each row sums to 1.
##' @return a square numeric matrix containing the average numbers of contact
##' persons recorded per day per survey participant in Germany,
##' potentially averaged over multiple \emph{row} (participant) age groups
##' and aggregated over the corresponding \emph{column} (contact) age groups.
##' @format The dataset \code{contactmatrix_POLYMOD} and its variants are all
##' square numeric matrices with 15 rows (participants) and 15 columns
##' (contacts), labelled with the corresponding age groups. There is an attribute
##' \code{"agedistri"}, a named numeric vector of length 15, which for the
##' \dQuote{_mossong_} and \dQuote{_POLYMOD_} variants gives the
##' age distribution of the German POLYMOD sample, and for the \code{_wallinga_}
##' variants gives the age distribution of Berlin, i.e.,
##' \code{prop.table(colSums(\link{pop2011}))}.
##' @source \code{contactmatrix_mossong} and \code{contactmatrix_mossong_physical}
##' are taken from the Supporting Information in Mossong et al. (2008):
##' the matrices from Table S5 (8.2), and the attached age distribution
##' from Table S2 (3.2).
##'
##' The corrected versions \code{contactmatrix_POLYMOD} and
##' \code{contactmatrix_POLYMOD_physical} were constructed from the raw POLYMOD
##' data available at
##' \url{https://www.researchgate.net/publication/232701632_POLYMOD_contact_survey_for_researchers}.
##' The reciprocal contact matrices \code{contactmatrix_wallinga} and
##' \code{contactmatrix_wallinga_physical} were estimated from these raw data
##' via the method of Wallinga et al. (2006).
##' @author Sebastian Meyer
##' @references
##' Meyer S and Held L (2017): Incorporating social contact data in
##' spatio-temporal models for infectious disease spread.
##' \emph{Biostatistics}, \bold{18} (2), 338-351.
##' \doi{10.1093/biostatistics/kxw051}
##'
##' Mossong et al. (2008):
##' Social contacts and mixing patterns relevant to the
##' spread of infectious diseases.
##' \emph{PLoS Medicine}, \bold{5} (3), e74.
##' \doi{http://dx.doi.org/10.1371/journal.pmed.0050074}
##'
##' Wallinga J, Teunis P and Kretzschmar M (2006):
##' Using data on social contacts to estimate age-specific transmission
##' parameters for respiratory-spread infectious agents.
##' \emph{American Journal of Epidemiology}, \bold{164} (10), 936-944.
##' \doi{10.1093/aje/kwj317}
##' @examples
##' ## contact matrix reported in Mossong et al (2008, Table S5)
##' (C_original <- contactmatrix(which = "mossong", grouping = NULL))
##' ## this simply returns the dataset 'contactmatrix_mossong'
##' stopifnot(identical(C_original, contactmatrix_mossong))
##'
##' ## with corrected numbers for the 70+ age group (the default)
##' C_corrected <- contactmatrix(which = "corrected", grouping = NULL)
##' ## this simply returns the dataset 'contactmatrix_POLYMOD'
##' stopifnot(identical(C_corrected, contactmatrix_POLYMOD))
##'
##' ## check for differences
##' C_original == round(C_corrected, 2)
##' ## compare entries of last row and last column
##' round(rbind(original = C_original[15,], corrected = C_corrected[15,]), 2)
##' round(cbind(original = C_original[,15], corrected = C_corrected[,15]), 2)
##'
##' ## contact matrix estimated to be reciprocal on the population level
##' C_reciprocal <- contactmatrix(which = "reciprocal", grouping = NULL)
##' ## this simply returns the dataset 'contactmatrix_wallinga'
##' ## (without its "overdisp" attribute)
##' stopifnot(all.equal(C_reciprocal, contactmatrix_wallinga, check.attributes=FALSE))
##'
##' ## check reciprocity
##' agedistriBE <- attr(C_reciprocal, "agedistri")
##' stopifnot(identical(agedistriBE, prop.table(colSums(pop2011))))
##' stopifnot(isSymmetric(C_reciprocal * agedistriBE, check.attributes=FALSE))
##'
##' ## visually compare raw to reciprocal contact matrix
##' if (require("gridExtra"))
##'     grid.arrange(plotC(C_corrected, main = "raw"),
##'                  plotC(C_reciprocal, main = "reciprocal"),
##'                  nrow = 1)
##'
##' ## select physical contacts and aggregate into 5 age groups
##' contactmatrix(type = "physical", grouping = c(1, 2, 7, 3, 2))
##'
##' ## the default 6 age groups, normalized to a transition matrix
##' contactmatrix(normalize = TRUE)
##'
##' ## reciprocity also holds for this grouping
##' (C6 <- contactmatrix(which = "reciprocal"))
##' stopifnot(isSymmetric(C6 * attr(C6, "agedistri"), check.attributes=FALSE))
##'
##' @export
contactmatrix <- function (which = c("corrected", "mossong", "reciprocal"),
                           type = c("all", "physical"),
                           grouping = c(1, 2, 2, 4, 4, 2),
                           normalize = FALSE)
{
    ## Note that "lazy-loaded datasets are not in the package's namespace
    ## so need to be accessed via ::"
    ## (http://cran.r-project.org/doc/manuals/r-release/R-exts.html)
    C <- getExportedValue("hhh4contacts", paste0("contactmatrix_",
        switch(match.arg(which), "corrected" = "POLYMOD",
               "mossong" = "mossong", "reciprocal" = "wallinga"),
        if (match.arg(type) == "physical") "_physical"))
    attr(C, "overdisp") <- NULL
    if (!is.null(grouping)) {
        weights <- attr(C, "agedistri")
        C <- aggregateC(C = C, grouping = grouping, weights = weights)
        attr(C, "agedistri") <- aggregateCountsArray(
            t(weights), dim = 2, grouping = grouping, sort = TRUE)[1L,]
    }
    if (normalize) C/rowSums(C) else C
}

##' @rdname contactmatrix
##' @format NULL
"contactmatrix_mossong"
##' @rdname contactmatrix
##' @format NULL
"contactmatrix_mossong_physical"
##' @rdname contactmatrix
##' @format NULL
"contactmatrix_POLYMOD"
##' @rdname contactmatrix
##' @format NULL
"contactmatrix_POLYMOD_physical"
##' @rdname contactmatrix
##' @format NULL
"contactmatrix_wallinga"
##' @rdname contactmatrix
##' @format NULL
"contactmatrix_wallinga_physical"


##' Aggregate a Contact Matrix
##'
##' The (age) groups of a contact matrix can be joined together by the
##' \code{grouping} argument, which first sums over contact groups (columns) and
##' then averages over the corresponding participant groups (rows), optionally
##' using weights such as the age distribution of the study participants.
##'
##' @inheritParams contactmatrix
##' @param C a square numeric contact matrix such as
##'     \code{\link{contactmatrix_POLYMOD}}.
##' @param ... specification of how to aggregate groups (alternative to using a
##'     named list as the \code{grouping} argument).
##' @param weights a named numeric vector containing the weights for the rows of
##'     \code{C}, typically the age distribution of the participants. The names
##'     are matched against \code{rownames(C)}. A value of \code{NULL} is
##'     interpreted as uniform weights.
##' @author Sebastian Meyer
##' @keywords manip
##' @importFrom stats setNames
##' @export
aggregateC <- function (C, grouping, ..., weights = NULL)
{
    ndn <- names(dimnames(C))
    stopifnot(!is.null(groups <- rownames(C)))
    if (missing(grouping)) {
        grouping <- list(...)
    }
    grouping <- grouping2list(grouping = grouping, groups = groups)
    if (is.null(weights)) {
        weights <- setNames(rep.int(1/length(groups), length(groups)), groups)
    } else {
        stopifnot(is.vector(weights, mode = "numeric"),
                  rownames(C) %in% names(weights))
    }

    ## aggregate C according to each of 'grouping'
    newgroups <- names(grouping)
    for (i in seq_along(grouping)) {
        C <- aggregateC1(C = C, groups = grouping[[i]],
                         weights = weights[grouping[[i]]],
                         name = newgroups[i])
    }
    names(dimnames(C)) <- ndn
    return(C)
}

##' @importFrom stats weighted.mean
aggregateC1 <- function (C, groups, weights, name, sort = TRUE)
{
    idxgroups <- match(groups, rownames(C))
    res <- C[-idxgroups, -idxgroups, drop = FALSE]

    ## mean number of contacts from other groups to 'union(groups)'
    ## is just the sum of the group-specific means
    res <- cbind(res, rowSums(C[-idxgroups, idxgroups, drop = FALSE]),
                 deparse.level = 0)

    ## mean number of contacts from the new joint group to the other groups
    ## is the weighted mean of the group-specific contact numbers
    res <- rbind(res,
        c(apply(C[idxgroups,-idxgroups,drop=FALSE], 2, weighted.mean, w=weights),
          weighted.mean(rowSums(C[idxgroups,idxgroups,drop=FALSE]), w=weights)),
        deparse.level = 0)

    ## set new name
    colnames(res)[ncol(res)] <- rownames(res)[nrow(res)] <- name

    ## sort columns
    if (sort) {
        ord <- order(rownames(res))
        res <- res[ord, ord, drop = FALSE]
    }

    return(res)
}


##' Expand the Contact Matrix over Regions
##'
##' This is simply the Kronecker product of the contact matrix \code{C}
##' with a matrix of ones of dimension \code{n} x \code{n}.
##' @param C a \code{\link{contactmatrix}}.
##' @param n the size of the secondary dimension to expand to.
##' @return a square matrix with \code{n*ncol(C)} rows and columns.
##' @examples
##' expandC(contactmatrix(), 2)
##' @export
expandC <- function (C, n)
{
    .kronecker(C, matrix(1, n, n))
}
