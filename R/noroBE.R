################################################################################
### Create "sts" objects for the Berlin norovirus data
###
### Copyright (C) 2014-2017 Sebastian Meyer
###
### This file is part of the R package "hhh4contacts",
### free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at http://www.r-project.org/Licenses/.
################################################################################

##' Create \code{"sts"} Objects from the Berlin Norovirus Data
##'
##' The function \code{noroBE()} creates an \code{"\link[surveillance:sts-class]{sts}"} object
##' based on the array of norovirus surveillance \code{counts}, the \code{map}
##' of Berlin's city district, and the \code{\link{pop2011}} data stored in the
##' package. This is the data analysed by Meyer and Held (2017).
##'
##' @param by character string determining the stratification, i.e., which units
##'     the resulting \code{"sts"} object should represent: \describe{
##' \item{"districts":}{aggregates \code{counts} and \code{pop2011} over the
##'     age groups and stores the matrix of adjacency orders from the \code{map}
##'     in the \code{neighbourhood} slot. The latter is obtained via
##'     \code{\link[surveillance]{nbOrder}(\link[surveillance]{poly2adjmat}(map), maxlag = 5)}.}
##' \item{"agegroups":}{aggregates \code{counts} and \code{pop2011} over the
##'     districts and stores the \code{\link{contactmatrix}()} in the
##'     \code{neighbourhood} slot, potentially also combining some age groups
##'     via the \code{agegroups} argument.}
##' \item{"all":}{retains both dimensions, either as a list of spatial
##'     \code{"sts"} objects per age group, or in a single \code{"sts"} object
##'     (see \code{flatten} below).}
##' \item{"none":}{creates the overall (univariate) time series of
##'     \code{rowSums(counts)}.}
##' }
##'
##' @param agegroups how the age groups in \code{counts} (and \code{pop2011})
##'     should be aggregated. Will be used as the \code{grouping} argument in
##'     \code{\link{aggregateCountsArray}} and \code{\link{contactmatrix}}.
##'     The default setting uses the six age groups of Meyer and Held (2017).
##' @param timeRange character vector of length two determining the time range
##'     of the \code{"sts"} object to generate. The two strings are matched
##'     against \code{dimnames(counts)[[1]]}, which ranges from
##'     \code{"2011-w01"} until \code{"2016-w30"}. The default value extracts
##'     four seasons (years) starting at \code{"2011-w27"}.
##' @param flatten logical indicating whether for \code{by = "all"} a single
##'     \code{"sts"} object should be returned where the observation unit is the
##'     interaction of district and age group (\dQuote{flattened} \code{counts}
##'     array, see \code{\link{as.data.frame.array}}). By default (\code{flatten
##'     = FALSE}), a list of district-based \code{"sts"} objects is returned,
##'     one for each age group.
##'
##' @format \describe{
##' \item{counts:}{an integer-valued array of norovirus surveillance counts
##' with labelled dimensions of size
##' 290 (\code{"week"}) x 12 (\code{"district"}) x 15 (\code{"agegroup"}).}
##' \item{map:}{a \code{"\link[sp:SpatialPolygonsDataFrame-class]{SpatialPolygonsDataFrame}"}
##' of length 12 with \code{row.names(map)} matching \code{colnames(counts)},
##' representing Berlin's city districts in longlat coordinates (WGS84).
##' The data slot contains the full \code{"NAME"}s of the city districts
##' as well as their \code{"POPULATION"}, i.e., \code{rowSums(pop2011)}.}
##' }
##' The function \code{noroBE()} returns an \code{"\link[surveillance:sts-class]{sts}"} object
##' generated from these data (and \code{\link{pop2011}}).
##' @source \describe{
##' \item{counts:}{based on norovirus surveillance counts retrieved from
##'   the SurvStat@RKI 2.0 online service (\url{https://survstat.rki.de})
##'   of Germany's public health institute, the Robert Koch Institute,
##'   as of 2016-09-08.}
##' \item{map:}{based on a KML file of Berlin's 97 local centres
##'   (\dQuote{Ortsteile}) downloaded from the Berlin Open Data repository at
##'   \url{https://daten.berlin.de/datensaetze/geometrien-der-ortsteile-von-berlin-juli-2012}
##'   as of 2014-11-12, published by
##'   \emph{Amt fuer Statistik Berlin-Brandenburg}
##'   (Statistical Office of Berlin-Brandenburg)
##'   under the \sQuote{CC BY 3.0 DE} license
##'   (\url{https://creativecommons.org/licenses/by/3.0/de/}).
##'   The \code{map} included here aggregates
##'   these local centres by city district.}
##' }
##' @author Sebastian Meyer
##' @references
##' Meyer S and Held L (2017): Incorporating social contact data in
##' spatio-temporal models for infectious disease spread.
##' \emph{Biostatistics}, \bold{18} (2), 338-351.
##' \doi{10.1093/biostatistics/kxw051}
##' @examples
##' ## the raw data
##' str(counts)
##' summary(map)
##'
##' ## district-specific time series
##' noroBEr <- noroBE(by = "districts")
##' plot(noroBEr)
##'
##' ## age group-specific time series
##' noroBEg <- noroBE(by = "agegroups")
##' plot(noroBEg)
##'
##' ## list of spatio-temporal surveillance counts, one for each age group
##' noroBErbyg <- noroBE(by = "all", flatten = FALSE)
##' plot(noroBErbyg[[1L]], par.list = list(oma=c(0,0,2,0)))
##' title(main = names(noroBErbyg)[1], outer = TRUE, line = -1)
##'
##' ## flattened "sts" object (the 'neighbourhood' only reflects spatial info)
##' noroBEall <- noroBE(by = "all", flatten = TRUE)
##' dev.new(width = 16, height = 7)
##' plot(noroBEall, par.list = list(
##'     xaxt = "n", mar = c(1,4,1,1), mfrow = c(ncol(noroBEg), ncol(noroBEr))
##' ))
##' @importFrom methods new
##' @importClassesFrom surveillance sts
##' @export
noroBE <- function (by = c("districts", "agegroups", "all", "none"),
                    agegroups = c(1, 2, 2, 4, 4, 2),
                    timeRange = c("2011-w27", "2015-w26"),
                    flatten = FALSE)
{
    by <- match.arg(by)

    ## subset counts to time range
    stopifnot(!is.na(
        timeRangeIdx <- match(timeRange, dimnames(hhh4contacts::counts)[[1L]])
    ))
    counts <- hhh4contacts::counts[timeRangeIdx[1L]:timeRangeIdx[2L], , , drop = FALSE]
    start <- as.integer(
        strsplit(dimnames(counts)[[1L]][1], "-w", fixed = TRUE)[[1L]]
    )

    ## population numbers
    pop2011 <- hhh4contacts::pop2011

    ## aggregate groups
    if (by %in% c("agegroups", "all") && !is.null(agegroups)) {
        counts <- aggregateCountsArray(counts = counts, dim = 3, grouping = agegroups)
        pop2011 <- aggregateCountsArray(counts = pop2011, dim = 2, grouping = agegroups)
    }

    ## neighbourhood structure
    neighbourhood <- if (by %in% c("districts", "all")) {
        map_nbOrder
    } else if (by == "agegroups") {
        contactmatrix(grouping = agegroups)
    } else { # no stratification
        matrix(NA, 1L, 1L)
    }

    ## constructor function
    makests <- function (observed, populationFrac, neighbourhood, ...) {
        new("sts", start = start, freq = 52,
            observed = observed, populationFrac = populationFrac,
            neighbourhood = neighbourhood, ...)
    }

    ## create "sts" object for the chosen level of stratification
    if (by == "none") {
        makests(
            observed = cbind(Berlin = rowSums(counts)),
            populationFrac = 1, neighbourhood = neighbourhood
        )
    } else if (by == "districts") {
        makests(
            observed = rowSums(counts, dims = 2L),
            populationFrac = prop.table(rowSums(pop2011)),
            map = hhh4contacts::map, neighbourhood = neighbourhood
        )
    } else if (by == "agegroups") {
        makests(
            observed = apply(counts, c(1,3), sum),
            populationFrac = prop.table(colSums(pop2011)),
            neighbourhood = neighbourhood
        )
    } else if (flatten) { # where districts vary faster than age groups
        ngroups <- dim(counts)[[3L]]
        nregions <- dim(counts)[[2L]]
        ## replicate 'neighbourhood' block ngroup times
        ## replicate columns
        neighbourhood <- rep.int(neighbourhood, ngroups)
        dim(neighbourhood) <- c(nregions, nregions * ngroups)
        ## replicate rows
        neighbourhood <- do.call("rbind", rep_len(
            list(as.name("neighbourhood")), ngroups))
        makests(
            observed = as.data.frame(counts),
            populationFrac = c(pop2011),
            neighbourhood = neighbourhood
        )
    } else { # list of by="districts" sts objects by age group
        sapply(dimnames(counts)[[3L]], function (g) {
            map <- hhh4contacts::map
            map$POPULATION <- pop2011[, g]
            makests(
                observed = counts[, , g],
                populationFrac = prop.table(pop2011[, g]),
                map = map, neighbourhood = neighbourhood
            )
        }, simplify = FALSE, USE.NAMES = TRUE)
    }
}

##' @rdname noroBE
##' @format NULL
"counts"

##' @rdname noroBE
##' @format NULL
##' @rawNamespace
##' ## 'map' uses external S4 class: force 'sp' (see WRE 1.1.6)
##' importClassesFrom(sp, SpatialPolygonsDataFrame)  # anything
##'
"map"


## hard-coded adjacency orders to avoid re-calculation (and "spdep" dependency)
## dput(nbOrder(poly2adjmat(hhh4contacts::map), maxlag = 5))
map_nbOrder <- structure(
    c(0L, 2L, 3L, 4L, 1L, 2L, 2L, 1L, 1L, 1L, 1L, 3L, 2L,
      0L, 1L, 2L, 1L, 1L, 1L, 2L, 3L, 2L, 1L, 1L, 3L, 1L, 0L, 1L, 2L,
      2L, 1L, 2L, 3L, 3L, 2L, 1L, 4L, 2L, 1L, 0L, 3L, 2L, 2L, 3L, 4L,
      4L, 3L, 1L, 1L, 1L, 2L, 3L, 0L, 2L, 1L, 1L, 2L, 2L, 1L, 2L, 2L,
      1L, 2L, 2L, 2L, 0L, 2L, 3L, 3L, 2L, 1L, 1L, 2L, 1L, 1L, 2L, 1L,
      2L, 0L, 1L, 2L, 3L, 2L, 2L, 1L, 2L, 2L, 3L, 1L, 3L, 1L, 0L, 1L,
      2L, 2L, 3L, 1L, 3L, 3L, 4L, 2L, 3L, 2L, 1L, 0L, 1L, 2L, 4L, 1L,
      2L, 3L, 4L, 2L, 2L, 3L, 2L, 1L, 0L, 1L, 3L, 1L, 1L, 2L, 3L, 1L,
      1L, 2L, 2L, 2L, 1L, 0L, 2L, 3L, 1L, 1L, 1L, 2L, 1L, 2L, 3L, 4L,
      3L, 2L, 0L),
    .Dim = c(12L, 12L),
    .Dimnames = list(
        c("chwi", "frkr", "lich", "mahe", "mitt", "neuk", "pank", "rein",
          "span", "zehl", "scho", "trko"),
        c("chwi", "frkr", "lich", "mahe", "mitt", "neuk", "pank", "rein",
          "span", "zehl", "scho", "trko")
    ))
