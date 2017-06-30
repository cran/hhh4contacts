################################################################################
### R script to reproduce the main results of
###
### Meyer and Held (2017): "Incorporating social contact data in
### spatio-temporal models for infectious disease spread."
### Biostatistics, 18 (2), pp. 338-351. DOI: 10.1093/biostatistics/kxw051
###
### Copyright (C) 2016-2017  Sebastian Meyer <seb.meyer@fau.de>
###
### This script can be redistributed and/or modified under the terms
### of the GNU General Public License, version 2 or later.
################################################################################

## Load the package dedicated to the manuscript.
## It contains the data and some specialized functions for age-structured
## spatio-temporal "hhh4" models based on the R package surveillance.
library("hhh4contacts")



### SURVEILLANCE DATA
## Load counts of norovirus gastroenteritis stratified by age group and district
## in Berlin, 2011--2015, as obtained from SurvStat@RKI.
###


## six age groups aggregated from original 5-year intervals
GROUPING <- c(1, 2, 2, 4, 4, 2)

## "sts" object where the units are "district.agegroup" (district varies faster)
noroBEall <- noroBE(by = "all", flatten = TRUE, agegroups = GROUPING)
## and the same as list of NGROUPS "sts" objects
noroBErbyg <- noroBE(by = "all", flatten = FALSE, agegroups = GROUPING)

## an "sts" object aggregated over all age groups
noroBEr <- noroBE(by = "districts")

## matrix of corresponding population counts (district x age group)
(popBErbyg <- aggregateCountsArray(pop2011, dim = 2, grouping = GROUPING))

## names of groups and districts
(DISTRICTS <- unique(stratum(noroBEall, 1)))
NDISTRICTS <- length(DISTRICTS)
(GROUPS <- unique(stratum(noroBEall, 2)))
NGROUPS <- length(GROUPS)



### plot the surveillance time series

## overall time series of counts
plot(noroBEall, type = observed ~ time,
     main = "Aggregated over all districts and age groups")

## group-specific time series
par(mfrow = c(2,3), las = 1)
for (g in GROUPS)
    plot(noroBErbyg[[g]], type = observed ~ time, main = g, ylim = c(0,130))

## district-specific time series
plot(noroBEr)



### plot disease incidence maps (mean yearly incidence)

## overall map
scalebar <- layout.scalebar(noroBErbyg[[1]]@map, corner = c(0.7, 0.9),
    scale = 10, labels = c(0, "10 km"), cex = 0.6, height = 0.02)
plot(noroBEr, type = observed ~ unit,
     main = "Aggregated over all age groups",
     population = rowSums(popBErbyg)/100000 * (nrow(noroBEall)/52), # per year
     labels = list(cex = 0.8), sp.layout = scalebar)

## group-specific maps
incidence_maps <- lapply(GROUPS, function (g) {
    stsObj <- noroBErbyg[[g]]
    plot(stsObj, type = observed ~ unit, main = g,
         population = popBErbyg[,g]/100000 * (nrow(stsObj)/52),
         labels = list(cex = 0.6), sp.layout = scalebar)
})
if (require("gridExtra")) {
    grid.arrange(grobs = incidence_maps, nrow = 2, ncol = 3)
} else {
    warning("install package \"gridExtra\" to plot all maps on one page")
    plot(incidence_maps[[1]])
}






### CONTACT DATA
## Load the age-structured contact matrix estimated from the German POLYMOD data
###


## aggregate to the same six age groups as the above counts
Cgrouped <- contactmatrix(
    which = "reciprocal", # estimated by the Wallinga et al (2006) method
    type = "all",         # alternatively: "physical" only
    grouping = GROUPING   # age group specification
)
Cgrouped
## the "agedistri" attribute corresponds to the population fractions in Berlin

## the row-normalized version
(Cgrouped_norm <- Cgrouped / rowSums(Cgrouped))

## for the no-mixing model: diagonal contact matrix (no mixing)
(Cgrouped_AR <- structure(diag(ncol(Cgrouped)), dimnames = dimnames(Cgrouped)))



### plot the contact matrix

## contact matrix with 5-year intervals
plotC(contactmatrix("reciprocal", grouping = NULL),
      grouping = GROUPING, scales = list(x = list(rot = 45, cex = 0.7)))

## aggregated to the six age groups
plotC(Cgrouped)

## power-adjustment of the contact matrix
powerC <- make_powerC(Cgrouped_norm, truncate = TRUE)
powerC(0.5)
if (require("gridExtra"))
    grid.arrange(grobs = lapply(c(0, 0.5, 1, 2), function (power)
        plotC(round(powerC(power), 7), main = bquote(kappa == .(power)),
              at = seq(0, 1, length.out = 15 + 2))))






### MODEL FITS
## We estimate various hhh4() models with spatial power-law weights,
## population gravity and (power-adjusted) age-structured contact matrix
###

## we want the power law to act on (o+1)
neighbourhood(noroBEall) <- neighbourhood(noroBEall) + 1


### setup time variables and indicators for group/district-specific effects

DATAt <- list(t = epoch(noroBEall) - 1,
              christmas = 1*(epochInYear(noroBEall) %in% c(52, 1)))

## setup a model matrix with group indicators
MMG <- sapply(GROUPS, function (g) {
    index <- which(stratum(noroBEall, which = 2) == g)
    res <- col(noroBEall)
    res[] <- res %in% index
    res
}, simplify = FALSE, USE.NAMES = TRUE)
str(MMG)

## setup model matrix with district indicators
MMR <- sapply(DISTRICTS, function (r) {
    index <- which(stratum(noroBEall, which = 1) == r)
    res <- col(noroBEall)
    res[] <- res %in% index
    res
}, simplify = FALSE, USE.NAMES = TRUE)
str(MMR)

## setup model matrix of group-specific seasonal terms
MMgS <- with(c(MMG, DATAt), unlist(lapply(
    X = GROUPS,
    FUN = function (g) {
        gIndicator <- get(g)
        res <- list(gIndicator * sin(2 * pi * t/52),
                    gIndicator * cos(2 * pi * t/52))
        names(res) <- paste0(c("sin", "cos"), "(2 * pi * t/52).", g)
        res
    }), recursive = FALSE, use.names = TRUE))
str(MMgS)


### specify the basic endemic model

## endemic formula: ~group + district + christmas + group:(sin+cos)
qGROUPS <- paste0("`", GROUPS, "`")
(FGRXgS <- reformulate(c(qGROUPS[-1], DISTRICTS[-1], "christmas",
                         paste0("`", names(MMgS), "`")),
                       intercept = TRUE))
control0 <- list(
    end = list(f = FGRXgS,
        offset = population(noroBEall) / rowSums(population(noroBEall))),
    family = factor(stratum(noroBEall, which = 2)), # group-specific dispersion
    data = c(MMG, MMR, DATAt, MMgS)
)

## fit the endemic-only model
ma0 <- hhh4(noroBEall, control0)


### add epidemic component with population gravity and POLYMOD contacts

## epidemic formula: ~group + district + log(pop)
(FGRpop <- reformulate(c(qGROUPS[-1], DISTRICTS[-1], "log(pop)"),
                       intercept = TRUE))

## fit the power-law model with the given contact matrix
ma_popPLC <- update(ma0,
    ne = list(
        f = FGRpop,
        weights = W_powerlaw(maxlag = 5, log = TRUE, normalize = FALSE,
                             initial = c("logd" = log(2))),
        scale = expandC(Cgrouped_norm, NDISTRICTS),
        normalize = TRUE),
    data = list(pop = population(noroBEall)/rowSums(population(noroBEall))))

## model summary
summary(ma_popPLC, maxEV = TRUE, reparamPsi = TRUE,
        amplitudeShift = TRUE, idx2Exp = TRUE)



### alternative spatial transmission weights
### (here still without power-adjustment of C)

## a group-specific power law
(rho <- coef(ma_popPLC, se = TRUE, idx2Exp = TRUE)["exp(neweights.logd)",])
PLgroups <- factor(stratum(noroBEall, which = 2))
levels(PLgroups)[2:3] <- paste0(levels(PLgroups)[2:3], collapse = " & ")
levels(PLgroups)
WPLgroups <- addGroups2WFUN(WFUN = ma_popPLC$control$ne$weights,
                            groups = PLgroups)
ma_popGPLC <- update(ma_popPLC, ne = list(weights = WPLgroups))
rho_groups <- coef(ma_popGPLC, se = TRUE, idx2Exp = TRUE)[
    paste0("exp(neweights.", levels(PLgroups), ")"),]
rho_groups

## unconstrained spatial transmission weights
ma_popNPC <- update(ma_popPLC,
    ne = list(weights = W_np(maxlag = 5, normalize = FALSE)))

## plot the power laws and the unconstrained weights
par(mfrow = c(1,1))
plot(1:5, (1:5)^-rho[1], type = "b", lwd = 3, pch = 19,
     xlab = "adjacency order", ylab = "weight", xaxt = "n")
axis(1, at = 1:5, labels = 0:4)
matlines(1:5, sapply(rho_groups[,1], function (rho) (1:5)^-rho),
         type = "l", lwd = 2, lty = 2:6, col = 2:6)
points(1:5, c(1, exp(coefW(ma_popNPC))), pch = 15, lwd = 2)
## Note: these were still models without power-adjustment of C,
##       so the results are different from what is reported in the paper



### fit models with alternative contact structures

## no mixing
ma_popPLAR <- update(ma_popPLC,
    ne = list(scale = expandC(Cgrouped_AR, NDISTRICTS)),
    use.estimates = FALSE)

## homogeneous mixing
ma_popPLhom <- update(ma_popPLC,
    ne = list(scale = NULL), # C = 1
    use.estimates = FALSE)

## AIC comparison
AIC(ma0, ma_popPLhom, ma_popPLAR, ma_popPLC)



### fit power-adjusted contact matrix via profile likelihood
### CAVE: this takes a while (approx. 3 minutes)

ma_popPLCpower <- fitC(ma_popPLC, Cgrouped, normalize = TRUE, truncate = TRUE)

## AIC comparison
AIC(ma_popPLC, ma_popPLCpower)

## model summary
summary(ma_popPLCpower, maxEV = TRUE, reparamPsi = TRUE,
        amplitudeShift = TRUE, idx2Exp = TRUE)



### plot fitted values

## overall fit
plotHHH4_fitted_groups(ma_popPLCpower,
    groups = stratum(noroBEall, which = 2), total = TRUE, pch = 20,
    legend.args = list(legend = c("from other age groups", "within age group", "endemic")))

## by age group
plotHHH4_fitted_groups(ma_popPLCpower,
    groups = stratum(noroBEall, which = 2), units = NULL, pch = 20, legend = 2,
    legend.args = list(legend = c("from other age groups", "within age group", "endemic")))

## by district
plotHHH4_fitted_groups(ma_popPLCpower,
    groups = factor(stratum(noroBEall, which = 1), levels = DISTRICTS),
    names = noroBEr@map@data[DISTRICTS,"NAME"], units = NULL,
    legend = 4, legend.args = list(cex = 0.8,
        legend = c("from other districts", "within district", "endemic")),
    pch = 20, ylim = c(0,45))



### plot endemic seasonality

set.seed(131015)  # confidence intervals involve simulations
plotHHH4_season_groups(ma_popPLCpower,
    component = "end", seasonStart = 27,
    col = c("#D53E4F", "#FC8D59", "#FEE08B", "#E6F598", "#99D594", "#3288BD"),
    xlab = "calendar week", ylab = "multiplicative effect",
    xaxt = "n", xaxs = "i", yaxs = "i", ylim = c(0, 5.5))
## add x-axis
weeks2mark <- seq.int(27, 27+52, by = 4)  # axTicks(1)
axis(1, at = 0:100, labels = FALSE, tcl = NA)
axis(1, at = weeks2mark, labels = weeks2mark %% 52)
