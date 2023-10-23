library("hhh4contacts")
library("surveillance")
data("measlesWeserEms")

## basic power-law model
WPL <- W_powerlaw(maxlag = 5, normalize = TRUE)
measlesModel <- list(
    end = list(f = addSeason2formula(~1), offset = population(measlesWeserEms)),
    ar = list(f = ~1),
    ne = list(f = ~1 + log(pop), weights = WPL),
    family = "NegBin1", data = list(pop = population(measlesWeserEms))
)
measlesFit <- hhh4(stsObj = measlesWeserEms, control = measlesModel)

## matrix of adjacency orders
nbmat <- neighbourhood(measlesWeserEms)

## fake group-specific power law (single group)
WPLgfake <- addGroups2WFUN(WPL, groups = factor(rep("d", ncol(nbmat))))
stopifnot(identical(WPL$w(0.5, nbmat), WPLgfake$w(0.5, nbmat)),
          identical(WPL$dw(0.5, nbmat), WPLgfake$dw(0.5, nbmat)[[1]]),
          identical(WPL$d2w(0.5, nbmat), WPLgfake$d2w(0.5, nbmat)[[1]]))
stopifnot(all.equal(
    measlesFit,
    update(measlesFit, ne = list(weights = WPLgfake), use.estimates = FALSE),
    ignore = "control"))

### uncomment below to check derivatives with multiple groups
### (time-consuming and analytical derivatives already verified)
## WPLgroups <- addGroups2WFUN(WPL, factor(substr(colnames(nbmat), 1, 4) == "0345"))
## pars_groups <- c(0.5, 2)
## dwnum <- sapply(seq_along(nbmat), function (i)
##    numDeriv::grad(function (wpar) WPLgroups$w(wpar, nbmat)[i], x = pars_groups))
## stopifnot(all.equal(dwnum[1,], c(WPLgroups$dw(pars_groups, nbmat)[[1]])),
##          all.equal(dwnum[2,], c(WPLgroups$dw(pars_groups, nbmat)[[2]])))
## d2wnum <- sapply(seq_along(nbmat), function (i)
##    numDeriv::hessian(function (wpar) WPLgroups$w(wpar, nbmat)[i], x = pars_groups))
## stopifnot(all.equal(d2wnum[1,], c(WPLgroups$d2w(pars_groups, nbmat)[[1]])),
##          all.equal(d2wnum[4,], c(WPLgroups$d2w(pars_groups, nbmat)[[3]])),
##          abs(c(d2wnum[c(2,3),]) - 0) < .Machine$double.eps)

### check score vector and Fisher info of all model parameters
## measlesModelGrouped <- modifyList(measlesModel, list(ne=list(weights=WPLgroups)))
## hhh4(measlesWeserEms, measlesModelGrouped, check.analyticals = TRUE)
