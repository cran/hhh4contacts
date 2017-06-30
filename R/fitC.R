################################################################################
### Profile likelihood inference for the power adjustment of the contact matrix
###
### Copyright (C) 2015-2017 Sebastian Meyer
###
### This file is part of the R package "hhh4contacts",
### free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at http://www.r-project.org/Licenses/.
################################################################################

##' Estimate the Power of the Contact Matrix in a \code{"hhh4"} Model
##'
##' The profile log-likelihood of the log(power) parameter of the contact matrix
##' (see \code{\link{powerC}}) is maximized using \code{\link{optim}}.
##' The \code{\link{hhh4}} fit for the optimal power value is returned with an
##' additional element \code{logpower} which holds information on the result of
##' the optimization.
##'
##' @param object a model fit of class \code{"\link{hhh4}"}.
##' @param C the contact matrix to use.
##' @param normalize,truncate see \code{\link{powerC}}.
##' @param optim.args a list to modify the default optimization parameters.
##' @param ... additional arguments for each run of \code{\link{update.hhh4}}.
##' @author Sebastian Meyer
##' @return an object of class \code{"fitC"}, which is an \code{"\link{hhh4}"}
##'     object with an additional element \code{logpower}.
##' @importFrom stats update optim
##' @importFrom utils modifyList
##' @export
fitC <- function (object, C, normalize = TRUE, truncate = TRUE,
                  optim.args = list(), ...)
{
    ## function that raises C to a power
    powerC <- make_powerC(C, normalize = normalize, truncate = truncate)

    ## refit the model for a given power
    update_power <- function (power) {
        Cpowered <- powerC(power)
        nExpand <- object$nUnit / ncol(Cpowered)
        update.args <- modifyList(list(...),
            list(object = quote(object),
                 ne = list(scale = expandC(Cpowered, nExpand))))
        do.call("update", update.args)
    }

    ## cache profile log-likelihood values during optimization
    .pllcache <- matrix(numeric(0L), 0L, 2L,
                        dimnames = list(NULL, c("logpower", "pll")))

    ## objective function: negative profile log-likelihood of log(power)
    negpll_logpower <- function (logpower) {
        fitHHH <- update_power(exp(logpower))
        pll <- if (fitHHH$convergence) fitHHH$loglikelihood else NA_real_
        .pllcache <<- rbind(.pllcache, c(logpower, pll))
        -pll
    }

    ## do the optimization
    if (getOption("warn") == 0) {
        ## non-convergence warnings from hhh4() should appear directly
        ## below the trace of the corresponding iteration, not at the end
        oopt <- options(warn = 1)
        on.exit(options(oopt))
    }
    optim.args <- modifyList(
        list(par = 0, fn = negpll_logpower, method = "BFGS",
             control = list(trace = 1, REPORT = 1), hessian = TRUE),
        optim.args
    )
    optim_logpower <- do.call("optim", optim.args)

    ## refit model with estimated power
    ## CAVE: cannot simply use the last result of update_power()
    ##       because of additional evaluations for the Hessian
    res <- update_power(power <- exp(logpower <- optim_logpower$par))

    ## attach log(power) results to the final "hhh4" fit
    res$logpower <- list(
        estimate = logpower,
        se = if (!is.null(optim_logpower$hessian))
                 ##try(sqrt(diag(solve(optim_logpower$hessian)))),
                 sqrt(1/optim_logpower$hessian[1L]),
        optim = optim_logpower,
        pll = .pllcache,
        Cpowered = powerC(power),
        args = list(C = C, normalize = normalize, truncate = truncate)
    )
    class(res) <- c("fitC", class(res))
    res
}

##' @export
logLik.fitC <- function (object, ...)
{
    val <- NextMethod("logLik")
    attr(val, "df") <- attr(val, "df") + 1L
    val
}

##' @export
coef.fitC <- function (object, se = FALSE, ...)
{
    coefs <- NextMethod("coef")
    power <- exp(object$logpower$estimate)
    attr(coefs, "power") <- if (se) {
        c("Estimate" = power, "Std. Error" = power * object$logpower$se)
    } else power
    coefs
}

##' @importFrom stats qnorm
##' @export
confint.fitC <- function (object, parm, level = 0.95, ...)
{
    ci <- NextMethod("confint")
    attr(ci, "power") <- with(object$logpower,
                              exp(estimate + c(1,-1) * qnorm((1 - level)/2) * se))
    names(attr(ci, "power")) <- colnames(ci)
    ci
}

##' @importFrom stats coef
##' @export
summary.fitC <- function (object, ...)
{
    ret <- NextMethod("summary")
    #ret$power <- attr(coef(object, se = TRUE), "power")
    ret$power <- c("Estimate" = attr(coef(object), "power"),
                   attr(confint.fitC(object), "power"))
    class(ret) <- c("summary.fitC", class(ret))
    ret
}

##' @export
print.summary.fitC <- function (x, ...)
{
    NextMethod("print")
    cat(sprintf("Power-adjusted C:  %.2f (95%% CI: %.2f to %.2f)\n\n",
                x$power[1L], x$power[2L], x$power[3L]))
    invisible(x)
}

##' @importFrom stats coef
##' @importFrom utils modifyList
##' @importFrom surveillance update.hhh4
##' @export
update.fitC <- function (object, ..., subset.upper = NULL,
                         use.estimates = object$convergence, evaluate = TRUE,
                         optim.args = list(), verbose = TRUE)
{
    if (!isTRUE(evaluate))
        stop("only 'evaluate = TRUE' is implemented for the \"fitC\" method")
    if (verbose)
        cat("Update the \"hhh4\" model ...\n")
    res <- update.hhh4(object, ..., subset.upper = subset.upper,
                       use.estimates = use.estimates)
    if (verbose)
        cat("Re-\"fitC\" (original power estimate: ",
            sprintf("%.2f", attr(coef(object), "power")), ") ...\n", sep = "")
    optim.args <- modifyList(
        c(list(control = list(trace = as.integer(verbose))),
          if (use.estimates) list(par = object$logpower$estimate)),
        optim.args)
    res <- fitC(res, # use original C specification
                C = object$logpower$args$C,
                normalize = object$logpower$args$normalize,
                truncate = object$logpower$args$truncate,
                optim.args = optim.args)
    if (verbose)
        cat("=> New estimate: ",
            sprintf("%.2f", attr(coef(res), "power")), "\n", sep = "")
    return(res)
}
