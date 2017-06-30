################################################################################
### Construct a function for matrix exponentiation via eigendecomposition
###
### Copyright (C) 2015-2017 Sebastian Meyer
###
### This file is part of the R package "hhh4contacts",
### free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at http://www.r-project.org/Licenses/.
################################################################################

##' Exponentiate a Matrix via Eigendecomposition
##'
##' Based on a (contact) matrix \code{C}, the function \code{make_powerC}
##' generates a function with a single argument \code{power} that returns
##' the input matrix raised to that power. Matrix exponentiation is thereby
##' defined via the eigendecomposition of \code{C} as
##' \eqn{C^{power} := E \Lambda^{power} E^{-1}}.
##' @param C a square numeric matrix.
##' @param normalize a logical indicating if \code{C} should be normalized
##' in advance such that all rows sum to 1 (becomes a transition matrix).
##' @param truncate a logical indicating whether to force entries in the
##' resulting matrix to be non-negative (by truncation at 0).
##' @return a function of the \code{power}
##' that returns the exponentiated matrix.
##' @examples
##' Cnorm <- contactmatrix(normalize = TRUE)
##' powerC <- make_powerC(Cnorm)
##' powerC(1)
##' powerC(0)
##' powers <- c(0, 0.5, 1, 2)
##' Cp <- lapply(powers, powerC)
##' if (require("gridExtra"))
##'     grid.arrange(
##'         grobs = mapply(plotC, C = Cp, main = paste("power =", powers),
##'                        SIMPLIFY = FALSE),
##'         nrow = 2, ncol = 2)
##'
##' ## truncation to enforce non-negative entries
##' powerC(0.2)  # some entries become negative for small powers
##' powerC0 <- make_powerC(Cnorm, truncate = TRUE)
##' powerC0(0.2)
##' @name powerC
##' @export
make_powerC <- function (C, normalize = FALSE, truncate = FALSE)
{
    EV <- eigen(if (normalize) C / rowSums(C) else C, symmetric = FALSE)
    if (complexEV <- is.complex(EV$values) || is.complex(EV$vectors))
        warning(if (normalize) "normalized ",
                "'C' has complex eigen values or vectors",
                immediate. = TRUE)

    ## template function (assuming complexEV and no truncation)
    powerC <- function (power = 1) {
        C[] <- EV$vectors %*% diag(EV$values^power) %*% solve(EV$vectors)
        if (all(Im(zapsmall(C)) == 0)) { # zap zero imaginary parts
            C <- Re(C)
        } else warning("powered C is complex")
        C
    }

    if (!complexEV) { # no need to check for zero imaginary parts
        body(powerC)[[3]] <- NULL
    }

    if (truncate)
        body(powerC)[[length(body(powerC))]] <-
            if (complexEV) { # truncation needs a real-valued result
                quote(if (is.complex(C)) C else pmax(C, 0))
            } else {
                quote(pmax(C, 0))
            }
    powerC
}
