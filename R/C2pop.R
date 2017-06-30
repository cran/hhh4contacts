################################################################################
### Adjust a contact matrix to obtain a desired stationary distribution
###
### Copyright (C) 2016 Leonhard Held, Sebastian Meyer
###
### This file is part of the R package "hhh4contacts",
### free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at http://www.r-project.org/Licenses/.
################################################################################

##' Stationary Distribution of a Transition Matrix
##'
##' This auxiliary function determines the stationary distribution
##' from a transition matrix.
##' @param P a transition matrix, i.e., a square matrix where all rows sum to 1.
##' @return the stationary probability vector.
##' @author Leonhard Held
##' @examples
##' Cgrouped_norm <- contactmatrix(normalize = TRUE)
##' Cgrouped_norm
##' (p <- stationary(Cgrouped_norm))
##' (Cpowered <- make_powerC(Cgrouped_norm)(1e6))
##' stopifnot(all.equal(Cpowered[1,], p))
##' @export
stationary <- function(P)
{
    I <- diag(nrow = nrow(P))
    colSums(solve(I-P+1))
    ## alternative calculation: prop.table(eigen(t(P))$vectors[,1])
}

##' Adapt a Transition Matrix to a Specific Stationary Distribution
##'
##' \strong{Experimental} Metropolis-Hastings algorithm, which tries
##' to adjust a transition matrix such that its stationary distribution
##' becomes approximately equal to a prespecified probability vector.
##' @param P a transition matrix, i.e., a square matrix where all rows sum to 1.
##' @param target the stationary probability vector to approximate.
##' @param niter the number of iterations of the MCMC algorithm
##' @return the adjusted transition matrix.
##' @author Leonhard Held
##' @seealso \code{\link{C2pop}} for an alternative method.
##' @examples
##' ## a row-normalized contact matrix
##' C <- matrix(c(0.8, 0.1, 0.1,
##'               0.2, 0.6, 0.2,
##'               0.1, 0.2, 0.7), byrow=TRUE, ncol=3, nrow=3)
##' stationary(C)
##' ## population fractions define the target distribution
##' popfracs <- c(0.4, 0.3, 0.3)
##' ## adapt 'C' to the given population fractions
##' Cpop <- adaptP(C, popfracs, niter = 50000)
##' stationary(Cpop)
##' ## this method increases the diagonal values of 'C'
##' round(C, 3)
##' round(Cpop, 3)
##' round(Cpop/C, 3)
##' @importFrom stats runif
##' @export
adaptP <- function(P, target, niter = 1000000)
{
    stopifnot(is.vector(target, mode = "numeric"),
              (dim <- nrow(P)) == length(target))
    sums <- matrix(0, nrow=dim, ncol=dim, dimnames=dimnames(P))
    yes <- 0
    no <- 0
    yesvec <- rep.int(0, dim)
    novec <- rep.int(0, dim)
    old <- new <- sample.int(dim, size=1, prob=target)
    for (i in 2:niter) {
        old <- new
        new <- sample.int(dim, size=1, prob=P[old,])
        acc <- (target[new]*P[new,old])/(target[old]*P[old,new])
        if (runif(1) <= acc) {
            yes <- yes + 1
            yesvec[old] <- yesvec[old] + 1
        } else {
            new <- old
            no <- no + 1
            novec[old] <- novec[old] + 1
        }
        sums[old,new] <- sums[old,new] + 1
    }
    res <- sums/.rowSums(sums, dim, dim)
    attr(res, "acc") <- yes/(yes+no)
    attr(res, "accvec") <- yesvec/(yesvec+novec)
    return(res)
}

##' Adapt a Contact Matrix to Population Fractions
##'
##' \strong{Experimental} function, which tries to adjust a given contact matrix
##' such that the stationary distribution of its row-normalized version (i.e.,
##' the transition matrix) becomes approximately equal to a prespecified
##' probability vector.
##' @param C a square numeric (contact) matrix.
##' @param target the stationary probability vector to approximate.
##' @param eps the tolerated mean absolute difference between the target
##'     probabilities and the stationary distribution of the adapted, normalized
##'     contact matrix.
##' @param iter.max maximum number of iterations (guard against infinite loop).
##' @return the adapted, normalized contact matrix.
##' @author Leonhard Held (original) and Sebastian Meyer (this implementation)
##' @seealso \code{\link{adaptP}} for an alternative method.
##' @examples
##' GROUPING <- c(1, 2, 2, 4, 4, 2)
##' C <- contactmatrix(grouping = GROUPING)
##' popBErbyg <- aggregateCountsArray(pop2011, dim = 2, grouping = GROUPING)
##' popfracs <- prop.table(colSums(popBErbyg))
##' ## adapt 'C' to the given population fractions
##' Cpop <- C2pop(C, popfracs)
##' ## compare the stationary distributions
##' compstat <- cbind(before = stationary(C/rowSums(C)), popBE = popfracs,
##'                   after = stationary(Cpop))
##' matplot(compstat, type="b", lty=1, ylim=c(0, max(compstat)),
##'         xlab="age group", ylab="population fraction")
##' ## compare the normalized contact matrices
##' print(plotC(C/rowSums(C), main="original", at=seq(0,0.6,length.out=17)),
##'       split=c(1,1,2,1), more=TRUE)
##' print(plotC(Cpop, main="adapted", at=seq(0,0.6,length.out=17)),
##'       split=c(2,1,2,1), more=FALSE)
##' @export
C2pop <- function (C, target, eps = 0.001, iter.max = 100)
{
    stopifnot(is.matrix(C), is.numeric(C), diff.default(dimC <- dim(C)) == 0)
    stopifnot(is.vector(target, mode = "numeric"),
              length(target) == (dim <- dimC[1L]))
    close_enough <- function (current) mean(abs(target - current)) <= eps
    normalize <- function (C) C/.rowSums(C, dim, dim)
    i <- 0
    while(!close_enough(current <- stationary(Cnorm <- normalize(C)))) {
        if (i == iter.max) {
            warning("reached maximum number of iterations")
            break
        }
        C <- C * matrix(target/current, byrow=TRUE, nrow=dim, ncol=dim)
        i <- i + 1
    }
    return(Cnorm)
}
