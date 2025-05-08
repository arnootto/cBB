#' Simulate Beta-Binomial Regression Data
#'
#' Simulates data from a mean-parameterized beta-binomial regression model.
#'
#' @param n Number of observations to simulate. Default is 100.
#' @param size Number of trials for the beta-binomial distribution. Default is 10.
#' @param beta Vector of regression coefficients for the mean parameter (mu). Default is c(1, 0.75, -1.5).
#' @param alpha Vector of regression coefficients for the dispersion parameter (sigma). Default is c(2, 5).
#'
#' @details The \code{rbb_regression} function simulates data from a beta-binomial regression model where the response follows a beta-binomial distribution with mean \eqn{\mu} and dispersion \eqn{\sigma}.]
#' The predictors \eqn{X} and \eqn{U} are generated as uniform random variables in [-2, 2].
#'
#' @return A data frame containing:
#'    \item{bby}{The simulated beta-binomial responses.}
#'    \item{x0, x1, ...,xm}{The design matrix columns for mu (including intercept).}
#'    \item{u0, u1, ...,um}{The design matrix columns for sigma (including intercept).}
#'
#' @import VGAM
#'
#' @export
rmbb_regression <- function(n = 100, size = 10, beta = c(1, 0.75, -1.5), alpha = c(2, 5)) {
  dX <- length(beta) - 1
  dU <- length(alpha) - 1

  X <- matrix(runif(n * dX, -2, 2), nrow = n, ncol = dX)
  Xstar <- cbind(rep(1, n), X)

  U <- matrix(runif(n * dU, -2, 2), nrow = n, ncol = dU)
  Ustar <- cbind(rep(1, n), U)

  lin.pred <- Xstar %*% beta
  cond.mu <- exp(lin.pred) / (1 + exp(lin.pred))

  lin.pred.sigma <- Ustar %*% alpha
  cond.sigma <- exp(lin.pred.sigma)  # log link for sigma

  nby <- VGAM::rbetabinom.ab(n, size, shape1 = cond.mu / cond.sigma, shape2 = (1 - cond.mu) / cond.sigma)
  out <- data.frame(cbind(nby, Xstar, Ustar))
  names(out) <- c("bby", paste("x", 0:dX, sep = ""), paste("u", 0:dU, sep = ""))
  return(out)
}
