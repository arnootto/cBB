#' Simulate Contaminated Beta-Binomial Regression Data
#'
#' Simulates data from a contaminated mean-parameterized beta-binomial regression model.
#'
#' @param n Number of observations to simulate. Default is 1000.
#' @param size Number of trials for the beta-binomial distribution. Default is 10.
#' @param beta Vector of regression coefficients for the mean parameter (mu). Default is c(1, 0.75).
#' @param alpha Vector of regression coefficients for the dispersion parameter (sigma). Default is c(log(0.5)).
#' @param gamma Vector of regression coefficients for the contamination proportion (delta). Default is c(log(0.05 / (1 - 0.05))).
#' @param lambda Vector of regression coefficients for the inflation parameter (eta). Default is c(log(1.2 - 1)).
#'
#' @details The \code{rcbb_regression} function simulates data from a contaminated beta-binomial regression model
#' The predictors \eqn{X}, \eqn{U}, \eqn{V}, and \eqn{Z} are generated as uniform random variables in [-2, 2].
#'
#' @return A data frame containing:
#'    \item{cbby}{The simulated contaminated beta-binomial responses.}
#'    \item{x0, x1, ...}{The design matrix columns for mu (including intercept).}
#'    \item{u0, u1, ...}{The design matrix columns for sigma (including intercept).}
#'    \item{v0, v1, ...}{The design matrix columns for delta (including intercept).}
#'    \item{z0, z1, ...}{The design matrix columns for eta (including intercept).}
#'
#' @import VGAM
#'
#' @export
rcmbb_regression <- function(n = 1000, size = 10, beta = c(1, 0.75), alpha = c(log(0.5)),
                            gamma = c(log(0.05 / (1 - 0.05))), lambda = c(log(1.2 - 1))) {
  dX <- length(beta) - 1  # Number of predictors for mu
  dU <- length(alpha) - 1  # Number of predictors for sigma
  dV <- length(gamma) - 1  # Number of predictors for delta
  dZ <- length(lambda) - 1  # Number of predictors for eta

  X <- matrix(runif(n * dX, -2, 2), nrow = n, ncol = dX)
  Xstar <- cbind(1, X)  # Include intercept

  U <- matrix(runif(n * dU, 2, -2), nrow = n, ncol = dU)
  Ustar <- cbind(1, U)  # Include intercept

  V <- matrix(runif(n * dV, -2, 2), nrow = n, ncol = dV)
  Vstar <- cbind(1, V)  # Include intercept

  Z <- matrix(runif(n * dZ, -2, 2), nrow = n, ncol = dZ)
  Zstar <- cbind(1, Z)  # Include intercept

  lin.pred.mu <- Xstar %*% beta
  cond.mu <- exp(lin.pred.mu) / (1 + exp(lin.pred.mu))

  lin.pred.sigma <- Ustar %*% alpha
  cond.sigma <- exp(lin.pred.sigma)  # log link for sigma

  lin.pred.delta <- Vstar %*% gamma
  cond.delta <- exp(lin.pred.delta) / (1 + exp(lin.pred.delta))  # logistic link for delta

  lin.pred.eta <- Zstar %*% lambda
  cond.eta <- exp(lin.pred.eta) + 1  # shifted log link for eta

  rbinom <- rbinom(n, size = 1, prob = cond.delta)

  rbad <- VGAM::rbetabinom.ab(n = sum(rbinom == 1), size = size,
                              shape1 = cond.mu[rbinom == 1] / (cond.eta[rbinom == 1] * cond.sigma[rbinom == 1]),
                              shape2 = (1 - cond.mu[rbinom == 1]) / (cond.eta[rbinom == 1] * cond.sigma[rbinom == 1]))

  rgood <- VGAM::rbetabinom.ab(n = sum(rbinom == 0), size = size,
                               shape1 = cond.mu[rbinom == 0] / cond.sigma[rbinom == 0],
                               shape2 = (1 - cond.mu[rbinom == 0]) / cond.sigma[rbinom == 0])

  rcBB <- numeric(n)
  rcBB[rbinom == 1] <- rbad
  rcBB[rbinom == 0] <- rgood

  cbby <- cbind(rcBB, Xstar, Ustar, Vstar, Zstar)

  colnames(cbby) <- c("cbby", paste("x", 0:dX, sep = ""), paste("u", 0:dU, sep = ""),
                      paste("v", 0:dV, sep = ""), paste("z", 0:dZ, sep = ""))
  return(cbby)
}
