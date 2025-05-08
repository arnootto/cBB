#' Maximum Likelihood Estimation of Beta-Binomial Regression
#'
#' Maximum likelihood estimation of the mean-parameterized beta-binomial regression model.
#'
#' @param formula An object of class 'formula': a symbolic description of the model to be fitted.
#' @param sigma.formula A formula for the dispersion parameter (sigma). Defaults to ~1.
#' @param data A mandatory data frame containing the variables in the model.
#' @param start Vector of initial values. If NULL, values produced by glm(family=Binomial) are used as initial values.
#' @param method Optimization method to be used. Default is "BFGS". Other options are "Nelder-Mead".
#' @param reltol Relative convergence tolerance in each optimization process. Defaults to 1e-15.
#' @param maxit The number of inner iterations in each optimization process. Defaults to 1000.
#' @param hessian Logical; if TRUE, computes the Hessian matrix for standard errors. Default is TRUE.
#'
#' @details The \code{ml.mbb} function fits the mean-parameterized beta-binomial regression model.
#'
#' @return A list of elements:
#'    \item{results}{A data frame with parameter estimates and standard errors.}
#'    \item{alpha}{Maximum likelihood estimate of the coefficients for sigma.}
#'    \item{beta}{Maximum likelihood estimates of the regression coefficients for mu.}
#'    \item{mu}{Fitted values of the mean parameter.}
#'    \item{sigma}{Fitted values of the dispersion parameter.}
#'    \item{X}{The design matrix for mu.}
#'    \item{U}{The design matrix for sigma.}
#'    \item{y}{The response variable.}
#'    \item{loglike}{The log-likelihood value at convergence.}
#'    \item{AIC}{Akaike Information Criterion (AIC) for the fitted model.}
#'    \item{BIC}{Bayesian Information Criterion (BIC) for the fitted model.}
#'    \item{HIC}{Hannan-Quinn Information Criterion (HIC) for the fitted model.}
#'
#' @import VGAM
#' @import stats
#'
#' @export
ml.mbb <- function(formula,sigma.formula =~1, data, start = NULL, method = "BFGS", reltol=1e-15, maxit=1000, hessian=T) {
  mf <- model.frame(formula, data)
  msigma <- model.frame(sigma.formula,data)
  mt <- attr(mf, "terms")
  y <- model.response(mf, "numeric")
  bbX <- model.matrix(formula, data = data)  #design matrix for mean
  sigmaU <- model.matrix(sigma.formula, data = data)  #design matrix for sigma
  mbb.reg.ml <- function(b.hat, X, U, y) {
    beta.hat <- b.hat[1:ncol(X)] #coefficients for the mean (mu)
    alpha.hat <- b.hat[(ncol(X)+1):(ncol(X)+ncol(U))] #coefficients for sigma

    xb.hat <- X %*% beta.hat  # mean regression
    mu.hat <- exp(xb.hat) / (1 + exp(xb.hat))


    ua.hat <- U %*% alpha.hat  # sigma regression
    sigma.hat <- exp(ua.hat)

    ll <- sum(dmbetabinom(x = y[,1],size = rowSums(y), mu = mu.hat, sigma = sigma.hat, log=T))
    return(ll)
  }

  if (is.null(start)) {
    fm <- glm(formula = formula, family = binomial, data = data)
    beta.hat <-coefficients(fm)
    xb.start= bbX %*% beta.hat
    sigma.hat= rep(log(1.1),ncol(sigmaU)) #link
    start <- c(beta.hat, sigma.hat)
  }

  fit <- optim(par = start,
               fn  = mbb.reg.ml,
               X   = bbX,
               U   = sigmaU,
               y   = y,
               method  = method,
               control = list(fnscale = -1, maxit = maxit, reltol = reltol),
               hessian = hessian
  )

  beta.hat <- fit$par[1:ncol(bbX)]  # coefficients for the mean
  alpha.hat <- fit$par[(ncol(bbX)+1):(ncol(bbX)+ncol(sigmaU))]

  xb.hat <- bbX %*% beta.hat #coefficients
  mu.hat <- exp(xb.hat)/(1+exp(xb.hat))

  ua.hat <- sigmaU %*% alpha.hat
  sigma.hat <- exp(ua.hat)

  lc=sum(dmbetabinom(x = y[,1],size = rowSums(y), mu = mu.hat,sigma = sigma.hat,log=T))

  if (hessian==T)
  { cov.mat <- base::solve(-fit$hessian)
  std.errors <- sqrt(diag(cov.mat))
  tvalue <- fit$par/std.errors
  pval <- 2*pt(abs(tvalue),df=nrow(data)-length(fit$par),lower.tail = F)
  results <- data.frame(
    Estimate = c(fit$par),
    `Std.Error` = std.errors,
    `t.value` = tvalue,
    `P.t` = pval
  )}
  if (hessian==F){
    results <- data.frame(
      Estimate = c(fit$par))
  }
  rownames(results) <- c(paste0("beta_", colnames(bbX)), paste0("alpha_", colnames(sigmaU)))

  AIC <- -2*lc+nrow(results)*2
  BIC <- -2*lc+nrow(results)*log(nrow(data))
  HIC <- -2*lc +nrow(results)*2*log(log(nrow(data)))
  return(list(
    results = results,
    beta = beta.hat,
    alpha = alpha.hat,
    mu = mu.hat,
    sigma = sigma.hat,
    X = bbX,
    U = sigmaU,
    y = y,
    loglike = lc,
    AIC = AIC,
    BIC = BIC,
    HIC = HIC
  ))
}
