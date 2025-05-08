#' Maximum Likelihood Estimation of Contaminated Beta-Binomial Regression
#'
#' Maximum likelihood estimation of the contaminated mean-parameterized beta-binomial regression model via direct optimization.
#'
#' @param formula An object of class 'formula': a symbolic description of the model for the mean parameter (mu).
#' @param sigma.formula A formula for the dispersion parameter (sigma). Defaults to ~1.
#' @param delta.formula A formula for the proportion of extreme values parameter (delta). Defaults to ~1.
#' @param eta.formula A formula for the inflation parameter (eta). Defaults to ~1.
#' @param data A mandatory data frame containing the variables in the model.
#' @param start Vector of initial values. If NULL, values are initialized using ml.mbb and with delta = 0.01 and eta = 1.1.
#' @param verbose Logical; if TRUE, prints progress information. Default is FALSE.
#' @param method Optimization method to be used. Default is "BFGS". Other options are "Nelder-Mead".
#' @param reltol Relative convergence tolerance in optimization. Defaults to 1e-15.
#' @param maxit The maximum number of iterations in optimization. Defaults to 1000.
#' @param hessian Logical; if TRUE, computes the Hessian matrix for standard errors. Default is TRUE.
#'
#' @details The \code{ml.cmbb} function fits a contaminated beta-binomial regression model
#'
#' @return A list of elements:
#'    \item{results}{A data frame with parameter estimates, standard errors, t-values, and p-values (if hessian = TRUE).}
#'    \item{beta}{Maximum likelihood estimates of the regression coefficients for mu.}
#'    \item{alpha}{Maximum likelihood estimates of the regression coefficients for sigma.}
#'    \item{gamma}{Maximum likelihood estimates of the regression coefficients for delta.}
#'    \item{lambda}{Maximum likelihood estimates of the regression coefficients for eta.}
#'    \item{mu}{Fitted values of the mean parameter.}
#'    \item{sigma}{Fitted values of the dispersion parameter.}
#'    \item{delta}{Fitted values of the contamination proportion.}
#'    \item{eta}{Fitted values of the inflation parameter.}
#'    \item{X}{The design matrix for mu.}
#'    \item{U}{The design matrix for sigma.}
#'    \item{V}{The design matrix for delta.}
#'    \item{Z}{The design matrix for eta.}
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
#'
ml.cmbb <- function(formula,sigma.formula=~1,delta.formula=~1,eta.formula=~1, data, start = NULL, method = "BFGS", reltol=1e-15, maxit=1000, hessian=T) {
  mf <- model.frame(formula, data)
  mt <- attr(mf, "terms")
  y <- model.response(mf, "numeric")
  betaX <- model.matrix(formula, data = data)  # design matrix for mean
  sigmaU <- model.matrix(sigma.formula, data = data)  # design matrix for sigma
  deltaV <- model.matrix(delta.formula, data = data)  # design matrix for delta
  etaZ <- model.matrix(eta.formula, data = data)
  cmbb.reg.ml <- function(b.hat, X, U, V, Z, y) {
    beta.hat <- b.hat[1:ncol(X)] #coefficients for the mean (mu)
    alpha.hat <- b.hat[(ncol(X)+1):(ncol(X)+ncol(U))] #coefficients for sigma
    gamma.hat <-b.hat[(ncol(X)+ncol(U)+1):(ncol(X)+ncol(U)+ncol(V))]
    lambda.hat <- b.hat[(ncol(X)+ncol(U)+ncol(V)+1):(ncol(X)+ncol(U)+ncol(V)+ncol(Z))]

    xb.hat <- X %*% beta.hat  # mean regression
    mu.hat <- exp(xb.hat)/(1+exp(xb.hat))

    ua.hat <- U %*% alpha.hat  # sigma regression
    sigma.hat <- exp(ua.hat)

    vg.hat <- V %*% gamma.hat  # delta regression
    delta.hat <- exp(vg.hat)/(1+exp(vg.hat))

    zl.hat <- Z %*% lambda.hat  # eta regression
    eta.hat <- exp(zl.hat)+1

    ll <- sum(dcbetabinom(x = y[,1],size = rowSums(y), mu = mu.hat, sigma = sigma.hat,delta = delta.hat, eta=eta.hat, log=T))
    return(ll)
  }

  if (is.null(start)) {
    # if initial param = NULL, use param est from bb
    fm <- ml.mbb(formula = formula, sigma.formula = sigma.formula, data = data, reltol = 1e-8, method = "Nelder-Mead")
    fm$loglike
    beta.hat <- fm$beta
    alpha.hat <- fm$alpha
    gamma.hat= rep(log(0.01/(1-0.01)),ncol(deltaV))
    lambda.hat= rep(log(1.1-1),ncol(etaZ))

    start <- c(beta.hat, alpha.hat,gamma.hat,lambda.hat)

  }


  fit <- optim(par = start,
               fn = cmbb.reg.ml,
               X = betaX,
               U = sigmaU,
               V = deltaV,
               Z = etaZ,
               y = y,
               method  = method,
               control = list(fnscale = -1, maxit = maxit, reltol = reltol),
               hessian = hessian
  )


  beta.hat <- fit$par[1:ncol(betaX)]  # coefficients for the mean
  alpha.hat <- fit$par[(ncol(betaX)+1):(ncol(betaX)+ncol(sigmaU))]
  gamma.hat <- fit$par[(ncol(betaX)+ncol(sigmaU)+1):(ncol(betaX)+ncol(sigmaU)+ncol(deltaV))]
  lambda.hat <- fit$par[(ncol(betaX)+ncol(sigmaU)+ncol(deltaV)+1):(ncol(betaX)+ncol(sigmaU)+ncol(deltaV)+ncol(etaZ))]

  xb.hat <- betaX %*% beta.hat #coefficients
  mu.hat <- exp(xb.hat)/(1+exp(xb.hat))

  ua.hat <- sigmaU %*% alpha.hat
  sigma.hat <- exp(ua.hat)

  vg.hat <- deltaV %*% gamma.hat
  delta.hat <- exp(vg.hat)/(1+exp(vg.hat))

  zl.hat <- etaZ %*% lambda.hat
  eta.hat <- exp(zl.hat)+1

  lc=sum(dcbetabinom(x = y[,1],size = rowSums(y), mu = mu.hat, sigma = sigma.hat,delta = delta.hat, eta=eta.hat, log=T))

  if (hessian==T)
  {cov.mat <- base::solve(-fit$hessian)
  std.errors <- sqrt(diag(cov.mat))
  tvalue <- fit$par/std.errors
  pval <- 2*pt(abs(tvalue),df=nrow(data)-length(fit$par),lower.tail = F)
  results <- data.frame(
    Estimate = c(fit$par),
    `Std.Error` = std.errors,
    `t.value` = tvalue,
    `P.t` = round(pval,5)
  )}
  if(hessian==F){
    results <- data.frame(
      Estimate = c(fit$par))
  }
  rownames(results) <- c(paste0("beta_", colnames(betaX)), paste0("alpha_", colnames(sigmaU)), paste0("gamma_", colnames(deltaV)), paste0("lambda_", colnames(etaZ)))

  AIC <- -2*lc+nrow(results)*2
  BIC <- -2*lc+nrow(results)*log(nrow(data))
  HIC <- -2*lc +nrow(results)*2*log(log(nrow(data)))
  return(list(
    results = results,
    beta = beta.hat,
    alpha = alpha.hat,
    gamma = gamma.hat,
    lambda = lambda.hat,
    mu = mu.hat,
    sigma = sigma.hat,
    delta = delta.hat,
    eta = eta.hat,
    X = betaX,
    U = sigmaU,
    V = deltaV,
    Z = etaZ,
    y = y,
    loglike = lc,
    AIC = AIC,
    BIC = BIC,
    HIC = HIC
  ))
}
