#' Maximum Likelihood Estimation of Contaminated Beta-Binomial Regression
#'
#' Maximum likelihood estimation of the contaminated mean-parameterized beta-binomial regression model via direct optimization.
#'
#' @param formula An object of class 'formula': a symbolic description of the model for the binomial probability of success parameter (pi).
#' @param sigma.formula A formula for the dispersion parameter (sigma). Defaults to ~1.
#' @param delta.formula A formula for the proportion of extreme values parameter (delta). Defaults to ~1.
#' @param eta.formula A formula for the inflation parameter (eta). Defaults to ~1.
#' @param data A mandatory data frame containing the variables in the model.
#' @param init Vector of initial values. If NULL, values are initialized using ml.mbb and with delta = 0.01 and eta = 1.1.
#' @param verbose Logical; if TRUE, prints progress information. Default is FALSE.
#' @param method Optimization method to be used. Default is "BFGS". Other options are "Nelder-Mead".
#' @param reltol Relative convergence tolerance in optimization. Defaults to 1e-15.
#' @param maxit The maximum number of iterations in optimization. Defaults to 10000.
#' @param hessian Logical; if TRUE, computes the Hessian matrix for standard errors. Default is TRUE.
#' @param EM Logical; if TRUE, the EM algorithm is used for maximum likelihood estimation. Default is TRUE.
#'
#'
#' @details The \code{ml.cmbb} function fits a contaminated beta-binomial regression model using maximum likelihood estimation.
#'
#' @return A list of elements:
#'    \item{results}{A data frame with parameter estimates, standard errors, t-values, and p-values (if hessian = TRUE).}
#'    \item{beta}{Maximum likelihood estimates of the regression coefficients for pi}
#'    \item{alpha}{Maximum likelihood estimates of the regression coefficients for sigma.}
#'    \item{gamma}{Maximum likelihood estimates of the regression coefficients for delta.}
#'    \item{lambda}{Maximum likelihood estimates of the regression coefficients for eta.}
#'    \item{mu}{Fitted values of the mean parameter.}
#'    \item{sigma}{Fitted values of the dispersion parameter.}
#'    \item{delta}{Fitted values of the contamination proportion.}
#'    \item{eta}{Fitted values of the inflation parameter.}
#'    \item{X}{The design matrix for pi.}
#'    \item{U}{The design matrix for sigma.}
#'    \item{V}{The design matrix for delta.}
#'    \item{Z}{The design matrix for eta.}
#'    \item{y}{The response variable.}
#'    \item{loglike}{The log-likelihood value at convergence.}
#'    \item{AIC}{Akaike Information Criterion (AIC) for the fitted model.}
#'    \item{BIC}{Bayesian Information Criterion (BIC) for the fitted model.}
#'    \item{HQIC}{Hannan-Quinn Information Criterion (HQIC) for the fitted model.}
#'
#' @import VGAM
#' @import stats
#'
#' @export
#'
ml.cmbb <- function (formula, sigma.formula = ~1, delta.formula = ~1, eta.formula = ~1,
                     data, init = NULL, method = "BFGS", reltol = 1e-15, maxit = 10000, hessian = TRUE, EM = TRUE)
{
  mf <- model.frame(formula, data)
  mt <- attr(mf, "terms")
  y <- model.response(mf, "numeric")
  X <- model.matrix(formula, data = data)
  U <- model.matrix(sigma.formula, data = data)
  V <- model.matrix(delta.formula, data = data)
  Z <- model.matrix(eta.formula, data = data)
  lc1 <- function(gamma.hat, V, w){
    vg.hat <- V %*% gamma.hat
    delta.hat <- exp(vg.hat)/(1 + exp(vg.hat))
    val <- ((1-w)*log(1-delta.hat)+w*log(delta.hat))
    return(sum(val))
  }
  lc2 <- function(par.hat, X,U,Z,y,w){
    m <- rowSums(y)
    beta.hat <- par.hat[1:ncol(X)]
    alpha.hat <- par.hat[(ncol(X) + 1):(ncol(X) + ncol(U))]
    lambda.hat <- par.hat[(ncol(X) + ncol(U) + 1):(ncol(X) + ncol(U)  + ncol(Z))]
    xb.hat <- X %*% beta.hat
    pi.hat <- exp(xb.hat)/(1 + exp(xb.hat))
    ua.hat <- U %*% alpha.hat
    sigma.hat <- exp(ua.hat)
    zl.hat <- Z %*% lambda.hat
    eta.hat <- exp(zl.hat) + 1
    val <- (1-w)*(lbeta(y[,1]+pi.hat/sigma.hat, m - y[,1]+(1-pi.hat)/sigma.hat) - lbeta(pi.hat/sigma.hat,(1-pi.hat)/sigma.hat))+  (w)*(lbeta(y[,1]+pi.hat/(eta.hat*sigma.hat), m - y[,1]+(1-pi.hat)/(eta.hat*sigma.hat)) - lbeta(pi.hat/(eta.hat*sigma.hat),(1-pi.hat)/(eta.hat*sigma.hat)))
    return(sum(val))
  }

  cmbb.reg.ml <- function(b.hat, X, U, V, Z, y) {
    beta.hat <- b.hat[1:ncol(X)]
    alpha.hat <- b.hat[(ncol(X) + 1):(ncol(X) + ncol(U))]
    gamma.hat <- b.hat[(ncol(X) + ncol(U) + 1):(ncol(X) + ncol(U) + ncol(V))]
    lambda.hat <- b.hat[(ncol(X) + ncol(U) + ncol(V) + 1):(ncol(X) + ncol(U) + ncol(V) + ncol(Z))]
    xb.hat <- X %*% beta.hat
    mu.hat <- exp(xb.hat)/(1 + exp(xb.hat))
    ua.hat <- U %*% alpha.hat
    sigma.hat <- exp(ua.hat)
    vg.hat <- V %*% gamma.hat
    delta.hat <- exp(vg.hat)/(1 + exp(vg.hat))
    zl.hat <- Z %*% lambda.hat
    eta.hat <- exp(zl.hat) + 1
    ll <- sum(dcbetabinom(x = y[, 1], size = rowSums(y),
                          mu = mu.hat, sigma = sigma.hat, delta = delta.hat,
                          eta = eta.hat, log = T))
    return(ll)
  }

  if (is.null(init)) {
    fm <- ml.mbb(formula = formula, sigma.formula = sigma.formula,
                 data = data, reltol = 1e-08, method = "Nelder-Mead")
    fm$loglike
    beta.hat <- fm$beta
    alpha.hat <- fm$alpha
    gamma.hat = rep(log(0.01/(1 - 0.01)), ncol(V))
    lambda.hat = rep(log(1.1 - 1), ncol(Z))
    start <- c(beta.hat, alpha.hat, gamma.hat, lambda.hat)
  } else {
    start <- init
  }

  if(EM==TRUE){
    ##########################################
    # EM algorithm
    ##########################################
    loglik <- -10000000000000
    converged <- FALSE
    iter <- 0
    pars <- start
    while(!converged && iter<maxit){
      iter <- iter +1
      old.loglik <- loglik
      beta   <- pars[1:ncol(X)]
      alpha  <- pars[(ncol(X)+1):(ncol(X)+ncol(U))]
      gamma  <- pars[(ncol(X)+ncol(U)+1):(ncol(X)+ncol(U)+ncol(V))]
      lambda <- pars[(ncol(X)+ncol(U)+ncol(V)+1):length(pars)]

      xb    <- X %*% beta
      pi    <- exp(xb)/(1+exp(xb))                     # mu = pi
      ua    <- U %*% alpha
      sigma <- exp(ua)
      vg    <- V %*% gamma
      delta <- exp(vg)/(1+exp(vg))
      zl    <- Z %*% lambda
      eta   <- exp(zl) + 1

      wi <- (delta)*dmbetabinom(x = y[, 1], size = rowSums(y), mu = pi, sigma = eta * sigma)/dcbetabinom(x = y[, 1], size = rowSums(y), mu = pi, sigma = sigma, delta = delta, eta = eta)



      Q1start <- gamma
      Q2start <- c(beta, alpha, lambda)

      Q1 <- optim(par = Q1start, fn = lc1, V = V, w=wi, method = method, control = list(fnscale = -1, maxit = maxit, reltol = reltol))
      Q2 <- optim(par = Q2start, fn = lc2, X = X, U = U, Z = Z, y = y, w=wi,  method = method, control = list(fnscale = -1, maxit = maxit, reltol = reltol))

      gamma.hat <- Q1$par
      beta.hat <- Q2$par[1:ncol(X)]
      alpha.hat <- Q2$par[(ncol(X) + 1):(ncol(X) + ncol(U))]
      lambda.hat <- Q2$par[(ncol(X) + ncol(U) + 1):(ncol(X) + ncol(U)  + ncol(Z))]
      xb.hat <- X %*% beta.hat
      mu.hat <- exp(xb.hat)/(1 + exp(xb.hat))
      ua.hat <- U %*% alpha.hat
      sigma.hat <- exp(ua.hat)
      vg.hat <- V %*% gamma.hat
      delta.hat <- exp(vg.hat)/(1 + exp(vg.hat))
      zl.hat <- Z %*% lambda.hat
      eta.hat <- exp(zl.hat) + 1
      loglik = sum(dcbetabinom(x = y[, 1], size = rowSums(y), mu = mu.hat,
                               sigma = sigma.hat, delta = delta.hat, eta = eta.hat,
                               log = T))
      if (abs(loglik - old.loglik) < reltol) {
        converged <- TRUE
      }
      # print(c(loglik, iter))
      pars <- c(beta.hat, alpha.hat, gamma.hat, lambda.hat)
    }
  }else{
    fit <- optim(par = start, fn = cmbb.reg.ml, X = X, U = U, V = V, Z = Z, y = y, method = method, control = list(fnscale = -1, maxit = maxit, reltol = reltol))
    beta.hat <- fit$par[1:ncol(X)]
    alpha.hat <- fit$par[(ncol(X) + 1):(ncol(X) + ncol(U))]
    gamma.hat <- fit$par[(ncol(X) + ncol(U) + 1):(ncol(X) + ncol(U) + ncol(V))]
    lambda.hat <- fit$par[(ncol(X) + ncol(U) + ncol(V) + 1):(ncol(X) + ncol(U) + ncol(V) + ncol(Z))]
    xb.hat <- X %*% beta.hat
    mu.hat <- exp(xb.hat)/(1 + exp(xb.hat))
    ua.hat <- U %*% alpha.hat
    sigma.hat <- exp(ua.hat)
    vg.hat <- V %*% gamma.hat
    delta.hat <- exp(vg.hat)/(1 + exp(vg.hat))
    zl.hat <- Z %*% lambda.hat
    eta.hat <- exp(zl.hat) + 1
    loglik = sum(dcbetabinom(x = y[, 1], size = rowSums(y), mu = mu.hat,
                             sigma = sigma.hat, delta = delta.hat, eta = eta.hat,
                             log = T))
    pars <- fit$par
  }

  if (hessian == T) {
    hess <- optimHess(pars, fn = cmbb.reg.ml, X = X, U = U, V = V, Z = Z, y = y, control = list(fnscale = -1, maxit = maxit, reltol = reltol))
    cov.mat <- base::solve(-hess)
    std.errors <- sqrt(diag(cov.mat))
    tvalue <- pars/std.errors
    pval <- 2 * pt(abs(tvalue), df = nrow(data) - length(pars),
                   lower.tail = F)
    results <- data.frame(Estimate = c(pars), Std.Error = std.errors,
                          t.value = tvalue, P.t = round(pval, 5))
  }
  if (hessian == F) {
    results <- data.frame(Estimate = c(pars))
  }
  rownames(results) <- c(paste0("beta_", colnames(X)),
                         paste0("alpha_", colnames(U)), paste0("gamma_",
                                                               colnames(V)), paste0("lambda_", colnames(Z)))
  AIC <- -2 * loglik + nrow(results) * 2
  BIC <- -2 * loglik + nrow(results) * log(nrow(data))
  HIC <- -2 * loglik + nrow(results) * 2 * log(log(nrow(data)))
  return(list(results = results, beta = beta.hat, alpha = alpha.hat,
              gamma = gamma.hat, lambda = lambda.hat, mu = mu.hat,
              sigma = sigma.hat, delta = delta.hat, eta = eta.hat,
              X = X, U = U, V = V, Z = Z, y = y, loglike = loglik,
              AIC = AIC, BIC = BIC, HQIC = HIC, EM=EM))
}
