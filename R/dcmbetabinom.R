#' Contaminated Mean-Parameterized Beta-Binomial Density
#'
#' Density function of a contaminated mean-parameterized beta-binomial distribution.
#'
#' @param x Vector of non-negative integer quantities.
#' @param size Number (or vector of length x) of trials (non-negative integer).
#' @param mu Analogue of the binomial probability of success parameter, where 0 < mu < 1.
#' @param sigma Dispersion parameter, where sigma > 0.
#' @param delta Proportion of extreme values, where 0 < delta < 1.
#' @param eta Inflation parameter, where eta > 1.
#' @param log Logical; if TRUE, probabilities are returned as log(p). Default is FALSE.
#'
#' @details The contaminated beta-binomial distribution with parameters mu (\eqn{\mu}), sigma (\eqn{\sigma}), delta (\eqn{\delta}), and eta (\eqn{\eta}) has a density given by:
#' \deqn{
#' f(x; m, \mu, \sigma, \delta, \eta) = (1 - \delta) \cdot \text{BetaBinom}(x; m, \mu, \sigma) + \delta \cdot \text{BetaBinom}(x; m, \mu, \eta \sigma)
#' }
#' where \eqn{\text{BetaBinom}(x; m, \mu, \sigma)} is the density of a beta-binomial distribution with:
#' - \eqn{n} = number of trials (\code{size}),
#' Constraints: \eqn{x = 0, 1, \dots, m}, \eqn{0 < \mu < 1}, \eqn{\sigma > 0}, \eqn{0 < \delta < 1}, \eqn{\eta > 1}.
#'
#' @return The density of the contaminated beta-binomial distribution.
#'
#' @import VGAM
#'
#' @export
dcbetabinom <- function (x, size, mu, sigma, delta, eta, log = FALSE)
{
  d = (1-delta)*dmbetabinom(x,size,mu = mu, sigma = sigma,log=F)+(delta)*dmbetabinom(x,size,mu=mu,sigma = sigma*eta,log=F)

  if (log == F) {
    return(d)
  } else {
    return(log(d))
  }
}
