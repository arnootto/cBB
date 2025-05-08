#' Mean-Parameterized Beta-Binomial Density
#'
#' Density function of a mean-parameterized beta-binomial distribution.
#'
#' @param x Vector of non-negative integer quantities.
#' @param size Number (or vector of length x) of trials (non-negative integer).
#' @param mu Analogue of the binomial probability of success parameter, where 0 < mu < 1.
#' @param sigma Dispersion parameter, where sigma > 0.
#' @param log Logical; if TRUE, probabilities are returned as log(p). Default is FALSE.
#'
#' @details The beta-binomial distribution with parameters mu (\eqn{\mu}) and sigma (\eqn{\sigma}) has a density given by:
#' \deqn{
#' f(x; n, \mu, \sigma) = \binom{m}{x} \frac{B(x + \frac{\mu}{\sigma}, m - x + \frac{1-\mu}{\sigma})}{B(\frac{\mu}{\sigma},\frac{1-\mu}{\sigma})}
#' }
#' where:
#' - \eqn{m} = number of trials (\code{size}),
#' - \eqn{B} is the beta function.
#' Constraints: \eqn{x = 0, 1, \dots, m}, \eqn{0 < \mu < 1}, \eqn{\sigma > 0}.
#'
#' @return The density of the beta-binomial distribution.
#'
#' @import VGAM
#'
#' @export
dmbetabinom <- function (x, size, mu, sigma, log = FALSE)
{  alpha = mu/(sigma)
beta  = (1-mu)/(sigma)
return(dbetabinom.ab(x= x, size =size, alpha, beta,log=log))
}
