#' Chen distribption
#'
#' Defines quantile function (ICDF), and quantile density function as well as probability function (CDF) for Chen distribution
#' Distribution range is \eqn{0<x<\infty}
#' Chen distribution QF, QDF, CDF and PDF are
#' \deqn{Q(u)=\left[\ln\left(1-\frac{\ln(1-u)}{\lambda}\right)\right]^{1/\beta}}
#' \deqn{q(u)=\frac{\ln\left(1-\frac{\ln(1-u)}{\lambda}\right)^{\frac{1}{\beta}-1}}{\beta(\ln(1-u)-\lambda)(u-1)}}
#' \deqn{F(x)=1-\exp(\lambda(1-\exp(x^\beta)))}
#' \deqn{f(x)=\lambda\beta x^{\beta-1}\exp[\lambda(1-\exp(x^\beta)+x^\beta)]}
#' @param p vector of probabilities
#' @param lambda shape parameter of Chen distribution, must be positive
#' @param beta shape parameter of Chen distribution, must be positive
#' @return vector
#' @references Chen Z (2000) A new two-parameter lifetime distribution with bathtub shape or increasing failure rate function, Statistics & Probability Letters, Volume 49, Issue 2, https://doi.org/10.1016/S0167-7152(00)00044-4.
#' @name chen
#'
#' @examples
#' qchen(0.1, 0.5, 1)
#' p <- runif(1e4)
#' x <- qchen(p, 0.5, 1)+qchen(p, 1, 0.5)
#' hist(x,30)
#' @export
qchen <- function(p, lambda, beta){
  stopifnot(lambda>0, beta>0)
  (log1p(-log1p(-p)/lambda))^(1/beta)
}

#' @rdname chen
#' @export
fchen <- function(p, lambda, beta){
  stopifnot(lambda>0, beta>0)
  num <- (log1p(-log1p(-p)/lambda))^(1/beta-1)
  den <- beta*(log1p(-p)-lambda)*(p-1)
}

#' @param log logical; if TRUE, log density is returnes. Default is FALSE
#' @rdname chen
#' @export
dqchen <- function(p, lambda, beta, log=FALSE){
  res <- fchen(p, lambda=lambda, beta=beta)
  if(log) return(ifelse(is.finite(res),-log(res),res))
  1/res
}

#' @param x numeric observations
#' @rdname chen
#' @export
pchen <- function(x, lambda, beta){
  stopifnot(lambda>0, beta>0)
  -expm1(lambda*(-expm1(x^beta)))
}

#' @rdname chen
#' @export
fchen <- function(x, lambda, beta, log=FALSE){
  res <- lambda*beta*x^(beta-1)*exp(lambda*(-expm1(x^beta)+x^beta))
  if (log) return(log(res))
  res
}

#' @param n numeric; number of samples to draw from Tukey Lambda distribution
#' @rdname chen
#' @export
rchen <- function(n, lambda, beta){
  qchen(runif(n), lambda, beta)
}
