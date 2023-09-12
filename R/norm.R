#' @keywords internal
#' @description
#' The Error Function is related to CDF of normal as \deqn{\text{erf}(x)=2\Phi(x\sqrt(2))-1}
#' The quantile function is the consequence of \deqn{\Phi^{-1}(p)=\sqrt{2}\text{erf}^{-1}(2p-1)}
#' @importFrom stats qnorm dnorm
#' @importFrom Rmpfr log1pexp log1mexp
perf <- function(x, lower.tail = TRUE, log.p = FALSE){
  p <- stats::pnorm(x * sqrt(2), lower.tail = lower.tail, log.p = log.p)
  if(log.p) return(log(exp(p)+expm1(p)))
  2 * p - 1
}
# https://stats.stackexchange.com/questions/187828/how-are-the-error-function-and-standard-normal-distribution-function-related
qerf <- function(p, lower.tail = TRUE, log.p = FALSE) {
  if (log.p) p <- Rmpfr::log1pexp(p)-log(2) else p <- (1+p)/2
  stats::qnorm(p, lower.tail = lower.tail, log.p = log.p)/sqrt(2)
}
#'

#' Derivative quantile functions for normal distribution
#'
#' Calculate QDF `fnorm` and DQF `dqnorm` for the normal distribution
#' @param p numeric probability
#' @param mean numeric parameter of normal distribution. Default is 0
#' @param sd numeric parameter of normal distribution. Default is 1
#' @param log logical, should the result be returned as a log. Default is FALSE
#' @param log.p logical; if TRUE, probabilities p are given as log(p). Default is FALSE
#' @param lower.tail logical; if TRUE (default), probabilities are \eqn{P[X\leq x]}, otherwise \eqn{P[X>x]}
#'
#' @return returns QDF (`fnorm`) and DQF (`dqnorm`) of normal distribution
#' @export
#' @importFrom stats qnorm dnorm
#' @rdname norm
#' @examples
#' p_grd <- make_pgrid()
#' fnorm(p_grd)
#' dqnorm(p_grd)
fnorm <- function(p, mean=0, sd=1, log=FALSE, lower.tail=TRUE, log.p=FALSE){
  z <- stats::qnorm(p, mean, sd, lower.tail=lower.tail, log.p=log.p)
  res <- 1/stats::dnorm(z, mean, sd, log=log)
  if(log) return(log(res))
  res
}

#' @rdname norm
#' @export
#' @importFrom stats qnorm dnorm
dqnorm <- function(p, mean=0, sd=1, log=FALSE, lower.tail=TRUE, log.p=FALSE){
  z <- stats::qnorm(p, mean, sd, lower.tail=lower.tail, log.p=log.p)
  res <- stats::dnorm(z, mean, sd, log=log)
  if(log) return(log(res))
  res
}

#' Half normal distribution
#'
#' @param p numeric probability
#' @param location numeric parameter of halfnormal distribution. Default is 0
#' @param scale numeric parameter of halfnormal distribution. Default is 1
#' @param log logical, should the result be returned as a log. Default is FALSE
#' @param log.p logical; if TRUE, probabilities p are given as log(p). Default is FALSE
#' @param lower.tail logical; if TRUE (default), probabilities are \eqn{P[X\leq x]}, otherwise \eqn{P[X>x]}

#' @return returns QDF (`fnorm`) and DQF (`dqnorm`) of normal distribution
#' @export
#' @importFrom stats qnorm dnorm
#' @rdname halfnorm
#' @examples
#' p_grd <- make_pgrid()
#' qhalfnorm(p_grd)

qhalfnorm <- function(p, location=0, scale=1, lower.tail=TRUE, log.p=FALSE){
  location+scale*sqrt(2)*qerf(p, lower.tail=lower.tail, log.p=log.p)
}

#' @param x numeric value
#' @export
#' @importFrom stats qnorm dnorm
#' @rdname halfnorm
phalfnorm <- function(x, location=0, scale=1, lower.tail=TRUE, log.p=FALSE){
  perf((x-location)/(scale*sqrt(2)), lower.tail=lower.tail, log.p=log.p)
}

#' @export
#' @importFrom stats qnorm dnorm
#' @rdname halfnorm
dhalfnorm <- function(x, location=0, scale=1, log=FALSE){
  res <- ifelse(x<=location, 0, sqrt(2/pi)/scale*exp(-1/2*((x-location)/scale)^2))
  if(log) return(log(res))
  res
}
