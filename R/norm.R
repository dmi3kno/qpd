#' @keywords internal
#' @importFrom stats qnorm dnorm
perf <- function(x) 2 * stats::pnorm(x * sqrt(2)) - 1
qerf <- function(p) stats::qnorm((p + 1)/2)/sqrt(2)
#'

#' Derivative quantile functions for normal distribution
#'
#' Calculate QDF `fnorm` and DQF `dqnorm` for the normal distribution
#' @param p numeric probability
#' @param mean numeric parameter of normal distribution. Default is 0
#' @param sd numeric parameter of normal distribution. Default is 1
#' @param log logical, should the result be returned as a log. Default is FALSE
#' @param log.p logical; if TRUE, probabilities p are given as log(p). Default is FALSE
#'
#' @return returns QDF (`fnorm`) and DQF (`dqnorm`) of normal distribution
#' @export
#' @importFrom stats qnorm dnorm
#' @rdname norm
#' @examples
#' p_grd <- make_pgrid()
#' fnorm(p_grd)
#' dqnorm(p_grd)
fnorm <- function(p, mean=0, sd=1, log=FALSE, log.p=FALSE){
  z <- stats::qnorm(p, mean, sd, log.p=log.p)
  res <- 1/stats::dnorm(z, mean, sd)
  if(log) return(log(res))
  res
}

#' @rdname norm
#' @export
#' @importFrom stats qnorm dnorm
dqnorm <- function(p, mean=0, sd=1, log=FALSE, log.p=FALSE){
  z <- stats::qnorm(p, mean, sd, log.p=log.p)
  res <- stats::dnorm(z, mean, sd)
  if(log) return(ifelse(is.finite(res),log(res),res))
  res
}

#' Half normal distribution
#'
#' @param p numeric probability
#' @param location numeric parameter of halfnormal distribution. Default is 0
#' @param scale numeric parameter of halfnormal distribution. Default is 1
#' @param log logical, should the result be returned as a log. Default is FALSE
#' @param log.p logical; if TRUE, probabilities p are given as log(p). Default is FALSE
#'
#' @return returns QDF (`fnorm`) and DQF (`dqnorm`) of normal distribution
#' @export
#' @importFrom stats qnorm dnorm
#' @rdname halfnorm
#' @examples
#' p_grd <- make_pgrid()
#' qhalfnorm(p_grd)

qhalfnorm <- function(p, location=0, scale=1, log.p=FALSE){
  location+scale*sqrt(2)*qerf(p)
}

#' @param x numeric value
#' @export
#' @importFrom stats qnorm dnorm
#' @rdname halfnorm
phalfnorm <- function(x, location=0, scale=1, log.p=FALSE){
  perf((x-location)/(scale*sqrt(2)))
}

#' @export
#' @importFrom stats qnorm dnorm
#' @rdname halfnorm
dhalfnorm <- function(x, location=0, scale=1, log.p=FALSE){
  ifelse(x<=location, 0, sqrt(2/pi)/scale*exp(-1/2*((x-location)/scale)^2))
}
