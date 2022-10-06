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
