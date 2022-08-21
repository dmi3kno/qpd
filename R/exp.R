#' Quantile derivative functions for exponential distributions
#'
#' @param u numeric probability vector
#' @param lambda parameter of exponential distribution
#'
#' @return quantile density, density quantile or quantile convexity for `fexp`, `dqexp` and `ffexp`, respectively
#' @rdname exp
#' @export
#'
#' @examples
#' fexp(0.1, 0.5)
#' dqexp(0.1, 0.5)
#' ffexp(0.1, 0.5)
fexp <- function(u, lambda){
  stopifnot(lambda>0)
  1/(lambda*(1-u))
}

#' @param log logical; if TRUE log density is returned. Default is FALSE
#' @rdname exp
#' @export
dqexp <- function(u, lambda, log=FALSE){
  stopifnot(lambda>0)
  res <- fexp(u,lambda)
  if(log) return(ifelse(is.finite(res),-log(res),res))
  1/res
}

#' @rdname exp
#' @export
ffexp <- function(u, lambda){
  stopifnot(lambda>0)
  1/(lambda*(1-u)^2)
}
