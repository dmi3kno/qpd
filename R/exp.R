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
  1/(lambda*(1-u))
}

#' @rdname exp
#' @export
dqexp <- function(u, lambda){
  lambda*(1-u)
}

#' @rdname exp
#' @export
ffexp <- function(u, lambda){
  1/(lambda*(1-u)^2)
}
