#' Generalized Normal distribution with skew parameter.
#'
#' @param x vector of quantiles. If `length(n)>1``, the length is taken to be the number required.
#' @param q1 minimum value
#' @param q3 maximum value
#' @param s skew parameter betweem -1 and 1. When skew parameter equal 0, distribution is normal.
#'
#' @return a vector of probabilities of length equal to `length(x)`.
#' @export
#' @rdname SkewNorm
#' @importFrom stats qnorm dnorm pnorm rnorm runif
#'
#' @examples
#' dSkewNorm(12, 10, 20, 0.5)
dSkewNorm <- function(x, q1, q3, s=0){
  # probability density function
  stopifnot(s>-1, s<1)
  q2 <- q1 + (s + 1) * (q3 - q1) / 2
  dMyerson(x, q1, q2, q3, alpha=stats::pnorm(-3))
}

#' @return a vector of exceedance probabilities of length equal to `length(x)`.
#' @export
#' @rdname SkewNorm
#' @examples
#' pSkewNorm(12, 10, 20, 0.5)
pSkewNorm <- function(x,q1, q3,s=0){
  # cumulative density
  stopifnot(s>-1, s<1)
  q2 <- q1 + (s + 1) * (q3 - q1) / 2
  pMyerson(x, q1, q2, q3, alpha=stats::pnorm(-3))
}

#' @param p numeric vector of probabilities
#' @return a vector of exceedance probabilities of length equal to `length(x)`.
#' @export
#' @rdname SkewNorm
#' @examples
#' qSkewNorm(12, 10, 20, 0.5)
qSkewNorm <- function(p,q1, q3,s=0){
  # cumulative density
  stopifnot(s>-1, s<1)
  q2 <- q1 + (s + 1) * (q3 - q1) / 2
  qMyerson(p, q1, q2, q3, alpha=stats::pnorm(-3))
}

#' @param n number of observations. If `length(n)>1``, the length is taken to be the number required.
#'
#' @return a length `n` vector of random values.
#' @export
#' @rdname SkewNorm
#' @examples
#' rSkewNorm(1, 10, 20, 0.5)
rSkewNorm <- function(n,q1,q3,s=0){
  # random number generator
  stopifnot(s>-1, s<1)
  q2 <- q1 + (s + 1) * (q3 - q1) / 2
  rMyerson(n, q1, q2, q3, alpha=pnorm(-3))
}


