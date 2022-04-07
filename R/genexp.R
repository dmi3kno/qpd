

#' Generalized exponential distribution
#'
#' Generalized exponential distribution with CDF F(x)=(1-exp(-lambda*x))^alpha
#'
#' @param x numeric vector of values
#' @param alpha parameter of generalized exponential distribution
#' @param lambda parameter of generalized exponential distribution
#' @rdname genexp
#'
#' @return cum
#'
#' @examples
#' pgenexp(1, 1, 1)
#' @export
pgenexp <- function(x, alpha, lambda){
  stopifnot(alpha>0)
  stopifnot(lambda>0)
  (1-exp(-lambda*x))^alpha
}

#' @rdname genexp
#' @export
dgenexp <- function(x, alpha, lambda){
  stopifnot(alpha>0)
  stopifnot(lambda>0)
  alpha*lambda*(1-exp(-lambda*x))^(alpha-1)*exp(-lambda*x)
}

#' @param p numeric vector of cumulative probabilities
#' @rdname genexp
#' @export
qgenexp <- function(p, alpha, lambda){
  stopifnot(alpha>0)
  stopifnot(lambda>0)
  inv_lambda <- 1/lambda
  inv_alpha <- 1/alpha
  -(inv_lambda)*log(1-p^inv_alpha)
}

#' @rdname genexp
#' @export
fgenexp <- function(p, alpha, lambda){
  stopifnot(alpha>0)
  stopifnot(lambda>0)
  inv_lambda <- 1/lambda
  inv_alpha <- 1/alpha
  p^(inv_alpha-1)/(alpha*lambda*(1-p^inv_alpha))
}

#' @param log logical; if TRUE, log density is returned. Default is FALSE
#' @rdname genexp
#' @export
dqgenexp <- function(p, alpha, lambda, log=FALSE){
  res <- fgenexp(p, alpha, lambda)
  if (log) return(log(1/res))
  1/res
}

#' @param n number of samples to draw
#' @rdname genexp
#' @export
rgenexp <- function(n, alpha, lambda){
  qgenexp(runif(n), alpha, lambda)
}

