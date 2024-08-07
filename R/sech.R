#' Hyperbolic secant distribution
#'
#' Hyperbolic secant distribution with CDF \eqn{F(x)=\frac{2}{\pi}\arctan(\exp(\frac{\pi}{2}\frac{(x-\mu)}{\sigma}))} where \eqn{\mu} and \eqn{\sigma} are location and scale, respectively.
#'
#' @param x numeric vector of values
#' @param location parameter of hyperbolic secant distribution
#' @param scale parameter of yperbolic secant distribution
#' @rdname sech
#'
#' @return cum
#'
#' @examples
#' psech(1, 1, 1)
#' @export
psech <- function(x, location=0, scale=1){
  stopifnot(scale>0)
  2/pi*atan(exp(pi/2*(x-location)/scale))
}

#' @rdname sech
#' @export
dsech <- function(x, location=0, scale=1){
  stopifnot(scale>0)
  1/2*sech(pi/2*(x-location)/scale)
}

#' @param p numeric vector of cumulative probabilities
#' @rdname sech
#' @export
qsech <- function(p, location=0, scale=1){
  stopifnot(scale>0)
  location+scale*(2/pi*log(tan(pi/2*p)))
}

#' @rdname sech
#' @export
fsech <- function(p, location=0, scale=1){
  #sec = 1 / cos x
  #sec^2(x)=1+tan^2(x)
  stopifnot(scale>0)
  hpip <- pi/2*p
  scale*(1+tan(hpip)^2)/tan(hpip)
}

#' @param log logical; if TRUE, log density is returned. Default is FALSE
#' @rdname sech
#' @export
dqsech <- function(p, location=0, scale=1, log=FALSE){
  res <- fsech(p, location, scale)
  if(log) return(ifelse(is.finite(res),-log(res),res))
  1/res
}

#' @param n number of samples to draw
#' @rdname sech
#' @export
rsech <- function(n, location=0, scale=1){
  qsech(runif(n), location, scale)
}

