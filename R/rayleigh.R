#' Two-parameter Rayleigh distribution
#'
#' @param x numeric vector of data
#' @param mu numeric location parameter of two-parameter Rayleigh distribution
#' @param lambda numeric scale parameter of two-parameter Rayleigh distribution
#' @rdname rayleigh
#' @return density, probability, quantile value or random variable from  two-parameter Rayleigh distribution
#' @export
#'
#' @examples
#' qrayleigh2(0.5, 0.5, 0.5)
drayleigh2 <- function(x, mu, lambda){
  xmm <- x-mu
  2*lambda*xmm*exp(-lambda*xmm^2)
}

#' @rdname rayleigh
#' @export
prayleigh2 <- function(x, mu, lambda){
  xmm <- x-mu
  1-exp(-lambda*xmm^2)
}

#' @param p numeric vector of probabilities
#' @rdname rayleigh
#' @export
qrayleigh2 <- function(p, mu, lambda){
  sqrt(-log1p(-p)/lambda)+mu
}

#' @param n numeric number of random values to draw
#' @rdname rayleigh
#' @export
rrayleigh2 <- function(n, mu, lambda){
  qrayleigh2(runif(n), mu, lambda)
}

#' @rdname rayleigh
#' @export
frayleigh2 <- function(p, mu, lambda){
  1/dqrayleigh2(p,mu,lambda)
}

#' @rdname rayleigh
#' @export
dqrayleigh2 <- function(p, mu, lambda){
  2*lambda*sqrt(-log1p(-p)/lambda)*(1-p)
}

#' @rdname rayleigh
#' @export
ffrayleigh2 <- function(p, mu, lambda){
  (2*log1p(-p)+1)/4*lambda^2*(-log1p(-p)/lambda)^(3/2)*(p-1)^2
}
