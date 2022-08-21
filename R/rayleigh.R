#' Two-parameter Rayleigh distribution
#'
#' @param x numeric vector of data
#' @param mu numeric location parameter of two-parameter Rayleigh distribution
#' @param lambda numeric scale parameter of two-parameter Rayleigh distribution
#' @rdname rayleigh2
#' @return density, probability, quantile value or random variable from  two-parameter Rayleigh distribution
#' @export
#'
#' @examples
#' qrayleigh2(0.5, 0.5, 0.5)
drayleigh2 <- function(x, mu, lambda){
  xmm <- x-mu
  2*lambda*xmm*exp(-lambda*xmm^2)
}

#' @rdname rayleigh2
#' @export
prayleigh2 <- function(x, mu, lambda){
  xmm <- x-mu
  1-exp(-lambda*xmm^2)
}

#' @param p numeric vector of probabilities
#' @rdname rayleigh2
#' @export
qrayleigh2 <- function(p, mu, lambda){
  sqrt(-log1p(-p)/lambda)+mu
}

#' @param n numeric number of random values to draw
#' @rdname rayleigh2
#' @export
rrayleigh2 <- function(n, mu, lambda){
  qrayleigh2(runif(n), mu, lambda)
}

#' @rdname rayleigh2
#' @export
frayleigh2 <- function(p, mu, lambda){
  1/dqrayleigh2(p,mu,lambda, log=FALSE)
}

#' @param log logical; if TRUE, log density is returned
#' @rdname rayleigh2
#' @export
dqrayleigh2 <- function(p, mu, lambda, log=FALSE){
  res <- 2*lambda*sqrt(-log1p(-p)/lambda)*(1-p)
  if(log) return(ifelse(is.finite(res),log(res),res))
  res
}

#' @rdname rayleigh2
#' @export
ffrayleigh2 <- function(p, mu, lambda){
  (2*log1p(-p)+1)/4*lambda^2*(-log1p(-p)/lambda)^(3/2)*(p-1)^2
}


#'Rayleigh distribution
#'
#' @param x numeric vector of data
#' @param sigma numeric scale parameter of the Rayleigh distribution
#' @rdname rayleigh
#' @return density, probability, quantile value or random variable from Rayleigh distribution
#' @export
#'
#' @examples
#' qrayleigh(0.5, 0.5)
drayleigh <- function(x, sigma){
  xbs2 <- x/sigma^2
  xbs2*exp(-xbs2*x/2)
}

#' @rdname rayleigh
#' @export
prayleigh <- function(x, sigma){
  xbs2 <- x/sigma^2
  1-exp(-xbs2*x/2)
}

#' @param p numeric vector of probabilities
#' @rdname rayleigh
#' @export
qrayleigh <- function(p, sigma){
  stopifnot(!any(is.na(p)))
  if(any(p>1)||any(p<0)) stop("p must be before 0 and 1!", call. = FALSE)
  sigma*sqrt(-2*log1p(-p))
}

#' @param n numeric number of random values to draw
#' @rdname rayleigh
#' @export
rrayleigh <- function(n, sigma){
  qrayleigh(runif(n), sigma)
}

#' @rdname rayleigh
#' @export
frayleigh <- function(p, sigma){
  stopifnot(!any(is.na(p)))
  if(any(p>1)||any(p<0)) stop("p must be before 0 and 1!", call. = FALSE)
  sigma/(sqrt(-2*log1p(-p))*(1-p))
}

#' @param log logical; if TRUE, log density is returned
#' @rdname rayleigh
#' @export
dqrayleigh <- function(p, sigma, log=FALSE){
  res <- frayleigh(p, sigma)
  if(log) return(ifelse(is.finite(res),-log(res),res))
  1/res
}

#' @rdname rayleigh
#' @export
ffrayleigh <- function(p, sigma){
  nu <- sigma*(2*log1p(-p)+1)
  de <- sqrt(-8*log1p(-p))*log1p(-p)*(p-1)^2
  nu/de
}




