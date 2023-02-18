#' Arctan Power distribution
#'
#' Defines quantile function (ICDF), and quantile density function as well as probability function (CDF) for Arctan Power distribution
#' @param p vector of probabilities
#' @param alpha shape parameter alpha of Arctan Power Distribution.
#' @param beta shape parameter beta of Arctan Power Distribution.
#' @return vector
#' @rdname atanpow
#' @references Nasiru S, Abubakari AG, Chesneau C. The Arctan Power Distribution: Properties, Quantile and Modal Regressions with Applications to Biomedical Data. Mathematical and Computational Applications. 2023; 28(1):25. https://doi.org/10.3390/mca28010025
#'
#' @examples
#' qatanpow(0.1, 2,2)
#' @export
qatanpow <- function(p, alpha=1, beta=1){
  stopifnot(alpha>0, beta>0)
  (tan(p*atan(alpha))/alpha)^(1/beta)
}


#' @rdname atanpow
#' @export
fatanpow <- function(p, alpha=1, beta=1){
  stopifnot(alpha>0, beta>0)
  aa <- atan(alpha)
  num <- aa*(1+tan(aa*p)^2)*qatanpow(p, alpha, beta)
  den <- beta*tan(aa*p)
  num/den
}

#' @param log logical; if TRUE, log density is returned. Default is FALSE
#' @rdname atanpow
#' @export
dqatanpow <- function(p, alpha=1, beta=1, log=FALSE){
  res <- fatanpow(p, alpha, beta)
  if(log) return(ifelse(is.finite(res),-log(res),res))
  1/res
}

#' @param x data
#' @rdname atanpow
#' @export
patanpow <- function(x, alpha=1, beta=1){
  stopifnot(alpha>0, beta>0)
  atan(alpha*x^beta)/atan(alpha)
}

#' @rdname atanpow
#' @export
datanpow <- function(x, alpha=1, beta=1){
  num <- alpha*beta*x^(beta-1)
  den <- atan(alpha)*(1+alpha^2*x^(2*beta))
  num/den
}



