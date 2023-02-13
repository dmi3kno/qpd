#'Tpkey lambda distribption
#'
#' Defines quantile function (ICDF), and quantile density function as well as probability function (CDF) for Tukey lambda distribution
#' @param p vector of probabilities
#' @param lambda parameter beta of Tukey lambda distribution (scale and shape).
#' @return vector
#' @name tlambda
#'
#' @examples
#' qtlambda(0.1, -0.5)
#' p <- runif(1e4)
#' x <- qtlambda(p, -0.5)+qtlambda(p, 0.5)
#' hist(x,30)
#' @export
qtlambda <- function(p, lambda){
  if(lambda==0) return(qlogis(p))
  1/lambda*(p^lambda-(1-p)^lambda)
}

#' @rdname tlambda
#' @export
ftlambda <- function(p, lambda){
  lam1 <- lambda-1
  if(lambda==0) return(flogis(p))
  p^lam1+(1-p)^lam1
}

#' @param log logical; if TRUE, log density is returnes. Default is FALSE
#' @rdname tlambda
#' @export
dqtlambda <- function(p, lambda, log=FALSE){
  res <- ftlambda(p, lambda=lambda)
  if(log) return(ifelse(is.finite(res),-log(res),res))
  1/res
}

#' @param q vector of quantiles
#' @param ... used by method
#' @param lower,upper the `stats::uniroot` lower and upper end points of the interval to be searched. Defaults are 0 and 1, respectively
#' @param tol the `stats::uniroot` desired accuracy (convergence tolerance). Default value 1e-06
#' @param silent the `base::try` argument. Default is TRUE
#' @param trace integer number passed to `stats::uniroot`; if positive, tracing information is produced. Higher values giving more details.
#' @importFrom stats uniroot
#' @include iqf.R
#' @rdname tlambda
#' @export
ptlambda <- iqf(qtlambda)

#' @param n numeric; number of samples to draw from Tukey Lambda distribution
#' @rdname tlambda
#' @export
rtlambda <- function(n, lambda){
  qtlambda(runif(n), lambda)
}
