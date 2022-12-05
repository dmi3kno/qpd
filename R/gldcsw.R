#' @keywords internal
sGLDcsw <- function(u, chi, xi){

  al <- 0.5*(0.5-xi)/sqrt(xi*(1-xi))
  bt <- 0.5*(chi/sqrt(1-chi^2))

  if(chi==0 && xi==0.5)  return(log(u) - log1p(-u))
  if(chi!=0){
    if(xi==0.5*(1+chi)) return(log(u) - 0.5/al*((1-u)^(2*al)-1))
    if(xi==0.5*(1-chi)) return(0.5/bt*(u^(2*bt)-1)-log1p(-u))
  }
  1/(al+bt)*(u^(al+bt)-1)-1/(al-bt)*((1-u)^(al-bt)-1)
}

#' Generalized Lambda Distribution (GLD) CSW parameterization
#'
#' Quantile function, quantile density, density quantile and inverse quantile
#' functions for GLD distribution with CSW parameterization.
#' `aGLDcsw_mean` and `aGLDcsw_median` are theoretical mean and median, which can be used
#'  for adjusting the quantile likelihood.
#' @param u numeric vector of probabilities
#' @param mu CSW GLD median parameter \eqn{\mu}
#' @param sg CSW GLD interquartile range parameter \eqn{\sg}
#' @param chi CSW GLD assymetry parameter \eqn{\chi}
#' @param xi CSW GLD steepness parameter \eqn{\lambda_4}
#'
#' @return quantiles, QDF, DQF, random samples or probabilities of GLD (CSW parameterization)
#' @rdname GLDcsw
#' @export
#'
#' @examples
#' p_grd <- make_pgrid()
#' qGLDcsw(p_grd, 1, 1, -1/8, 1/32)
qGLDcsw <- function(u, mu, sg, chi, xi){
  stopifnot("chi must be in c(-1,1)!"=((chi>=-1)&&(chi<=1)))
  stopifnot("xi must be in c(0,1)!"=((xi>=0)&&(xi<=1)))
  al <- 0.5*(0.5-xi)/sqrt(xi*(1-xi))
  bt <- 0.5*(chi/sqrt(1-chi^2))

  lb <- (-Inf); ub <- Inf
  if (xi<0.5*(1+chi)) lb <- (-1/(al+bt))
  if (xi<0.5*(1-chi)) ub <- (1/(al-bt))

  iqr <- sGLDcsw(0.75, chi, xi)-sGLDcsw(0.25, chi, xi)
  med <- sGLDcsw(0.5, chi, xi)

  s <-   sGLDcsw(u, chi, xi)
  s[u==0] <- lb
  s[u==1] <- ub

  mu + sg*(s-med)/iqr
}

#' @rdname GLDcsw
#' @export
fGLDcsw <- function(u, mu, sg, chi, xi){
  al <- 0.5*(0.5-xi)/sqrt(xi*(1-xi))
  bt <- 0.5*(chi/sqrt(1-chi^2))
  iqr <- sGLDcsw(0.75, chi, xi)-sGLDcsw(0.25, chi, xi)

  sg/iqr*(u^(al+bt-1)+(1-u)^(al-bt-1))
}

#' @param log should the log density be returned. Default=FALSE
#' @rdname GLDcsw
#' @export
dqGLDcsw <- function(u, mu, sg, chi, xi, log=FALSE){
  res <- fGLDcsw(u, mu, sg, chi, xi)
  if(log) return(ifelse(is.finite(res),-log(res),res))
  1/res
}

#' @rdname GLDcsw
#' @export
#' @param n numeric number of samples to draw
rGLDcsw <- function(n, mu, sg, chi, xi){
  qGLDcsw(runif(n), mu, sg, chi, xi)
}

#' @param q vector of quantiles
#' @param ... used by method
#' @param lower,upper the `stats::uniroot` lower and upper end points of the interval to be searched. Defaults are 0 and 1, respectively
#' @param tol the `stats::uniroot` desired accuracy (convergence tolerance). Default value 1e-06
#' @param silent the `base::try` argument. Default is TRUE
#' @param trace integer number passed to `stats::uniroot`; if positive, tracing information is produced. Higher values giving more details.
#' @rdname GLDcsw
#' @importFrom stats uniroot
#' @include iqf.R
#' @export
pGLDcsw <- iqf(qGLDcsw)

