#' Govindarajulu probability distribution
#'
#' Defines quantile function (ICDF), probability density and quantile density function as well as probability function (CDF) for Govindarajulu distribution
#' @param p vector of probabilities
#' @param sg parameter sigma of Govindarajulu distribution
#' @param bt parameter beta of Govindarajulu distribution
#'
#' @return vector
#' @name govindarajulu
#'
#' @examples
#' qgovindarajulu(0.1, 1, 5)
#' @export
qgovindarajulu <- function(p, sg, bt){
  stopifnot(sg>0 & bt>0)
  stopifnot(p>=0 & p<=1 | is.na(p))
  bp1 <- bt+1
  sg*(bp1*p^bt-bt*p^bp1)
}

#' @rdname govindarajulu
#' @export
fgovindarajulu <- function(p, sg, bt){
  stopifnot(sg>0 & bt>0)
  stopifnot(all(p>=0 & p<=1 || is.na(p)))
  bp1 <- bt+1
  #sg*bt*bp1*p^bt*(1-p) # wrong formula in Nair et al 2020
  sg*bt*bp1*p^(bt-1)*(1-p)
}

#' @param log logical; if TRUE, log density is returned. Default is FALSE
#' @rdname govindarajulu
#' @export
dqgovindarajulu <- function(p, sg, bt, log=FALSE){
  stopifnot(sg>0 & bt>0)
  res <- 1/fgovindarajulu(p, sg, bt)
  if(log) return(ifelse(is.finite(res),log(res),res))
  res
}

#' @param q vector of quantiles
#' @param tol tolerance value for optimization. Default value 1e-06
#' @rdname govindarajulu
#' @importFrom stats uniroot
#' @export
pgovindarajulu <- function(q, sg, bt, tol=1e-06){
  stopifnot(sg>0 & bt>0)

  afun <- function(x, p) {x - qgovindarajulu(p, sg, bt)}
  ps <- sapply(q, function(.q) {
    tmp_ps <- NULL
    tmp_ps <- try(stats::uniroot(afun, lower=0, upper = 1, x=.q, tol = tol), silent=TRUE)
    ifelse(is.null(tmp_ps) || inherits(tmp_ps, "try-error"), NA, tmp_ps$root)
    #tmp_ps
  })

  ps[ps < 0] <- 0
  ps[ps>sg] <- 1

  ps[!is.finite(ps)] <- NA
  ps
}

#' @param x numeric vector
#' @rdname govindarajulu
#' @export
dgovindarajulu <- function(x, sg, bt, tol=1e-06){
  tmp <- 1 / (sg*bt*(bt+1))

  p <- pgovindarajulu(x, sg, bt, tol)
  p[x<0] <- NA
  p[x>sg] <- NA
  f <- tmp * p^(1-bt)/(1-p)
  f[! is.finite(f)] <- NA
  f
}
