#' Flattened Generalized Logistic Distribution (FLD)
#'
#' Defines quantile function (ICDF), probability density and quantile density function as well as probability function (CDF) for flattened generalized logistic distribution
#' @param p vector of probabilities
#' @param bt parameter beta of FGLD distribution (scale). Should be non-negative.
#' @param k parameter k of FGLD distribution (shape). Should be non-negative.
#' @param dlt parameter delta of FGLD distribution(mixing parameter). Should be within the interval (0,1), default is 0.5
#'
#' @return vector
#' @name fld
#'
#' @examples
#' qfld(0.1, 0.5, 0.3, 0.5)
#' # centered fld
#' p <- runif(1e4)
#' x <- qfld(p, 0.25, 1)-qfld(0.5, 0.25, 1)
#' hist(x,30)
#' @export
qfld<- function(p, bt, k, dlt=0.5){
  stopifnot("Beta parameter should be non-negative!"=(bt>=0))
  stopifnot("Delta parameter should be between 0 and 1!"=(dlt>=0 && dlt<=1))
  stopifnot("k parameter should be non-negative!"=(k>=0))

  return(bt*(1-dlt)*log(p)-dlt*log(1-p)+k*p)
}

#' @rdname fld
#' @export
ffld <- function(p, bt, k, dlt=0.5){
  stopifnot("Beta parameter should be non-negative!"=(bt>=0))
  stopifnot("Delta parameter should be between 0 and 1!"=(dlt>=0 && dlt<=1))
  stopifnot("k parameter should be non-negative!"=(k>=0))
  return(bt*(1-dlt)/p+dlt/(1-p)+k)
}

#' @rdname fld
#' @export
dqfld <- function(p, bt, k, dlt=0.5){
  stopifnot("Beta parameter should be non-negative!"=(bt>=0))
  stopifnot("Delta parameter should be between 0 and 1!"=(dlt>=0 && dlt<=1))
  stopifnot("k parameter should be non-negative!"=(k>=0))
  return(1/ffld(p, bt, dlt, k))
}

#' @param q vector of quantiles
#' @param tol tolerance value for optimization. Default value 1e-06
#' @rdname fld
#' @importFrom stats uniroot
#' @export
pfld <- function(q, bt, k, dlt=0.5, tol=1e-06){
  stopifnot("Beta parameter should be non-negative!"=(bt>=0))
  stopifnot("Delta parameter should be between 0 and 1!"=(dlt>=0 && dlt<=1))
  stopifnot("k parameter should be non-negative!"=(k>=0))

  afun <- function(x, p) {x - qfld(p, bt, dlt, k)}
  ps <- sapply(q, function(.q) {
    tmp_ps <- NULL
    tmp_ps <- try(stats::uniroot(afun, lower=0, upper = 1, x=.q, tol = tol), silent=TRUE)
    ifelse(is.null(tmp_ps) || inherits(tmp_ps, "try-error"), NA, tmp_ps$root)
    #tmp_ps
  })

  ps[ps < 0] <- 0
  ps[ps > 1] <- 1

  ps[!is.finite(ps)] <- NA
  ps
}
