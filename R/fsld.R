#' Flattened Skew-Logistic Distribution (FSLD)
#'
#' Defines quantile function (ICDF), probability density and quantile density function as well as probability function (CDF) for generalized flattened logistic distribution
#' @param p vector of probabilities
#' @param bt parameter beta of GFLD distribution (scale). Should be non-negative.
#' @param k parameter k of GFLD distribution (shape). Should be non-negative.
#' @param dlt parameter delta of GFLD distribution(mixing parameter). Should be within the interval (0,1), default is 0.5
#' @param a location parameter alpha of GFLD distribution(location parameter), default is 0
#'
#' @return vector
#' @name fsld
#'
#' @examples
#' qfsld(0.1, 0.5, 0.3, 0.5)
#' # centered gfld
#' p <- runif(1e4)
#' x <- qfsld(p, 0.25, 1)-qfsld(0.5, 0.25, 1)
#' hist(x,30)
#' @export
qfsld<- function(p, bt, k, dlt=0.5, a=0){
  stopifnot("Beta parameter should be non-negative!"=(bt>=0))
  stopifnot("Delta parameter should be between 0 and 1!"=(dlt>=0 && dlt<=1))
  stopifnot("k parameter should be non-negative!"=(k>=0))
  Qslogis <- (1-dlt)*log(p)-dlt*log(1-p)
  Qunf <- k*p
  return(a+bt*(Qslogis+Qunf))
}

#' @rdname fsld
#' @export
ffsld <- function(p, bt, k, dlt=0.5, a=0){
  stopifnot("Beta parameter should be non-negative!"=(bt>=0))
  stopifnot("Delta parameter should be between 0 and 1!"=(dlt>=0 && dlt<=1))
  stopifnot("k parameter should be non-negative!"=(k>=0))
  return(bt*((1-dlt)/p+dlt/(1-p)+k))
}

#' @param log logical; if TRUE, log density is returnes. Default is FALSE
#' @rdname fsld
#' @export
dqfsld <- function(p, bt, k, dlt=0.5, a=0, log=FALSE){
  res <- 1/ffsld(p, bt=bt, k=k, dlt=dlt, a=a)
  if(log) return(ifelse(is.finite(res),log(res),res))
  res
}

#' @param n numeric; number of samples to draw from FLD distribution
#' @rdname fsld
#' @export
rfsld <- function(n, bt, k, dlt=0.5, a=0){
  qfsld(runif(n), bt, k, dlt=dlt, a=a)
}

#' @param q vector of quantiles
#' @param tol tolerance value for optimization. Default value 1e-06
#' @rdname fsld
#' @importFrom stats uniroot
#' @export
pfsld <- function(q, bt, k, dlt=0.5, a=0, tol=1e-06){
  stopifnot("Beta parameter should be non-negative!"=(bt>=0))
  stopifnot("Delta parameter should be between 0 and 1!"=(dlt>=0 && dlt<=1))
  stopifnot("k parameter should be non-negative!"=(k>=0))

  afun <- function(x, p) {x - qfsld(p, bt=bt, k=k, dlt=dlt, a=a)}
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

#' Flattened Logistic Distribution (FLD)
#'
#' Defines quantile function (ICDF), probability density and quantile density function as well as probability function (CDF) for flattened logistic distribution
#' @param p vector of probabilities
#' @param bt parameter beta of GFLD distribution (scale). Should be non-negative.
#' @param k parameter k of GFLD distribution (shape). Should be non-negative.
#' @param a location parameter alpha of GFLD distribution(location parameter), default is 0
#'
#' @return vector
#' @name fld
#'
#' @examples
#' qfld(0.1, 0.5, 0.3, 0.5)
#' # centered gfld
#' p <- runif(1e4)
#' x <- qfld(p, 0.25, 1)-qfld(0.5, 0.25, 1)
#' hist(x,30)
#' @export
qfld<- function(p, bt, k, a=0){
  stopifnot("Beta parameter should be non-negative!"=(bt>=0))
  stopifnot("k parameter should be non-negative!"=(k>=0))
  Qlogis <- log(p/(1-p))
  Qunf <- k*p
  return(a+bt*(Qlogis+Qunf))
}

#' @rdname fld
#' @export
ffld <- function(p, bt, k, a=0){
  stopifnot("Beta parameter should be non-negative!"=(bt>=0))
  stopifnot("k parameter should be non-negative!"=(k>=0))
  # bt*(k*p^2-k*p-1)/((p-1)*p)
  omp <- 1-p
  bt*(omp/p*(p/omp^2+1/omp)+k)
}

#' @param log logical; if TRUE, log density is returnes. Default is FALSE
#' @rdname fld
#' @export
dqfld <- function(p, bt, k, a=0, log=FALSE){
  res <- 1/ffld(p, bt=bt, k=k, a=a)
  if(log) return(ifelse(is.finite(res),log(res),res))
  res
}

#' @param n numeric; number of samples to draw from FLD distribution
#' @rdname fld
#' @export
rfld <- function(n, bt, k, a=0){
  qfld(runif(n), bt, k, a=a)
}

#' @param q vector of quantiles
#' @param tol tolerance value for optimization. Default value 1e-06
#' @rdname fld
#' @importFrom stats uniroot
#' @export
pfld <- function(q, bt, k, a=0, tol=1e-06){
  stopifnot("Beta parameter should be non-negative!"=(bt>=0))
  stopifnot("k parameter should be non-negative!"=(k>=0))

  afun <- function(x, p) {x - qfld(p, bt=bt, k=k, a=a)}
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

