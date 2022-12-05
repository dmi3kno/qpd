#' Quantile function for Gilchrist's skew-logistic distribution
#'
#' Calculate QF, QDF, DQF `qslogis`, `flogis`, `dqlogis` for the skew-logistic distribution
#' @param p numeric probability
#' @param location numeric parameter of logistic distribution. Default is 0
#' @param scale numeric parameter of logistic distribution. Default is 1
#' @param dlt parameter delta of SLD distribution(mixing parameter). Should be within the interval (0,1), default is 0.5
#' @param log logical, should the result be returned as a log. Default is FALSE
#'
#' @return returns QF (`qslogis`), QDF (`fslogis`), DQF (`dqslogis`) or random variate of skew-logistic distribution
#' @export
#' @rdname slogis
#' @examples
#' p_grd <- make_pgrid()
#' fslogis(p_grd)
#' dqslogis(p_grd)
qslogis <- function(p, location=0, scale=1, dlt=0.5){
  location+scale*(dlt*(-log1p(-p)+(1-dlt)*log(p)))
}

#' @rdname slogis
#' @export
fslogis <- function(p, location=0, scale=1, dlt=0.5, log=FALSE){
  res <- scale*((1-dlt)/p+dlt/(1-p))
  if(log) return(log(res))
  res
}

#' @rdname slogis
#' @export
dqslogis <- function(p, location=0, scale=1, dlt=0.5, log=FALSE){
  res <- fslogis(p, location, scale, dlt)
  if(log) return(ifelse(is.finite(res),-log(res),res))
  1/res
}

#' @rdname slogis
#' @param n number of random variates to generate
#' @export
rslogis <- function(n, location, scale, dlt){
  qslogis(runif(n), location, scale, dlt)
}



#' Derivative quantile functions for logistic distribution
#'
#' Calculate QDF `flogis` and DQF `dqlogis` for the logistic distribution
#' @param p numeric probability
#' @param location numeric parameter of logistic distribution. Default is 0
#' @param scale numeric parameter of logistic distribution. Default is 1
#' @param log logical, should the result be returned as a log. Default is FALSE
#'
#' @return returns QDF (`flogis`) and DQF (`dqlogis`) of logistic distribution
#' @export
#' @importFrom stats qlogis dlogis
#' @rdname logis
#' @examples
#' p_grd <- make_pgrid()
#' flogis(p_grd)
#' dqlogis(p_grd)
flogis <- function(p, location=0, scale=1, log=FALSE){
  ifelse(p<0 | p>1, NA_real_, p) # replace invalid inputs with NA
  res <- -scale/(p*(p-1))
  if(log) return(log(res))
  res
}

#' @rdname logis
#' @export
#' @importFrom stats qlogis dlogis
dqlogis <- function(p, location=0, scale=1, log=FALSE){
  ifelse(p<0 | p>1, NA_real_, p) # replace invalid inputs with NA
  res <- flogis(p, location, scale)
  if(log) return(ifelse(is.finite(res),-log(res),res))
  1/res
}
