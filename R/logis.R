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
  z <- stats::qlogis(p, location, scale)
  res <- 1/stats::dlogis(z, location, scale)
  if(log) return(log(res))
  res
}

#' @rdname logis
#' @export
#' @importFrom stats qlogis dlogis
dqlogis <- function(p, location=0, scale=1, log=FALSE){
  ifelse(p<0 | p>1, NA_real_, p) # replace invalid inputs with NA
  z <- stats::qlogis(p, location, scale)
  res <- stats::dlogis(z, location, scale)
  if(log) return(ifelse(is.finite(res),log(res),res))
  res
}
