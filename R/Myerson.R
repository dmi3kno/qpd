#' Density function for Generalized Lognormal distribution (aka Myerson).
#'
#' @param x,q vector of quantiles. 
#' @param p vector of probabilities
#' @param n number of observations. If `length(n)>1``, the length is taken to be the number required.
#' @param q1,q2,q3 numeric values representing bottom quantile, middle quantile (50th percentile) 
#' and top quantile. Quantiles are assumed to be symmetrical.
#' @param alpha numerical fixed proportion of distribution below the bottom quantile. 
#' Default value is 0.25
#'
#' @return a vector of probabilities of length equal to `length(x)`
#' @name Myerson
#' @export
#' @importFrom stats qnorm dnorm pnorm rnorm runif
#'
#' @examples
#' dMyerson(0.25, 10, 20, 40)
dMyerson <- function(x, q1,q2,q3, alpha=0.25){
  # Probability density function
  br <- (q3-q2)/(q2-q1)

  if(br==1)
    return(stats::dnorm(x, mean=q2, sd=(q3-q2)/stats::qnorm(1-alpha)))

  q <- stats::qnorm(1-alpha)*(log((x-q2)*(br-1)/(q3-q2)+1)/log(br))
  num <- stats::qnorm(1-alpha)*(br-1)
  den <- ((q3-q2)+(x-q2)*(br-1))*log(br)
  num/den*stats::dnorm(q, mean=0, sd=1)
}

#' @rdname Myerson
#'
#' @return a vector of exceedance probabilities of length equal to `length(x)`.
#' @export
#' @importFrom stats qnorm dnorm pnorm rnorm runif
#'
#' @examples
#' pMyerson(0.25, 10, 20, 40)
pMyerson <- function(q, q1,q2,q3, alpha=0.25){
  # cumulative density
  br <- (q3-q2)/(q2-q1)

  if(br==1)
    return(stats::pnorm(q, mean=q2, sd=(q3-q2)/stats::qnorm(1-alpha)))

  nom <- log(1+(q-q2)*(br-1)/(q3-q2))
  den <- log(br)
  qs <- stats::qnorm(1-alpha)*(nom/den)
  stats::pnorm(qs, mean=0, sd=1)
}


#' @rdname Myerson
#' @export
#'
#' @examples
#' qMyerson(0.25, 10, 20, 40)
qMyerson <- function(p, q1,q2,q3, alpha=0.25){

  br <- (q3-q2)/(q2-q1)
  
  if(br==1)
    return(stats::qnorm(p, mean=q2, sd=(q3-q2)/stats::qnorm(1-alpha)))
  
  q2+(q3-q2)*(br^stats::qnorm(p, mean = 0, sd=1/stats::qnorm(1-alpha))-1)/(br-1)
}

#' @param lower,upper vectors of lower and upper bounds of distribution
#' @rdname Myerson
#' @export
#' @importFrom stats qnorm dnorm pnorm rnorm runif
#'
#' @examples
#' rMyerson(1, 10, 20, 40)
rMyerson <- function(n,q1,q2,q3,lower=-Inf, upper=Inf, alpha=0.25){
  # random number generator
  if(length(n)>1) n <-length(n)
  x <- stats::runif(n)

  res <- qMyerson(x, q1, q2, q3, alpha)
  pmin(pmax(res, lower), upper)
}


