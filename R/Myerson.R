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
  r <- (q3-q2)
  b <- r/(q2-q1)
  qn1ma <-stats::qnorm(1-alpha)
  if(b==1)
    return(stats::dnorm(x, mean=q2, sd=r/qn1ma))

  psi <- qn1ma*(log(1+(x-q2)*(b-1)/r)/log(b))
  num <- qn1ma*(b-1)
  den <- (r+(x-q2)*(b-1))*log(b)
  num/den*stats::dnorm(psi)
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
  r <- (q3-q2)
  b <- r/(q2-q1)
  qn1ma <-stats::qnorm(1-alpha)

  if(b==1)
    return(stats::pnorm(q, mean=q2, sd=r/qn1ma))

  nom <- log(1+(q-q2)*(b-1)/r)
  den <- log(b)
  psi <- qn1ma*(nom/den)
  stats::pnorm(psi)
}

#' @rdname Myerson
#' @export
#'
#' @examples
#' qMyerson(0.25, 10, 20, 40)
qMyerson <- function(p, q1,q2,q3, alpha=0.25){
  r <- (q3-q2)
  b <- r/(q2-q1)
  qn1ma <-stats::qnorm(1-alpha)
  k <- stats::qnorm(p)/qn1ma

  if(b==1)
    return(q2+r*k)

  q2+r*(b^k-1)/(b-1)
}

#' @rdname Myerson
#' @export
#'
#' @examples
#' fMyerson(0.25, 10, 20, 40)
fMyerson <- function(p, q1,q2,q3, alpha=0.25){
  r <- (q3-q2)
  b <- r/(q2-q1)
  qn1ma <-stats::qnorm(1-alpha)
  k <- stats::qnorm(p)/qn1ma

  if(b==1)
    return(r*fnorm(p)/qn1ma)

  r*log(b)*(b^k)*fnorm(p)/(b-1)/qn1ma
}

#' @rdname Myerson
#' @export
#'
#' @examples
#' fMyerson(0.25, 10, 20, 40)
dqMyerson <- function(p, q1,q2,q3, alpha=0.25){
  1/fMyerson(p, q1,q2,q3, alpha)
}

#' @rdname Myerson
#' @export
#' @importFrom stats qnorm dnorm pnorm rnorm runif
#'
#' @examples
#' rMyerson(1, 10, 20, 40)
rMyerson <- function(n,q1,q2,q3,alpha=0.25){
  # random number generator
  if(length(n)>1) n <-length(n)
  qMyerson(stats::runif(n), q1, q2, q3, alpha)
}


