#' Density function for Generalized Lognormal distribution (aka probit Myerson).
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

#' @rdname Myerson
#' @export
#'
#' @examples
#' dMyerson(15, 10, 20, 40)
dMyerson <- function(x, q1,q2,q3, alpha=0.25){
  # Probability density function
  r <- (q3-q2)
  b <- r/(q2-q1)
  qn1ma <-stats::qnorm(1-alpha)
  if(b==1)
    return(stats::dnorm(x, mean=q2, sd=r/qn1ma))

  psi <- qn1ma*(log1p((x-q2)*(b-1)/r)/log(b))
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

  nom <- log1p((q-q2)*(b-1)/r)
  den <- log(b)
  psi <- qn1ma*(nom/den)
  stats::pnorm(psi)
}

#' Density function for logit Myerson distribution (Myerson distribution with logit base function).
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
#' @name logitMyerson
#' @export
#' @importFrom stats qnorm dnorm pnorm rnorm runif
#'
#' @examples
#' qlogitMyerson(0.25, 10, 20, 40)
qlogitMyerson <- function(p, q1,q2,q3, alpha=0.25){
  r <- (q3-q2)
  b <- r/(q2-q1)
  qn1ma <-logit(1-alpha)
  k <- logit(p)/qn1ma

  if(b==1)
    return(q2+r*k)

  q2+r*(b^k-1)/(b-1)
}

#' @rdname logitMyerson
#' @export
#'
#' @examples
#' flogitMyerson(0.25, 10, 20, 40)
flogitMyerson <- function(p, q1,q2,q3, alpha=0.25){
  r <- (q3-q2)
  b <- r/(q2-q1)
  qn1ma <-logit(1-alpha)
  k <- logit(p)/qn1ma

  if(b==1)
    return(r/(p*(1-p)*qn1ma))

  r*log(b)*(b^k)/(p*(1-p)*(b-1)*qn1ma)
}

#' @rdname logitMyerson
#' @export
#'
#' @examples
#' dqlogitMyerson(0.25, 10, 20, 40)
dqlogitMyerson <- function(p, q1,q2,q3, alpha=0.25){
  1/flogitMyerson(p, q1,q2,q3, alpha)
}

#' @rdname logitMyerson
#' @export
#' @importFrom stats qnorm dnorm pnorm rnorm runif
#'
#' @examples
#' rlogitMyerson(1, 10, 20, 40)
rlogitMyerson <- function(n,q1,q2,q3,alpha=0.25){
  # random number generator
  if(length(n)>1) n <-length(n)
  qlogitMyerson(stats::runif(n), q1, q2, q3, alpha)
}


#' @rdname logitMyerson
#'
#' @return a vector of exceedance probabilities of length equal to `length(x)`.
#' @export
#' @importFrom stats qnorm dnorm pnorm rnorm runif
#'
#' @examples
#' plogitMyerson(0.25, 10, 20, 40)
plogitMyerson <- function(q, q1,q2,q3, alpha=0.25){
  # cumulative density
  r <- (q3-q2)
  b <- r/(q2-q1)
  qn1ma <-logit(1-alpha)
  if(b==1)
    return(invlogit(qn1ma*(q-q2)/r))

  nom <- log1p((q-q2)*(b-1)/r)
  den <- log(b)
  psi <- qn1ma*(nom/den)
  invlogit(psi)
}


#' @rdname logitMyerson
#' @export
#'
#' @examples
#' dlogitMyerson(15, 10, 20, 40)
dlogitMyerson <- function(x, q1,q2,q3, alpha=0.25){
  # Probability density function
  r <- (q3-q2)
  b <- r/(q2-q1)
  qn1ma <-logit(1-alpha)
  psi <- qn1ma*(x-q2)/r
  if(b==1)
    return(stats::dlogis(psi))

  psi <- qn1ma*(log(1+(x-q2)*(b-1)/r)/log(b))
  num <- qn1ma*(b-1)
  den <- (r+(x-q2)*(b-1))*log(b)
  num/den*stats::dlogis(psi)
}


#' Density function for Cauchy Myerson distribution (Myerson distribution with Cauchy base function).
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
#' @name cauchyMyerson
#' @export
#' @importFrom stats qnorm dnorm pnorm rnorm runif
#'
#' @examples
#' qcauchyMyerson(0.25, 10, 20, 40)
qcauchyMyerson <- function(p, q1,q2,q3, alpha=0.25){
  r <- (q3-q2)
  b <- r/(q2-q1)
  qn1ma <-tan(pi*(0.5-alpha))
  k <- tan(pi*(p-0.5))/qn1ma

  if(b==1)
    return(q2+r*k)

  q2+r*(b^k-1)/(b-1)
}

#' @rdname cauchyMyerson
#' @export
#'
#' @examples
#' fcauchyMyerson(0.25, 10, 20, 40)
fcauchyMyerson <- function(p, q1,q2,q3, alpha=0.25){
  r <- (q3-q2)
  b <- r/(q2-q1)
  qn1ma <-tan(pi*(0.5-alpha))
  k <- tan(pi*(p-0.5))/qn1ma

  if(b==1)
    return(r*pi*(1+tan(pi*(p-0.5))^2)/qn1ma)

  r*log(b)*(b^k)*pi*(1+tan(pi*(p-0.5))^2)/((b-1)*qn1ma)
}

#' @rdname cauchyMyerson
#' @export
#'
#' @examples
#' dqcauchyMyerson(0.25, 10, 20, 40)
dqcauchyMyerson <- function(p, q1,q2,q3, alpha=0.25){
  1/fcauchyMyerson(p, q1,q2,q3, alpha)
}

#' @rdname cauchyMyerson
#' @export
#' @importFrom stats qnorm dnorm pnorm rnorm runif
#'
#' @examples
#' rcauchyMyerson(1, 10, 20, 40)
rcauchyMyerson <- function(n,q1,q2,q3,alpha=0.25){
  # random number generator
  if(length(n)>1) n <-length(n)
  qcauchyMyerson(stats::runif(n), q1, q2, q3, alpha)
}


#' @rdname cauchyMyerson
#'
#' @return a vector of exceedance probabilities of length equal to `length(x)`.
#' @export
#' @importFrom stats qnorm dnorm pnorm rnorm runif
#'
#' @examples
#' pcauchyMyerson(0.25, 10, 20, 40)
pcauchyMyerson <- function(q, q1,q2,q3, alpha=0.25){
  # cumulative density
  r <- (q3-q2)
  b <- r/(q2-q1)
  qn1ma <-tan(pi*(0.5-alpha))
  if(b==1)
    return(stats::pcauchy(qn1ma*(q-q2)/r))

  nom <- log1p((q-q2)*(b-1)/r)
  den <- log(b)
  psi <- qn1ma*(nom/den)
  stats::pcauchy(psi)
}


#' @rdname cauchyMyerson
#' @export
#'
#' @examples
#' dcauchyMyerson(15, 10, 20, 40)
dcauchyMyerson <- function(x, q1,q2,q3, alpha=0.25){
  # Probability density function
  r <- (q3-q2)
  b <- r/(q2-q1)
  qn1ma <-tan(pi*(0.5-alpha))
  psi <- qn1ma*(x-q2)/r
  if(b==1)
    return(stats::dcauchy(psi))

  psi <- qn1ma*(log(1+(x-q2)*(b-1)/r)/log(b))
  num <- qn1ma*(b-1)
  den <- (r+(x-q2)*(b-1))*log(b)
  num/den*stats::dcauchy(psi)
}


#' Density function for Hyperbolic Secant Myerson distribution (Myerson distribution with hyperbolic secant base function).
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
#' @name sechMyerson
#' @export
#' @importFrom stats qnorm dnorm pnorm rnorm runif
#'
#' @examples
#' qsechMyerson(0.25, 10, 20, 40)
qsechMyerson <- function(p, q1,q2,q3, alpha=0.25){
  r <- (q3-q2)
  b <- r/(q2-q1)
  qn1ma <-2/pi*log(tan(pi/2*(1-alpha))) #sech QF
  k <- 2/pi*log(tan(pi/2*p))/qn1ma

  if(b==1)
    return(q2+r*k)

  q2+r*(b^k-1)/(b-1)
}

#' @rdname sechMyerson
#' @export
#'
#' @examples
#' fsechMyerson(0.25, 10, 20, 40)
fsechMyerson <- function(p, q1,q2,q3, alpha=0.25){
  r <- (q3-q2)
  b <- r/(q2-q1)
  qn1ma <-2/pi*log(tan(pi/2*(1-alpha)))
  k <- 2/pi*log(tan(pi/2*p))/qn1ma

  if(b==1)
    return(r*pi*(1+tan(pi*p/2)^2)/(pi*qn1ma*tan(pi*p/2)))

  r*log(b)*(b^k)*pi*(1+tan(pi*p/2)^2)/((b-1)*pi*qn1ma*tan(pi*p/2))
}

#' @rdname sechMyerson
#' @export
#'
#' @examples
#' dqsechMyerson(0.25, 10, 20, 40)
dqsechMyerson <- function(p, q1,q2,q3, alpha=0.25){
  1/fsechMyerson(p, q1,q2,q3, alpha)
}

#' @rdname sechMyerson
#' @export
#' @importFrom stats qnorm dnorm pnorm rnorm runif
#'
#' @examples
#' rsechMyerson(1, 10, 20, 40)
rsechMyerson <- function(n,q1,q2,q3,alpha=0.25){
  # random number generator
  if(length(n)>1) n <-length(n)
  qsechMyerson(stats::runif(n), q1, q2, q3, alpha)
}


#' @rdname sechMyerson
#'
#' @return a vector of exceedance probabilities of length equal to `length(x)`.
#' @export
#' @importFrom stats qnorm dnorm pnorm rnorm runif
#'
#' @examples
#' psechMyerson(0.25, 10, 20, 40)
psechMyerson <- function(q, q1,q2,q3, alpha=0.25){
  # cumulative density
  r <- (q3-q2)
  b <- r/(q2-q1)
  qn1ma <-2/pi*log(tan(pi/2*(1-alpha)))
  if(b==1)
    return(2/pi*atan(exp(pi/2*(qn1ma*(q-q2)/r))))

  nom <- log1p((q-q2)*(b-1)/r)
  den <- log(b)
  psi <- qn1ma*(nom/den)
  2/pi*atan(exp(pi/2*psi))
}


#' @rdname sechMyerson
#' @export
#'
#' @examples
#' dsechMyerson(15, 10, 20, 40)
dsechMyerson <- function(x, q1,q2,q3, alpha=0.25){
  # Probability density function
  r <- (q3-q2)
  b <- r/(q2-q1)
  qn1ma <-2/pi*log(tan(pi/2*(1-alpha)))
  psi <- qn1ma*(x-q2)/r
  if(b==1)
    return(1/2*sech(pi/2*psi))

  psi <- qn1ma*(log1p((x-q2)*(b-1)/r)/log(b))
  num <- qn1ma*(b-1)
  den <- (r+(x-q2)*(b-1))*log(b)
  num/den*1/2*sech(pi/2*psi)
}







