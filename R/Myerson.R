#' @keywords internal
sMyerson_rba <- function(p,r,b,alpha){
  bt <- r/b
  qn1ma <-stats::qnorm(1-alpha)
  k <- stats::qnorm(p)/qn1ma
  if(bt==1) return(r*k)
  r*(bt^k-1)/(bt-1)
}

fMyerson_rba <- function(p,r,b,alpha,log){
  bt <- r/b
  qn1ma <-stats::qnorm(1-alpha)
  k <- stats::qnorm(p)/qn1ma
  if(bt==1){
    res <- r*fnorm(p)/qn1ma
  } else {
    res <- r*log(bt)*(bt^k)*fnorm(p)/(bt-1)/qn1ma
  }
  if(log) return(log(res))
  res
}

#' Density function for Generalized Lognormal distribution (aka probit Myerson).
#'
#' @param x,q vector of quantiles.
#' @param p vector of probabilities
#' @param n number of observations. If `length(n)>1``, the length is taken to be the number required.
#' @param q1,q2,q3 numeric values representing bottom quantile, middle quantile (50th percentile)
#' and top quantile. Quantiles are assumed to be symmetrical.
#' @param alpha numerical fixed proportion of distribution below the bottom quantile.
#' Default value is 0.25
#' @param log logical; if TRUE, log density is returned
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
  b <- (q2-q1)
  q2+sMyerson_rba(p,r,b,alpha)
}

#' @rdname Myerson
#' @export
#'
#' @examples
#' fMyerson(0.25, 10, 20, 40)
fMyerson <- function(p, q1,q2,q3, alpha=0.25, log=FALSE){
  r <- (q3-q2)
  b <- (q2-q1)
  fMyerson_rba(p,r,b, alpha,log)
}

#' @rdname Myerson
#' @export
#'
#' @examples
#' fMyerson(0.25, 10, 20, 40)
dqMyerson <- function(p, q1,q2,q3, alpha=0.25, log=FALSE){
  res <- fMyerson(p, q1,q2,q3, alpha, log=FALSE)
  if(log) return(ifelse(is.finite(res),-log(res),res))
  1/res
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
dMyerson <- function(x, q1,q2,q3, alpha=0.25, log=FALSE){
  # Probability density function
  r <- (q3-q2)
  b <- (q2-q1)
  bt <- r/b
  qn1ma <-stats::qnorm(1-alpha)
  if(bt==1){
    res <- stats::dnorm(x, mean=q2, sd=r/qn1ma)
  }else{
    psi <- qn1ma*(log1p((x-q2)*(bt-1)/r)/log(bt))
    num <- qn1ma*(bt-1)
    den <- (r+(x-q2)*(bt-1))*log(bt)
    res <- num/den*stats::dnorm(psi)
  }
  if(log) return(log(res))
  res
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
  b <- (q2-q1)
  bt <- r/b
  qn1ma <-stats::qnorm(1-alpha)

  if(bt==1)
    return(stats::pnorm(q, mean=q2, sd=r/qn1ma))

  nom <- log1p((q-q2)*(bt-1)/r)
  den <- log(bt)
  psi <- qn1ma*(nom/den)
  stats::pnorm(psi)
}

########################### LOGIT MYERSON ###################################

slogitMyerson_rba <- function(p,r,b,alpha){
  bt <- r/b
  qn1ma <-logit(1-alpha)
  k <- logit(p)/qn1ma
  if(bt==1) return(r*k)
  r*(bt^k-1)/(bt-1)
}

flogitMyerson_rba<- function(p,r,b,alpha,log){
  bt <- r/b
  qn1ma <-logit(1-alpha)
  k <- logit(p)/qn1ma
  if(bt==1){
    res <- r/(p*(1-p)*qn1ma)
  } else {
    res <- r*log(bt)*(bt^k)/(p*(1-p)*(bt-1)*qn1ma)
  }

  if(log) return(log(res))
  res
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
#' @param log logical; if TRUE, log density is returned
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
  b <- (q2-q1)
  q2+slogitMyerson_rba(p,r,b,alpha)
}

#' @rdname logitMyerson
#' @export
#'
#' @examples
#' flogitMyerson(0.25, 10, 20, 40)
flogitMyerson <- function(p, q1,q2,q3, alpha=0.25, log=FALSE){
  r <- (q3-q2)
  b <- (q2-q1)
  flogitMyerson_rba(p,r,b,alpha,log)
}

#' @rdname logitMyerson
#' @export
#'
#' @examples
#' dqlogitMyerson(0.25, 10, 20, 40)
dqlogitMyerson <- function(p, q1,q2,q3, alpha=0.25, log=FALSE){
   res <- flogitMyerson(p, q1,q2,q3, alpha, log=FALSE)
   if(log) return(ifelse(is.finite(res),-log(res),res))
   1/res
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
  b <- (q2-q1)
  bt <- r/b
  qn1ma <-logit(1-alpha)
  if(bt==1)
    return(invlogit(qn1ma*(q-q2)/r))

  nom <- log1p((q-q2)*(bt-1)/r)
  den <- log(bt)
  psi <- qn1ma*(nom/den)
  invlogit(psi)
}


#' @rdname logitMyerson
#' @export
#'
#' @examples
#' dlogitMyerson(15, 10, 20, 40)
dlogitMyerson <- function(x, q1,q2,q3, alpha=0.25, log=FALSE){
  # Probability density function
  r <- (q3-q2)
  b <- (q2-q1)
  bt <- r/b
  qn1ma <-logit(1-alpha)
  psi <- qn1ma*(x-q2)/r
  if(bt==1){
    res <- stats::dlogis(psi)
  }else{
    psi <- qn1ma*(log(1+(x-q2)*(bt-1)/r)/log(bt))
    num <- qn1ma*(bt-1)
    den <- (r+(x-q2)*(bt-1))*log(bt)
    res <- num/den*stats::dlogis(psi)
  }

  if(log) return(log(res))
  res
}

####################### CAUCHY MYERSON ###################################
scauchyMyerson_rba <- function(p,r,b,alpha){
  bt <- r/b
  qn1ma <-tan(pi*(0.5-alpha))
  k <- tan(pi*(p-0.5))/qn1ma
  if(bt==1) return(r*k)
  r*(bt^k-1)/(bt-1)
}

fcauchyMyerson_rba <- function(p,r,b,alpha,log){
  bt <- r/b
  qn1ma <-tan(pi*(0.5-alpha))
  k <- tan(pi*(p-0.5))/qn1ma
  if(bt==1){
    res <- r*pi*(1+tan(pi*(p-0.5))^2)/qn1ma
  } else {
    res <- r*log(bt)*(bt^k)*pi*(1+tan(pi*(p-0.5))^2)/((bt-1)*qn1ma)
  }

  if(log) return(log(res))
  res
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
#' @param log logical; if TRUE, log density is returned
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
  b <- (q2-q1)
  q2+scauchyMyerson_rba(p,r,b,alpha)
}

#' @rdname cauchyMyerson
#' @export
#'
#' @examples
#' fcauchyMyerson(0.25, 10, 20, 40)
fcauchyMyerson <- function(p, q1,q2,q3, alpha=0.25, log=FALSE){
  r <- (q3-q2)
  b <- (q2-q1)
  fcauchyMyerson_rba(p,r,b,alpha,log)
}

#' @rdname cauchyMyerson
#' @export
#'
#' @examples
#' dqcauchyMyerson(0.25, 10, 20, 40)
dqcauchyMyerson <- function(p, q1,q2,q3, alpha=0.25, log=FALSE){
  res <- fcauchyMyerson(p, q1,q2,q3, alpha, log=FALSE)
  if(log) return(ifelse(is.finite(res),-log(res),res))
  1/res
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
  b <- (q2-q1)
  bt <- r/b
  qn1ma <-tan(pi*(0.5-alpha))
  if(bt==1)
    return(stats::pcauchy(qn1ma*(q-q2)/r))

  nom <- log1p((q-q2)*(bt-1)/r)
  den <- log(bt)
  psi <- qn1ma*(nom/den)
  stats::pcauchy(psi)
}


#' @rdname cauchyMyerson
#' @export
#'
#' @examples
#' dcauchyMyerson(15, 10, 20, 40)
dcauchyMyerson <- function(x, q1,q2,q3, alpha=0.25, log=FALSE){
  # Probability density function
  r <- (q3-q2)
  b <- (q2-q1)
  bt <- r/b
  qn1ma <-tan(pi*(0.5-alpha))
  psi <- qn1ma*(x-q2)/r
  if(bt==1){
    res <- stats::dcauchy(psi)
  }else{
    psi <- qn1ma*(log(1+(x-q2)*(bt-1)/r)/log(bt))
    num <- qn1ma*(bt-1)
    den <- (r+(x-q2)*(bt-1))*log(bt)
    res <- num/den*stats::dcauchy(psi)
  }

  if(log) return(log(res))
  res
}

######################### SECH MYERSON #######################################

ssechMyerson_rba <- function(p,r,b,alpha){
  bt <- r/b
  qn1ma <-2/pi*log(tan(pi/2*(1-alpha))) #sech QF
  k <- 2/pi*log(tan(pi/2*p))/qn1ma
  if(bt==1) return(r*k)
  r*(bt^k-1)/(bt-1)
}

fsechMyerson_rba <- function(p,r,b,alpha,log){
  bt <- r/b
  qn1ma <-2/pi*log(tan(pi/2*(1-alpha))) #sech QF
  k <- 2/pi*log(tan(pi/2*p))/qn1ma
  if(bt==1){
    res <- r*pi*(1+tan(pi*p/2)^2)/(pi*qn1ma*tan(pi*p/2))
  } else {
    res <- r*log(bt)*(bt^k)*pi*(1+tan(pi*p/2)^2)/((bt-1)*pi*qn1ma*tan(pi*p/2))
  }
  if(log) return(log(res))
  res
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
#' @param log logical; if TRUE, log density is returned
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
  b <- (q2-q1)
  q2+ssechMyerson_rba(p,r,b,alpha)
}

#' @rdname sechMyerson
#' @export
#'
#' @examples
#' fsechMyerson(0.25, 10, 20, 40)
fsechMyerson <- function(p, q1,q2,q3, alpha=0.25, log=FALSE){
  r <- (q3-q2)
  b <- (q2-q1)
  qn1ma <-2/pi*log(tan(pi/2*(1-alpha)))
  k <- 2/pi*log(tan(pi/2*p))/qn1ma
  fsechMyerson_rba(p,r,b,alpha,log)
}

#' @rdname sechMyerson
#' @export
#'
#' @examples
#' dqsechMyerson(0.25, 10, 20, 40)
dqsechMyerson <- function(p, q1,q2,q3, alpha=0.25, log=FALSE){
  res <- fsechMyerson(p, q1,q2,q3, alpha, log=FALSE)
  if(log) return(ifelse(is.finite(res),-log(res),res))
  1/res
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
  b <- (q2-q1)
  bt <- r/b
  qn1ma <-2/pi*log(tan(pi/2*(1-alpha)))
  if(bt==1)
    return(2/pi*atan(exp(pi/2*(qn1ma*(q-q2)/r))))

  nom <- log1p((q-q2)*(bt-1)/r)
  den <- log(bt)
  psi <- qn1ma*(nom/den)
  2/pi*atan(exp(pi/2*psi))
}


#' @rdname sechMyerson
#' @export
#'
#' @examples
#' dsechMyerson(15, 10, 20, 40)
dsechMyerson <- function(x, q1,q2,q3, alpha=0.25, log=FALSE){
  # Probability density function
  r <- (q3-q2)
  b <- (q2-q1)
  bt <- r/b
  qn1ma <-2/pi*log(tan(pi/2*(1-alpha)))
  psi <- qn1ma*(x-q2)/r
  if(bt==1) {
    res <- 1/2*sech(pi/2*psi)
  } else {
    psi <- qn1ma*(log1p((x-q2)*(bt-1)/r)/log(bt))
    num <- qn1ma*(bt-1)
    den <- (r+(x-q2)*(bt-1))*log(bt)
    res <- num/den*1/2*sech(pi/2*psi)
  }
  if(log) return(log(res))
  res
}







