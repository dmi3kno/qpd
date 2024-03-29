########################### MYERSON ###################################
#' Quantile function for Modified Myerson distribution.
#'
#' @param p vector of probabilities
#' @param q1,q2,q3 numeric values representing bottom quantile, middle quantile (50th percentile)
#' and top quantile. Quantiles are assumed to be symmetrical.
#' @param alpha numerical fixed proportion of distribution below the bottom quantile.
#' Default value is 0.25
#' @param sfun function; standard vectorized quantile function with a single depth parameter
#' @param log,log.p logical; if TRUE, probabilities p are given as log(p)
#' @param lower.tail logical; if TRUE (default), probabilities are \eqn{P[X\leq x]}, otherwise, \eqn{P[X > x]}
#'
#' @return a vector of quantiles of length equal to `length(x)`
#' @name MMyerson
#' @export
#' @importFrom stats qnorm
#'
#' @examples
#' pgrd <- make_pgrid()
#' qMMyerson(pgrd, 10, 20, 40, alpha=0.1, sfun=stats::qlogis)

qMMyerson <- function(p, q1,q2,q3, alpha=0.25, sfun=stats::qnorm, lower.tail=TRUE, log.p=FALSE){
  if(log.p) p <- exp(p)
  if(!lower.tail) p <- 1-p
  rho <- (q3-q2)
  b <- rho/(q2-q1)
  k <- sfun(p)/sfun(1-alpha)
  if(b==1) return(q2+rho*k)
  q2+rho*(b^k-1)/(b-1)
}

#' @keywords internal
sMyerson_ba <- function(p,b,alpha){
  qn1ma <-stats::qnorm(1-alpha)
  k <- stats::qnorm(p)/qn1ma
  if(b==1) return(k)
  (b^k-1)/(b-1)
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
#' @param log,log.p logical; if TRUE, probabilities p are given as log(p)
#' @param lower.tail logical; if TRUE (default), probabilities are \eqn{P[X\leq x]}, otherwise, \eqn{P[X > x]}
#'
#'
#' @return a vector of probabilities of length equal to `length(x)`
#' @name Myerson
#' @export
#' @importFrom stats qnorm dnorm pnorm rnorm runif
#'
#' @examples
#' qMyerson(0.25, 10, 20, 40)
qMyerson <- function(p, q1,q2,q3, alpha=0.25, lower.tail=TRUE, log.p=FALSE){
  if(log.p) p <- exp(p)
  if(!lower.tail) p <- 1-p
  rho <- (q3-q2)
  b <- rho/(q2-q1)
  q2+rho*sMyerson_ba(p,b,alpha)
}

#' @rdname Myerson
#' @export
#'
#' @examples
#' fMyerson(0.25, 10, 20, 40)
fMyerson <- function(p, q1,q2,q3, alpha=0.25, log=FALSE, lower.tail=TRUE, log.p=FALSE){
  if(log.p) p <- exp(p)
  if(!lower.tail) p <- 1-p
  rho <- (q3-q2)
  b <- rho/(q2-q1)
  qn1ma <-stats::qnorm(1-alpha)
  k <- stats::qnorm(p)/qn1ma
  if(b==1){
    res <- rho*fnorm(p)/qn1ma
  } else {
    res <- rho*log(b)*(b^k)*fnorm(p)/(b-1)/qn1ma
  }
  if(log) return(log(res))
  res
}

#' @rdname Myerson
#' @export
#'
#' @examples
#' fMyerson(0.25, 10, 20, 40)
dqMyerson <- function(p, q1,q2,q3, alpha=0.25, log=FALSE, lower.tail=TRUE, log.p=FALSE){
  if(log.p) p <- exp(p)
  if(!lower.tail) p <- 1-p
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
  rho <- (q3-q2)
  b <- rho/(q2-q1)
  qn1ma <-stats::qnorm(1-alpha)
  if(b==1){
    res <- stats::dnorm(x, mean=q2, sd=rho/qn1ma)
  }else{
    psi <- qn1ma*(log1p((x-q2)*(b-1)/rho)/log(b))
    num <- qn1ma*(b-1)
    den <- (rho+(x-q2)*(b-1))*log(b)
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
pMyerson <- function(q, q1,q2,q3, alpha=0.25, lower.tail=TRUE, log.p=FALSE){
  # cumulative density
  rho <- (q3-q2)
  b <- rho/(q2-q1)
  qn1ma <-stats::qnorm(1-alpha)

  if(b==1){
    p <- stats::pnorm(q, mean=q2, sd=rho/qn1ma)
    if(!lower.tail) p <- 1-p
    if(log.p) p <- log(p)
    return(p)
  }

  nom <- log1p((q-q2)*(b-1)/rho)
  den <- log(b)
  psi <- qn1ma*(nom/den)
  p <- stats::pnorm(psi)
  if(!lower.tail) p <- 1-p
  if(log.p) p <- log(p)
  return(p)
}

########################### LOGIT MYERSON ###################################

slogitMyerson_ba <- function(p,b,alpha){
  qn1ma <-logit(1-alpha)
  k <- logit(p)/qn1ma
  if(b==1) return(k)
  (b^k-1)/(b-1)
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
#' @param log,log.p logical; if TRUE, probabilities p are given as log(p)
#' @param lower.tail logical; if TRUE (default), probabilities are \eqn{P[X\leq x]}, otherwise, \eqn{P[X > x]}
#'
#'
#' @return a vector of probabilities of length equal to `length(x)`
#' @name logitMyerson
#' @export
#' @importFrom stats qnorm dnorm pnorm rnorm runif
#'
#' @examples
#' qlogitMyerson(0.25, 10, 20, 40)
qlogitMyerson <- function(p, q1,q2,q3, alpha=0.25, lower.tail=TRUE, log.p=FALSE){
  if(log.p) p <- exp(p)
  if(!lower.tail) p <- 1-p
  rho <- (q3-q2)
  b <- rho/(q2-q1)
  q2+rho*slogitMyerson_ba(p,b,alpha)
}

#' @rdname logitMyerson
#' @export
#'
#' @examples
#' flogitMyerson(0.25, 10, 20, 40)
flogitMyerson <- function(p, q1,q2,q3, alpha=0.25, log=FALSE, lower.tail=TRUE, log.p=FALSE){
  if(log.p) p <- exp(p)
  if(!lower.tail) p <- 1-p
  rho <- (q3-q2)
  b <- rho/(q2-q1)
  qn1ma <-logit(1-alpha)
  k <- logit(p)/qn1ma
  if(b==1){
    res <- rho*1/(p*(1-p)*qn1ma)
  } else {
    res <- rho*log(b)*(b^k)/(p*(1-p)*(b-1)*qn1ma)
  }

  if(log) return(log(res))
  res
}

#' @rdname logitMyerson
#' @export
#'
#' @examples
#' dqlogitMyerson(0.25, 10, 20, 40)
dqlogitMyerson <- function(p, q1,q2,q3, alpha=0.25, log=FALSE, lower.tail=TRUE, log.p=FALSE){
  if(log.p) p <- exp(p)
  if(!lower.tail) p <- 1-p
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
rlogitMyerson <- function(n,q1,q2,q3,alpha=0.25, lower.tail=TRUE, log.p=FALSE){
  # random number generator
  if(length(n)>1) n <-length(n)
  qlogitMyerson(stats::runif(n), q1, q2, q3, alpha)
}


#' @rdname logitMyerson
#'
#' @return a vector of exceedance probabilities of length equal to `length(x)`.
#' @export
#' @importFrom stats qnorm dnorm pnorm rnorm runif plogis qlogis
#'
#' @examples
#' plogitMyerson(0.25, 10, 20, 40)
plogitMyerson <- function(q, q1,q2,q3, alpha=0.25, lower.tail=TRUE, log.p=FALSE){
  # cumulative density
  rho <- (q3-q2)
  b <- rho/(q2-q1)
  qn1ma <- stats::qlogis(1-alpha)
  if(b==1){
      p <- stats::plogis(qn1ma*(q-q2)/rho)
      if(!lower.tail) p <- 1-p
      if(log.p) p <- log(p)
      return(p)
      }

  nom <- log1p((q-q2)*(b-1)/rho)
  den <- log(b)
  psi <- qn1ma*(nom/den)
  p <- stats::plogis(psi)
  if(!lower.tail) p <- 1-p
  if(log.p) p <- log(p)
  return(p)
}


#' @rdname logitMyerson
#' @export
#'
#' @examples
#' dlogitMyerson(15, 10, 20, 40)
dlogitMyerson <- function(x, q1,q2,q3, alpha=0.25, log=FALSE){
  # Probability density function
  rho <- (q3-q2)
  b <- rho/(q2-q1)
  qn1ma <-logit(1-alpha)
  psi <- qn1ma*(x-q2)/rho
  if(b==1){
    res <- stats::dlogis(psi)
  }else{
    psi <- qn1ma*(log1p((x-q2)*(b-1)/rho)/log(b))
    num <- qn1ma*(b-1)
    den <- (rho+(x-q2)*(b-1))*log(b)
    res <- num/den*stats::dlogis(psi)
  }

  if(log) return(log(res))
  res
}

####################### CAUCHY MYERSON ###################################
scauchyMyerson_ba <- function(p,b,alpha){
  qn1ma <-tan(pi*(0.5-alpha))
  k <- tan(pi*(p-0.5))/qn1ma
  if(b==1) return(k)
  (b^k-1)/(b-1)
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
#' @param log,log.p logical; if TRUE, probabilities p are given as log(p)
#' @param lower.tail logical; if TRUE (default), probabilities are \eqn{P[X\leq x]}, otherwise, \eqn{P[X > x]}
#'
#'
#' @return a vector of probabilities of length equal to `length(x)`
#' @name cauchyMyerson
#' @export
#' @importFrom stats qnorm dnorm pnorm rnorm runif
#'
#' @examples
#' qcauchyMyerson(0.25, 10, 20, 40)
qcauchyMyerson <- function(p, q1,q2,q3, alpha=0.25, lower.tail=TRUE, log.p=FALSE){
  if(log.p) p <- exp(p)
  if(!lower.tail) p <- 1-p
  rho <- (q3-q2)
  b <- rho/(q2-q1)
  q2+rho*scauchyMyerson_ba(p,b,alpha)
}

#' @rdname cauchyMyerson
#' @export
#'
#' @examples
#' fcauchyMyerson(0.25, 10, 20, 40)
fcauchyMyerson <- function(p, q1,q2,q3, alpha=0.25, log=FALSE, lower.tail=TRUE, log.p=FALSE){
  if(log.p) p <- exp(p)
  if(!lower.tail) p <- 1-p
  rho <- (q3-q2)
  b <- rho/(q2-q1)
  qn1ma <-tan(pi*(0.5-alpha))
  k <- tan(pi*(p-0.5))/qn1ma
  if(b==1){
    res <- rho*pi*(1+tan(pi*(p-0.5))^2)/qn1ma
  } else {
    res <- rho*log(b)*(b^k)*pi*(1+tan(pi*(p-0.5))^2)/((b-1)*qn1ma)
  }

  if(log) return(log(res))
  res
}

#' @rdname cauchyMyerson
#' @export
#'
#' @examples
#' dqcauchyMyerson(0.25, 10, 20, 40)
dqcauchyMyerson <- function(p, q1,q2,q3, alpha=0.25, log=FALSE, lower.tail=TRUE, log.p=FALSE){
  if(log.p) p <- exp(p)
  if(!lower.tail) p <- 1-p
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
pcauchyMyerson <- function(q, q1,q2,q3, alpha=0.25, lower.tail=TRUE, log.p=FALSE){
  # cumulative density
  rho <- (q3-q2)
  b <- rho/(q2-q1)
  qn1ma <-tan(pi*(0.5-alpha))
  if(b==1){
      p <- stats::pcauchy(qn1ma*(q-q2)/rho)
      if(!lower.tail) p <- 1-p
      if(log.p) p <- log(p)
      return(p)
      }

  nom <- log1p((q-q2)*(b-1)/rho)
  den <- log(b)
  psi <- qn1ma*(nom/den)
  p  <- stats::pcauchy(psi)
  if(!lower.tail) p <- 1-p
  if(log.p) p <- log(p)
  p
}


#' @rdname cauchyMyerson
#' @export
#'
#' @examples
#' dcauchyMyerson(15, 10, 20, 40)
dcauchyMyerson <- function(x, q1,q2,q3, alpha=0.25, log=FALSE){
  # Probability density function
  rho <- (q3-q2)
  b <- rho/(q2-q1)
  qn1ma <-tan(pi*(0.5-alpha))
  psi <- qn1ma*(x-q2)/rho
  if(b==1){
    res <- stats::dcauchy(psi)
  }else{
    psi <- qn1ma*(log1p((x-q2)*(b-1)/rho)/log(b))
    num <- qn1ma*(b-1)
    den <- (rho+(x-q2)*(b-1))*log(b)
    res <- num/den*stats::dcauchy(psi)
  }

  if(log) return(log(res))
  res
}

######################### SECH MYERSON #######################################

ssechMyerson_ba <- function(p,b,alpha){
  qn1ma <-2/pi*log(tan(pi/2*(1-alpha))) #sech QF
  k <- 2/pi*log(tan(pi/2*p))/qn1ma
  if(b==1) return(k)
  (b^k-1)/(b-1)
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
#' @param log,log.p logical; if TRUE, probabilities p are given as log(p)
#' @param lower.tail logical; if TRUE (default), probabilities are \eqn{P[X\leq x]}, otherwise, \eqn{P[X > x]}
#'
#'
#' @return a vector of probabilities of length equal to `length(x)`
#' @name sechMyerson
#' @export
#' @importFrom stats qnorm dnorm pnorm rnorm runif
#'
#' @examples
#' qsechMyerson(0.25, 10, 20, 40)
qsechMyerson <- function(p, q1,q2,q3, alpha=0.25, lower.tail=TRUE, log.p=FALSE){
  if(log.p) p <- exp(p)
  if(!lower.tail) p <- 1-p
  rho <- (q3-q2)
  b <- rho/(q2-q1)
  q2+rho*ssechMyerson_ba(p,b,alpha)
}

#' @rdname sechMyerson
#' @export
#'
#' @examples
#' fsechMyerson(0.25, 10, 20, 40)
fsechMyerson <- function(p, q1,q2,q3, alpha=0.25, log=FALSE, lower.tail=TRUE, log.p=FALSE){
  if(log.p) p <- exp(p)
  if(!lower.tail) p <- 1-p
  rho <- (q3-q2)
  b <- rho/(q2-q1)
  qn1ma <-qsech(1-alpha) #sech QF
  k <- qsech(p)/qn1ma
  if(b==1){
    res <- rho*fsech(p)/qn1ma
  } else {
    res <- rho*log(b)*(b^k)*fsech(p)/(b-1)/qn1ma
  }
  if(log) return(log(res))
  res
}

#' @rdname sechMyerson
#' @export
#'
#' @examples
#' dqsechMyerson(0.25, 10, 20, 40)
dqsechMyerson <- function(p, q1,q2,q3, alpha=0.25, log=FALSE, lower.tail=TRUE, log.p=FALSE){
  if(log.p) p <- exp(p)
  if(!lower.tail) p <- 1-p
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
psechMyerson <- function(q, q1,q2,q3, alpha=0.25, lower.tail=TRUE, log.p=FALSE){
  # cumulative density
  rho <- (q3-q2)
  b <- rho/(q2-q1)
  qn1ma <-qsech(1-alpha)
  if(b==1){
      p <- 2/pi*atan(exp(pi/2*(qn1ma*(q-q2)/rho)))
      if(!lower.tail) p <- 1-p
      if(log.p) p <- log(p)
      return(p)
  }

  nom <- log1p((q-q2)*(b-1)/rho)
  den <- log(b)
  psi <- qn1ma*(nom/den)
  p <- 2/pi*atan(exp(pi/2*psi))
  if(!lower.tail) p <- 1-p
  if(log.p) p <- log(p)
  p
}


#' @rdname sechMyerson
#' @export
#'
#' @examples
#' dsechMyerson(15, 10, 20, 40)
dsechMyerson <- function(x, q1,q2,q3, alpha=0.25, log=FALSE){
  # Probability density function
  rho <- (q3-q2)
  b <- rho/(q2-q1)
  qn1ma <-qsech(1-alpha)
  psi <- qn1ma*(x-q2)/rho
  if(b==1) {
    res <- 1/2*sech(pi/2*psi)
  } else {
    psi <- qn1ma*(log1p((x-q2)*(b-1)/rho)/log(b))
    num <- qn1ma*(b-1)
    den <- (rho+(x-q2)*(b-1))*log(b)
    res <- num/den*1/2*sech(pi/2*psi)
  }
  if(log) return(log(res))
  res
}

######################### TUKEY MYERSON #######################################
# Q(u)= 1/lambda*(u^lambda-(1-p)^lambda) ,for lamba!=0
# Q(u)=logit(p) for lambda==0
stukeyMyerson_ba <- function(p,b, alpha, tlambda){
  qn1ma <-qtlambda(1-alpha, tlambda)
  k <- qtlambda(p, tlambda)/qn1ma # here the scale parameters 1/lambda cancel out
  if(b==1) return(k)
  (b^k-1)/(b-1)
}

#' Density function for Tukey-Myerson distribution (Myerson distribution with Tukey lambda base function).
#'
#' @param q vector of quantiles.
#' @param p vector of probabilities
#' @param n number of observations. If `length(n)>1``, the length is taken to be the number required.
#' @param q1,q2,q3 numeric values representing bottom quantile, middle quantile (50th percentile)
#' and top quantile. Quantiles are assumed to be symmetrical.
#' @param alpha numerical fixed proportion of distribution below the bottom quantile.
#' Default value is 0.25
#' @param tlambda numerical Tukey lambda `lambda` parameter
#' @param log,log.p logical; if TRUE, probabilities p are given as log(p)
#' @param lower.tail logical; if TRUE (default), probabilities are \eqn{P[X\leq x]}, otherwise, \eqn{P[X > x]}
#'
#'
#' @return a vector of probabilities of length equal to `length(x)`
#' @name tukeyMyerson
#' @export
#' @importFrom stats qnorm dnorm pnorm rnorm runif
#'
#' @examples
#' qtukeyMyerson(0.25, 10, 20, 40)
qtukeyMyerson <- function(p, q1,q2,q3, alpha=0.25, tlambda=0, lower.tail=TRUE, log.p=FALSE){
  if(log.p) p <- exp(p)
  if(!lower.tail) p <- 1-p
  if(tlambda==0) return(qlogitMyerson(p,q1,q2,q3,alpha))
  rho <- (q3-q2)
  b <- rho/(q2-q1)
  q2+rho*stukeyMyerson_ba(p,b,alpha, tlambda)
}

#' @rdname tukeyMyerson
#' @export
#'
#' @examples
#' flogitMyerson(0.25, 10, 20, 40)
ftukeyMyerson <- function(p, q1,q2,q3, alpha=0.25, tlambda=0, log=FALSE, lower.tail=TRUE, log.p=FALSE){
  if(log.p) p <- exp(p)
  if(!lower.tail) p <- 1-p
  if(tlambda==0) return(flogitMyerson(p,q1=q1,q2=q2,q3=q3,alpha=alpha, log=log))
  rho <- (q3-q2)
  b <- rho/(q2-q1)
  qn1ma <-qtlambda(1-alpha, tlambda)
  k <- qtlambda(p, tlambda)/qn1ma # here the scale parameters 1/lambda cancel out
  if(b==1){
    res <- rho*ftlambda(p, tlambda)/qn1ma
  } else {
    res <- rho*log(b)*(b^k)*ftlambda(p, tlambda)/(b-1)/qn1ma
  }

  if(log) return(log(res))
  res
}

#' @rdname tukeyMyerson
#' @export
#'
#' @examples
#' dqlogitMyerson(0.25, 10, 20, 40)
dqtukeyMyerson <- function(p, q1,q2,q3, alpha=0.25, tlambda=0, log=FALSE, lower.tail=TRUE, log.p=FALSE){
  if(log.p) p <- exp(p)
  if(!lower.tail) p <- 1-p
  res <- ftukeyMyerson(p, q1,q2,q3, alpha, tlambda=tlambda, log=FALSE)
  if(log) return(ifelse(is.finite(res),-log(res),res))
  1/res
}

#' @rdname tukeyMyerson
#' @export
#' @importFrom stats qnorm dnorm pnorm rnorm runif
#'
#' @examples
#' rlogitMyerson(1, 10, 20, 40)
rtukeyMyerson <- function(n,q1,q2,q3,alpha=0.25, tlambda=0){
  # random number generator
  if(length(n)>1) n <-length(n)
  qtukeyMyerson(stats::runif(n), q1, q2, q3, alpha, tlambda)
}

#' @param ... used by method
#' @param lower,upper the `stats::uniroot` lower and upper end points of the interval to be searched. Defaults are 0 and 1, respectively
#' @param tol the `stats::uniroot` desired accuracy (convergence tolerance). Default value 1e-06
#' @param silent the `base::try` argument. Default is TRUE
#' @param trace integer number passed to `stats::uniroot`; if positive, tracing information is produced. Higher values giving more details.
#' @include iqf.R
#' @rdname tukeyMyerson
#'
#' @return a vector of exceedance probabilities of length equal to `length(x)`.
#' @export
#' @importFrom stats uniroot
#'
#' @examples
#' ptukeyMyerson(0.25, 10, 20, 40)
ptukeyMyerson <- iqf(qtukeyMyerson)

