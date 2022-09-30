#' Quantile function of Johnson Distribution System (JDS) SU distribution
#' The transformation function is hyperbolic sine transform
#' @keywords internal
qJDS_SU <- function(p, zeta, lambda, sigma, gamma){
  zeta+lambda*sinh(sigma*(stats::qnorm(p)+gamma))
}

#' Quantile function of Johnson Distribution System (JDS) SB distribution
#' The transformation function is inverse logit T(y)=exp(y)/(1+exp(y))
#' @keywords internal
qJDS_SB <- function(p, zeta, lambda, sigma, gamma){
  num <- lambda*exp(sigma*(stats::qnorm(p)+gamma))
  denom <- 1 + exp(sigma*(stats::qnorm(p)+gamma))
  zeta+num/denom
}

#' The J-QPD-B distribution
#'
#' Density, distribution function, quantile function and random generation for the
#' Johnson Quantile Parametrized (Bounded) distribution parametrized by three symmetrical
#' quantiles (`q1`, `q2` and `q3`) and an `alpha` argument, representing the proportion of
#' density below the bottom quantile. All functions are vectorized.
#'
#' @details The distribution is created by applying the inverse probit transform T(x)=Phi(x) to the Johnson SU distribution
#'
#' @param x,q vector of quantiles
#' @param p vector of probabilities
#' @param q1,q2,q3 vectors of values, corresponding to lower, median and upper (symmetrical) quantiles
#' @param lower,upper vectors of lower and upper bounds of distribution
#' @param alpha vector of proportions of probability density under the lower bound
#' (or above the upper bound, since the quantiles are symmetrical)
#'
#' @return vector of values
#' @name JQPDB
#' @export
#' @importFrom stats qnorm dnorm pnorm rnorm
#'
#' @examples
#' # should result in c(0,1,5,12,20)
#' qJQPDB(c(0, 0.05, 0.5, 0.95, 1), 1, 5, 12, 0, 20, alpha=0.05)
qJQPDB <- function(p, q1, q2, q3, lower, upper, alpha=0.1){
  stopifnot(all(lower < upper))
  stopifnot(all(q1<q2 & q2<q3))
  stopifnot(alpha>0 & alpha<0.5)
  small_c <- stats::qnorm(1-alpha)

  u_m_l <- upper-lower

  L <- stats::qnorm((q1-lower)/u_m_l)
  B <- stats::qnorm((q2-lower)/u_m_l)
  H <- stats::qnorm((q3-lower)/u_m_l)
  sgn <- sign(L+H-2*B)

  if (sgn==0)
      return(lower + u_m_l*stats::pnorm(B+0.5*(H-L)/small_c*stats::qnorm(p)))

  sigma <- (1/small_c)*acosh(0.5*(H-L)/pmin(B-L, H-B))
  lambda <- (H-L)/sinh(2*sigma*small_c)
  zeta  <- 0.5*(L*(1+sgn)+H*(1-sgn)) # idea from Matlab code in Hadlock thesis

  lower + u_m_l*stats::pnorm(qJDS_SU(p, zeta, lambda, sigma, gamma=sgn*small_c))
}


#' @rdname JQPDB
#'
#' @export
#' @importFrom stats qnorm dnorm pnorm rnorm
#'
#' @examples
#' # should return c(0.00, 0.05, 0.50, 0.95, 1.00)
#' pJQPDB(c(0, 1, 5, 12, 20), 1,5,12, 0, 20, alpha=0.05)
pJQPDB <- function(q, q1, q2, q3, lower, upper, alpha=0.1){
  stopifnot(all(lower < upper))
  stopifnot(all(q1<q2 & q2<q3))
  stopifnot(alpha>0 & alpha<0.5)
  small_c <- stats::qnorm(1-alpha)

  u_m_l <- upper-lower

  L <- stats::qnorm((q1-lower)/u_m_l)
  B <- stats::qnorm((q2-lower)/u_m_l)
  H <- stats::qnorm((q3-lower)/u_m_l)
  sgn <- sign(L+H-2*B)

  q_xmlb_uml <- stats::qnorm((q-lower)/u_m_l)

  if (sgn==0)
    return(stats::pnorm((2*small_c/(H-L))*(-B+q_xmlb_uml)) )

  sigma <- (1/small_c)*acosh(0.5*(H-L)/pmin(B-L, H-B))
  lambda <- (H-L)/sinh(2*sigma*small_c)
  zeta  <- 0.5*(L*(1+sgn)+H*(1-sgn)) # idea from Matlab code in Hadlock thesis

  stats::pnorm((1/sigma)*asinh((1/lambda)*(q_xmlb_uml-zeta))-sgn*small_c)
}

#' @rdname JQPDB
#'
#' @export
#' @importFrom stats qnorm dnorm pnorm rnorm
#' @examples
#' # should return vector with first and last element equal to NaN
#' dJQPDB(c(0, 1, 5, 12, 20), 1,5,12, 0, 20, alpha=0.05)
dJQPDB <- function(x, q1, q2, q3, lower, upper, alpha=0.1){
  stopifnot(all(lower < upper))
  stopifnot(all(q1<q2 & q2<q3))
  stopifnot(alpha>0 & alpha<0.5)
  small_c <- stats::qnorm(1-alpha)

  u_m_l <- upper-lower

  L <- stats::qnorm((q1-lower)/u_m_l)
  B <- stats::qnorm((q2-lower)/u_m_l)
  H <- stats::qnorm((q3-lower)/u_m_l)
  sgn <- sign(L+H-2*B)

  q_xmlb_uml <- stats::qnorm((x-lower)/u_m_l)

  if (sgn==0)
      return(2*small_c/((H-L)*u_m_l)*
               (1/stats::dnorm(q_xmlb_uml))*
               stats::dnorm(2*small_c/(H-L)*(-B+q_xmlb_uml))
             )

  sigma <- (1/small_c)*acosh(0.5*(L-H)/pmax(B-H, L-B))
  lambda <- (H-L)/sinh(2*sigma*small_c)
  zeta  <- 0.5*(L*(1+sgn)+H*(1-sgn)) # idea from Matlab code in Hadlock thesis

  (1/sigma)*
    (1/u_m_l)*
    stats::dnorm(-sgn*small_c+(1/sigma)*asinh((1/lambda)*(-zeta+q_xmlb_uml))) *
    (1/stats::dnorm(q_xmlb_uml))/
    sqrt((lambda^2)+((-zeta+q_xmlb_uml))^2)
}


#' The J-QPD-S distribution
#'
#' Density, distribution function, quantile function and random generation for the
#' Johnson Quantile Parametrized (Semibounded) distribution parametrized by three symmetrical
#' quantiles (`q1`, `q2` and `q3`) and an `alpha` argument, representing the proportion of
#' density below the bottom quantile. All functions are vectorized.
#'
#' @details The distribution is created by applying the exponential transform T(x)=exp(x) to the Johnson SU distribution
#'
#' @param q  vector of quantiles
#' @param p vector of probabilities
#' @param q1,q2,q3 vectors of values, corresponding to lower, median and upper (symmetrical) quantiles
#' @param lower vector of lower bounds of distribution
#' @param alpha vector of proportions of probability density under the lower bound
#' (or above the upper bound, since the quantiles are symmetrical)
#'
#' @return vector of values
#' @name JQPDS
#' @export
#' @importFrom stats qlnorm plnorm dlnorm qnorm dnorm pnorm rnorm
#'
#' @examples
#' # should result in c(0,1,5,12,20)
#' qJQPDB(c(0, 0.05, 0.5, 0.95, 1), 1, 5, 12, 0, 20, alpha=0.05)
#' qJQPDS(0.6, 20,50,90)
qJQPDS <- function(p, q1, q2, q3, lower=0, alpha=0.1){
  stopifnot(all(q1<q2 & q2<q3))
  stopifnot(alpha>0 & alpha<0.5)
  small_c <- stats::qnorm(1-alpha)

  L <- log(q1-lower)
  B <- log(q2-lower)
  H <- log(q3-lower)
  sgn <- sign(L+H-2*B)

  theta  <- 0.5*(exp(L)*(1+sgn)+exp(H)*(1-sgn)) # idea from Matlab code in Hadlock thesis

  if (sgn==0)
    return(stats::qlnorm(p,meanlog=log(theta), sdlog=(H-B)/small_c))

  sigma <- (1/small_c)*sinh(acosh(0.5*(H-L)/pmin(B-L, H-B)))
  lambda <- 1/(sigma*small_c)*pmin(H-B, B-L)

  lower + theta*exp(lambda*sinh(asinh(sigma*stats::qnorm(p))+asinh(sgn*small_c*sigma)))
}


#' @rdname JQPDS
#'
#' @export
#' @importFrom stats dlnorm plnorm qnorm dnorm pnorm rnorm
#'
#' @examples
#' # should return c(0.00, 0.05, 0.50, 0.95, 1.00)
#' pJQPDB(c(0, 1, 5, 12, 20), 1,5,12, 0, 20, alpha=0.05)
pJQPDS<- function(q, q1, q2, q3, lower=0, alpha=0.1){
  stopifnot(all(q1<q2 & q2<q3))
  stopifnot(alpha>0 & alpha<0.5)
  small_c <- stats::qnorm(1-alpha)

  L <- log(q1-lower)
  B <- log(q2-lower)
  H <- log(q3-lower)
  sgn <- sign(L+H-2*B)

  theta  <- 0.5*(exp(L)*(1+sgn)+exp(H)*(1-sgn)) # idea from Matlab code in Hadlock thesis

  if (sgn==0)
    return(stats::plnorm(q-lower,meanlog=log(theta), sdlog=(H-B)/small_c))

  sigma <- (1/small_c)*sinh(acosh(0.5*(H-L)/pmin(B-L, H-B)))
  lambda <- 1/(sigma*small_c)*pmin(H-B, B-L)

  stats::pnorm((1/sigma)*sinh(asinh((1/lambda)*log((q-lower)/theta))-asinh(sgn*small_c*sigma))) #CHECK
}


#' @param x  vector of observations
#' @rdname JQPDS
#' @export
#' @importFrom stats dlnorm plnorm qnorm dnorm pnorm rnorm
#' @examples
#' # should return vector with first and last element equal to NaN
#' dJQPDB(c(0, 1, 5, 12, 20), 1,5,12, 0, 20, alpha=0.05)
dJQPDS <- function(x, q1, q2, q3, lower=0, alpha=0.1){
  stopifnot(all(q1<q2 & q2<q3))
  stopifnot(alpha>0 & alpha<0.5)
  small_c <- stats::qnorm(1-alpha)

  L <- log(q1-lower)
  B <- log(q2-lower)
  H <- log(q3-lower)
  sgn <- sign(L+H-2*B)

  theta  <- 0.5*(exp(L)*(1+sgn)+exp(H)*(1-sgn)) # idea from Matlab code in Hadlock thesis

  if (sgn==0)
    return(stats::dlnorm(x-lower, meanlog=log(theta), sdlog=(H-B)/small_c))

  sigma <- (1/small_c)*sinh(acosh(0.5*(H-L)/pmin(B-L, H-B)))
  lambda <- 1/(sigma*small_c)*pmin(H-B, B-L)

  par0 <- log((x-lower)/theta)
  par1 <- asinh(small_c*sgn*sigma)-asinh(1/lambda*par0)
  stats::dnorm(sinh(par1)/sigma)*
                cosh(par1)/
                ((x-lower)*sigma*lambda*sqrt(1+(par0/lambda)^2))

}



#' The J-QPD-S Type II distribution
#'
#' Density, distribution function, quantile function and random generation for the
#' Johnson Quantile Parametrized (Semibounded) Type II distribution parametrized by three symmetrical
#' quantiles (`q1`, `q2` and `q3`) and an `alpha` argument, representing the proportion of
#' density below the bottom quantile. All functions are vectorized.
#'
#' @details The distribution is created by applying the exponential transform T(x)=exp(x) to the Johnson SU distribution
#'
#' @param q  vector of quantiles
#' @param p vector of probabilities
#' @param q1,q2,q3 vectors of values, corresponding to lower, median and upper (symmetrical) quantiles
#' @param lower vector of lower bounds of distribution
#' @param alpha vector of proportions of probability density under the lower bound
#' (or above the upper bound, since the quantiles are symmetrical)
#'
#' @return vector of values
#' @name JQPDS2
#' @export
#' @importFrom stats qlnorm plnorm dlnorm qnorm dnorm pnorm rnorm
#'
#' @examples
#'
#' qJQPDS2(0.6, 20,50,90)
qJQPDS2 <- function(p, q1, q2, q3, lower=0, alpha=0.1){
  stopifnot(all(q1<q2 & q2<q3))
  stopifnot(alpha>0 & alpha<0.5)
  small_c <- stats::qnorm(1-alpha)

  L <- q1-lower
  B <- q2-lower
  H <- q3-lower
  s1 <- B/H
  s2 <- L/B
  sgn <- sign(log(s2/s1))

  theta  <- 0.5*(L*(1+sgn)+H*(1-sgn)) # idea from Matlab code in Hadlock thesis

  if (sgn==0)
    return(lower + B*(s1^(-stats::qnorm(p)/small_c)))

  sigma <- (1/small_c)*acosh(0.5*log(s1*s2)/log(pmax(s1, s2)))
  lambda <- (-log(pmax(s1, s2)))/sinh(sigma*small_c)

  lower + theta*exp(lambda*sinh(sigma*(stats::qnorm(p)+sgn*small_c)))
}


#' @rdname JQPDS2
#'
#' @export
#' @importFrom stats dlnorm plnorm qnorm dnorm pnorm rnorm
#'
#' @examples
#' pJQPDS2(15, 20,50,90)
pJQPDS2 <- function(q, q1, q2, q3, lower=0, alpha=0.1){
  stopifnot(all(q1<q2 & q2<q3))
  stopifnot(alpha>0 & alpha<0.5)
  small_c <- stats::qnorm(1-alpha)

  L <- q1-lower
  B <- q2-lower
  H <- q3-lower
  s1 <- B/H
  s2 <- L/B
  sgn <- sign(log(s2/s1))

  theta  <- 0.5*(L*(1+sgn)+H*(1-sgn)) # idea from Matlab code in Hadlock thesis

  if (sgn==0)
    return(stats::pnorm(-small_c/log(s1))*log((q-lower)/B))

  sigma <- (1/small_c)*acosh(0.5*log(s1*s2)/log(pmax(s1, s2)))
  lambda <- (-log(pmax(s1, s2)))/sinh(sigma*small_c)

  stats::pnorm((1/sigma)*asinh((1/lambda)*log((q-lower)/theta))-sgn*small_c)
}


#' @param x  vector of observations
#' @rdname JQPDS2
#' @export
#' @importFrom stats dlnorm plnorm qnorm dnorm pnorm rnorm
#' @examples
#' # should return vector with first and last element equal to NaN
#' dJQPDS2(15, 20,50,90)
dJQPDS2 <- function(x, q1, q2, q3, lower=0, alpha=0.1){
  stopifnot(all(q1<q2 & q2<q3))
  stopifnot(alpha>0 & alpha<0.5)
  small_c <- stats::qnorm(1-alpha)

  L <- q1-lower
  B <- q2-lower
  H <- q3-lower
  s1 <- B/H
  s2 <- L/B
  sgn <- sign(log(s2/s1))

  theta  <- 0.5*(L*(1+sgn)+H*(1-sgn)) # idea from Matlab code in Hadlock thesis

  if (sgn==0)
    return(stats::dnorm(-small_c/log(s1)*log((x-lower)/B)))

  sigma <- (1/small_c)*acosh(0.5*log(s1*s2)/log(pmax(s1, s2)))
  lambda <- (-log(pmax(s1, s2)))/sinh(sigma*small_c)

  (1/sigma)*stats::dnorm((1/sigma)*asinh((1/lambda)*log((x-lower)/theta))-sgn*small_c)*
    ((((1/lambda)*log((x-lower)/theta)^2)+1)^(-0.5))*
    (1/((x-lower)*lambda))

}
