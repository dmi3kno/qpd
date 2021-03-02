#' The J-QPD-B distribution
#' 
#' Density, distribution function, quantile function and random generation for the 
#' Johnson Quantile Parametrized (Bounded) distribution parametrized by three symmetrical
#' quantiles (`q1`, `q2` and `q3`) and an `alpha` argument, representing the proportion of
#' density below the bottom quantile. All functions are vectorized.
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
  n <- sign(L+H-2*B)

  if (n==0)
      return(lower + u_m_l*stats::pnorm(B+0.5*(H-L)/small_c*stats::qnorm(p)))

  sigma <- (1/small_c)*acosh(0.5*(H-L)/pmin(B-L, H-B))
  lambda <- (H-L)/sinh(2*sigma*small_c)
  zeta  <- 0.5*(L*(1+n)+H*(1-n)) # idea from Matlab code in Hadlock thesis

  lower + u_m_l*stats::pnorm(zeta+lambda*sinh(sigma*(stats::qnorm(p)+n*small_c)))
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
  n <- sign(L+H-2*B)
  
  q_xmlb_uml <- stats::qnorm((q-lower)/u_m_l)

  if (n==0)
    return(stats::pnorm((2*small_c/(H-L))*(-B+q_xmlb_uml)) )

  sigma <- (1/small_c)*acosh(0.5*(H-L)/pmin(B-L, H-B))
  lambda <- (H-L)/sinh(2*sigma*small_c)
  zeta  <- 0.5*(L*(1+n)+H*(1-n)) # idea from Matlab code in Hadlock thesis

  stats::pnorm((1/sigma)*asinh((1/lambda)*(q_xmlb_uml-zeta))-n*small_c)
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
  n <- sign(L+H-2*B)
  
  q_xmlb_uml <- stats::qnorm((x-lower)/u_m_l)

  if (n==0)
      return(2*small_c/((H-L)*u_m_l)*
               (1/stats::dnorm(q_xmlb_uml))*
               stats::dnorm(2*small_c/(H-L)*(-B+q_xmlb_uml))
             )

  sigma <- (1/small_c)*acosh(0.5*(L-H)/pmax(B-H, L-B))
  lambda <- (H-L)/sinh(2*sigma*small_c)
  zeta  <- 0.5*(L*(1+n)+H*(1-n)) # idea from Matlab code in Hadlock thesis

  (1/sigma)*
    (1/u_m_l)*
    stats::dnorm(-n*small_c+(1/sigma)*asinh((1/lambda)*(-zeta+q_xmlb_uml))) *
    (1/stats::dnorm(q_xmlb_uml))/
    sqrt((lambda^2)+((-zeta+q_xmlb_uml)^2))
}
