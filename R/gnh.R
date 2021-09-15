#' Functions for g-and-h distribution
#'
#' @param p numeric vector of probabilities
#' @param A numeric location parameter of g-and-h distribution
#' @param B positive numeric scale parameter of g-and-h distribution
#' @param C numeric parameter of g-and-h distribution. Default is 0.8
#' @param g numeric skeweness parameter of g-and-h distribution
#' @param h non-negative numeric kurtosis parameter of g-and-h distribution
#' @rdname gnh
#' @param zscale is the probability input provided on scale of `z` values, i.e. as `z=qnorm(p)`. Default is FALSE
#'
#' @return vector of quantiles for g-and-h distribution
#' @export
#'
#' @examples
#' qgnh(0.1, 3,1,0.8, 2,0.5)
qgnh <- function(p,A,B,C=0.8,g,h, zscale=FALSE){
  stopifnot(B>0 && h>=0)
  if(zscale){z <-p} else {z <- qnorm(p)}
  res <- A+B*(1+C*tanh(g*z/2))*z*exp(h*z^2/2)
  res
}

#' @param log logical should the result be returned as log(). Default is FALSE
#'
#' @return vector of quantiles for g-and-h distribution
#' @rdname gnh
#' @export
fgnh <- function(p,A,B,C=0.8,g,h, log=FALSE, zscale=FALSE){
  # from gk paper
  stopifnot(B>0 && h>=0)
  if(zscale){z <-p} else {z <- qnorm(p)}
  dQdz <- B*exp(h/2*z^2)*((1+C*tanh(g*z/2))*(1+h*z^2)+C*g*z/(2*cosh(g*z/2)^2))
  res <- dQdz*fnorm(p)
  if(log) return(log(res))
  res
}

#' @rdname gnh
#' @export
dqgnh <- function(p,A,B,C=0.8,g,h, log=FALSE, zscale=FALSE){
  stopifnot(B>0 && h>=0)
  res <- fgnh(p, A,B, C, g, h, log=FALSE, zscale)
  if (log) return(log(1/res))
  1/res
}
#' @param n numeric number of samples to draw
#' @rdname gnh
#' @export
#' @importFrom stats runif
rgnh <- function(n,A,B,C=0.8,g,h){
  qgnh(stats::runif(n), A,B,C,g,h)
}
#' @param q quantile value for which the corresponding cumulative probability value should be found.
#' @param tol numeric tolerance balue for approximating CDF. Default 1e-15
#' @param maxiter numeric maximum number of iteration
#' @param n_grid integer size of helper grid to be passed to `make_pgrid`. Default is 50
#' @param s_grid integer beta shape of helper grid to be passed to `make_pgrid`. Default is 5
#' @export
#' @rdname gnh
#' @importFrom stats pnorm qnorm uniroot
pgnh <- function(q, A,B,C=0.8, g,h, n_grid=100L, s_grid=5L, tol=1e-15, maxiter=1e3){
  afun <- function(z, q, A,B, C, g, h) {q - qgnh(z,A,B,C,g,h, zscale=TRUE)}
  stopifnot(B>0 && h>=0)
  mtol <- .Machine$double.eps
  p_grd <- make_pgrid(n=n_grid, s=s_grid)
  q_grd <- qgnh(p_grd, A,B,C, g,h, zscale=FALSE)
  idx_lower <- findInterval(q, q_grd, all.inside = TRUE)
  idx_upper <- idx_lower+1L
  int_lower <- stats::qnorm(p_grd[idx_lower])
  int_upper <- stats::qnorm(p_grd[idx_upper])

  zs <- mapply(function(.q, .il, .iu) {
    tmp_zs <- NULL
    tmp_zs <- stats::uniroot(afun, q=.q, A=A, B=B, C=C, g=g, h=h,
                             interval=c(.il,.iu), extendInt="downX",
                             check.conv=TRUE, tol = tol, maxiter = maxiter)
    if(is.null(tmp_zs)) res <- NA else res <- tmp_zs$root
    res
  }, q, int_lower, int_upper)

  ps <- stats::pnorm(zs)

  ps[!is.finite(ps)] <- NA_real_
  ps
}

#' @export
#' @rdname gnh
#' @examples
#' is_gnh_valid(A=5, B=5, C=0.8, g=0.5, h=0.5)
is_gnh_valid <- function(A, B, C=0.8, g, h, n_grid=100L, s_grid=2L){
  grd <- make_pgrid(n_grid, s_grid)
  all(fgnh(grd, A, B, C, g, h)>0)
}

