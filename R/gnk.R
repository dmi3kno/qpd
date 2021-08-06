#' Functions for g-and-k distribution
#'
#' @param p numeric vector of probabilities
#' @param A numeric location parameter of g-and-k distribution
#' @param B positive numeric scale parameter of g-and-k distribution
#' @param C numeric parameter of g-and-k distribution. Default is 0.8
#' @param g numeric skeweness parameter of g-and-k distribution
#' @param k non-negative numeric kurtosis parameter of g-and-k distribution
#' @param zscale is the probability input provided on scale of `z` values, i.e. as `z=qnorm(p)`. Default is FALSE
#'
#' @return vector of quantiles for g-and-k distribution
#' @export
#' @rdname gnk
#'
#' @examples
#' qgnk(0.1, 3,1,0.8, 2,0.5)
qgnk <- function(p,A,B,C=0.8,g,k, zscale=FALSE){
  stopifnot(B>0 && k>=-0.5)
  if(zscale){z <-p} else {z <- qnorm(p)}
  res <- A+B*(1+C*tanh(g*z/2))*z*(1+z^2)^k
  res
}

#' @param log logical should the result be returned as log(). Default is FALSE
#'
#' @rdname gnk
#' @export
#'
fgnk <- function(p,A,B,C=0.8,g,k, log=FALSE, zscale=FALSE){
  stopifnot(B>0 && k>=-0.5)
  if(zscale){z <-p} else {z <- qnorm(p)}
  dQdz <- B*(1+z^2)^k*((1+C*tanh(g*z/2))*(1+(2*k+1)*z^2)/(1+z^2)+(C*g*z)/(2*cosh(g*z/2)^2))
  res <- dQdz * fnorm(p)
  if (log) return(log(res))
  res
}

#' @rdname gnk
#' @export
dqgnk <- function(p,A,B,C=0.8,g,k, log=FALSE, zscale=FALSE){
  stopifnot(B>0 && k>=-0.5)
  res <- fgnk(p, A,B, C, g, k, log=FALSE, zscale)
  if (log) return(log(1/res))
  1/res
}

#' @param n numeric number of samples to draw
#' @rdname gnk
#' @export
rgnk <- function(n, A,B,C=0.8,g,k){
  qgnk(runif(n), A,B,C,g,k)
}

#' @param q quantile value for which the corresponding cumulative probability value should be found.
#' @param tol numeric tolerance balue for approximating CDF. Default 1e-15
#' @param maxiter numeric maximum number of iteration
#' @param n_grid integer size of helper grid to be passed to `make_pgrid`. Default is 50
#' @param s_grid integer beta shape of helper grid to be passed to `make_pgrid`. Default is 5
#' @importFrom stats pnorm qnorm uniroot
#' @rdname gnk
#' @export
pgnk <- function(q, A,B,C=0.8, g, k, n_grid=100L, s_grid=5L, tol=1e-15, maxiter=1e3, zscale=FALSE){
  afun <- function(z, q, A,B, C, g, k) {q - qgnk(z,A,B,C,g,k, zscale=TRUE)}
  stopifnot(B>0 && k>=-0.5)
  mtol <- .Machine$double.eps
  p_grd <- make_pgrid(n=n_grid, s=s_grid)
  q_grd <- qgnk(p_grd, A,B,C, g,k, zscale=FALSE)
  idx_lower <- findInterval(q, q_grd, all.inside = TRUE)
  idx_upper <- idx_lower+1L
  int_lower <- stats::qnorm(p_grd[idx_lower])
  int_upper <- stats::qnorm(p_grd[idx_upper])
  zs <- mapply(function(.q, .il, .iu) {
    tmp_zs <- NULL
    tmp_zs <- stats::uniroot(afun, q=.q, A=A, B=B, C=C, g=g, k=k,
                             interval=c(.il,.iu), extendInt="downX", check.conv=TRUE, tol = tol, maxiter = maxiter)
    if(is.null(tmp_zs)) res <- NA else res <- tmp_zs$root
    res
  }, q, int_lower, int_upper)
  ps <- stats::pnorm(zs)

  ps[!is.finite(ps)] <- NA
  ps
}

#' @export
#' @rdname gnk
#' @importFrom stats uniroot
#' @examples
#' is_gnk_valid(A=5, B=5, C=0.8, g=0.5, k=-0.3) #FALSE
#' is_gnk_valid(A=5, B=5, C=0.8, g=0.5, k=-0.1) #TRUE
is_gnk_valid <- function(A, B, C=0.8, g, k, n_grid=100L, s_grid=2L){
  grd <- make_pgrid(n_grid, s_grid)
  all(fgnk(grd, A, B, C, g, k, zscale = FALSE)>0)
}


