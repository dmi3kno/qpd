#' Generalized Lambda Distribution (GLD) RS parameterization
#'
#' @param u numeric vector of probabilities
#' @param l1 GLD parameter \eqn{\lambda_1}
#' @param l2 GLD parameter \eqn{\lambda_2}
#' @param l3 GLD parameter \eqn{\lambda_3}
#' @param l4 GLD parameter \eqn{\lambda_4}
#'
#' @return quantiles, QDF, DQF, random samples or probabilities of GLD
#' @rdname gldrs
#' @export
#'
#' @examples
#' p_grd <- make_pgrid()
#' qgldrs(p_grd, 1, -1, -1/8, -1/32)
qgldrs <- function(u, l1, l2, l3, l4){
  res <- l1+(u^l3-(1-u)^l4)/l2
  res
}
#' @rdname gldrs
#' @export
#' @param log should the log density be returned. Default=FALSE
fgldrs <- function(u, l1, l2, l3, l4, log=FALSE){
  res <- (l3*u^(l3-1)+l4*(1-u)^(l4-1))/l2
  if(log) return(log(res))
  res
}
#' @rdname gldrs
#' @export
dqgldrs <- function(u, l1, l2, l3, l4, log=FALSE){
  res <- fgldrs(u, l1, l2, l3, l4, log=FALSE)
  if(log) return(log(1/res))
  1/res
}

#' @rdname gldrs
#' @export
#' @param n numeric number of samples to draw
rgldrs <- function(n, l1, l2, l3, l4){
  qgldrs(runif(n), l1, l2, l3, l4)
}

#' @param q quantile value for which the corresponding cumulative probability value should be found.
#' @param tol numeric tolerance balue for approximating CDF. Default 1e-15
#' @param maxiter numeric maximum number of iteration
#' @param n_grid integer size of helper grid to be passed to `make_pgrid`. Default is 50
#' @param s_grid integer beta shape of helper grid to be passed to `make_pgrid`. Default is 5
#' @rdname gldrs
#' @export
pgldrs <- function(q, l1, l2, l3, l4, n_grid=50L, s_grid=5L, tol=1e-15, maxiter=1e3){

  afun <- function(u, q, l1, l2, l3, l4) {q - qgldrs(u,l1, l2, l3, l4)}
  p_grd <- sort(c(tol, qpd::make_pgrid(n=n_grid, s=s_grid), 1-tol))
  q_grd <- qgldrs(p_grd, l1, l2, l3, l4)
  idx_lower <- findInterval(q, q_grd, all.inside = TRUE)
  idx_upper <- idx_lower+1L
  int_lower <- p_grd[idx_lower]
  int_upper <- p_grd[idx_upper]
  ps <- mapply(function(.q, .il, .iu) {
    tmp_us <- NULL
    tmp_us <- stats::uniroot(afun, q=.q, l1=l1, l2=l2, l3=l3, l4=l4,
                             interval=c(.il, .iu), extendInt="no", check.conv=TRUE, tol = tol, maxiter = maxiter)
    if(is.null(tmp_us)) res <- NA else res <- tmp_us$root
    #tmp_ps
  },  q, int_lower, int_upper)
  ps <- pmin(1, pmax(0, ps))

  ps[!is.finite(ps)] <- NA
  ps
}

#' @export
#' @rdname gldrs
is_gldrs_valid <- function(l1, l2, l3, l4, n_grid=100L, s_grid=2L){
  grd <- make_pgrid(n_grid, s_grid)
  all(fgldrs(grd, l1,l2,l3,l4)>0)
}

