#' @keywords internal
wakeby_params_valid <- function(alpha, beta, gamma, delta, xi){
  c1 <- (beta+delta>=0)
  c2 <- (beta+delta>0 || (beta==0 && delta==0 && gamma==0))
  c3 <- TRUE; if(gamma>0) c3 <-(delta>0)
  c4 <- (gamma>=0)
  c5 <- (alpha+gamma >= 0)
  all(c(c1, c2, c3, c4, c5))
}

#' Wakeby distribution
#'
#' @param u numeric vector of probabilities
#' @param alpha Wakeby scale parameter \eqn{\alpha}
#' @param beta Wakeby shape parameter \eqn{\beta}
#' @param gamma Wakeby scale parameter \eqn{\gamma}
#' @param delta Wakeby shape parameter \eqn{\delta}
#' @param xi Wakeby location parameter \eqn{\xi}
#' @details Support \eqn{\xi} to \eqn{\infty} if  \eqn{\delta \geq 0, \gamma>0} or \eqn{\xi} to \eqn{\xi+(\alpha/\beta)-(\gamma/\delta)} otherwise.
#' @return quantiles, QDF, DQF, random samples or probabilities of Wakeby disribution
#' @rdname wakeby
#' @export
#'
#' @examples
#' p_grd <- make_pgrid()
#' q_grd <- qwakeby(p_grd, 5, 3, 0.1, 0.2, 0)
#' f_grd <- fwakeby(p_grd, 5, 3, 0.1, 0.2, 0)
#' dq_grd <- dqwakeby(p_grd, 5, 3, 0.1, 0.2, 0)
#' d_grd <- dwakeby(q_grd,5, 3, 0.1, 0.2, 0)
#' plot(q_grd,dq_grd, type="l", lwd=1, lty=2)
#' lines(q_grd,d_grd, col="firebrick")
#' is_wakeby_valid(10,5,1,0.3,0)
qwakeby <- function(u, alpha, beta, gamma, delta, xi){
  stopifnot(wakeby_params_valid(alpha, beta, gamma, delta, xi))
  xi+alpha/beta*(1-(1-u)^beta)-gamma/beta*(1-(1-u)^(-delta))
}

#' @param log should the log density be returned
#' @rdname wakeby
#' @export
fwakeby <- function(u, alpha, beta, gamma, delta, xi, log=FALSE){
  stopifnot(wakeby_params_valid(alpha, beta, gamma, delta, xi))
  res <- alpha*(1-u)^(beta-1)+gamma*(1-u)^(-delta-1)
  if(log) return(log(res))
  res
}

#' @rdname wakeby
#' @export
dqwakeby <- function(u, alpha, beta, gamma, delta, xi, log=FALSE){
   res <- fwakeby(u, alpha, beta, gamma, delta, xi, log=FALSE)
   if(log) return(log(1/res))
   1/res
}

#' @param n numeric number of samples to draw
#' @rdname wakeby
#' @export
rwakeby <- function(n, alpha, beta, gamma, delta, xi){
  qwakeby(runif(n), alpha, beta, gamma, delta, xi)
}


#' @param q quantile value for which the corresponding cumulative probability value should be found.
#' @param tol numeric tolerance value for approximating CDF. Default 1e-15
#' @param maxiter numeric maximum number of iteration
#' @param n_grid integer size of helper grid to be passed to `make_pgrid`. Default is 50
#' @param s_grid integer beta shape of helper grid to be passed to `make_pgrid`. Default is 5
#' @rdname wakeby
#' @export
pwakeby <- function(q, alpha, beta, gamma, delta, xi, n_grid=50L, s_grid=5L, tol=1e-15, maxiter=1e3){

  afun <- function(u, q, alpha, beta, gamma, delta, xi) {q - qwakeby(u, alpha, beta, gamma, delta, xi)}
  p_grd <- sort(c(tol, qpd::make_pgrid(n=n_grid, s=s_grid), 1-tol))
  q_grd <- qwakeby(p_grd, alpha, beta, gamma, delta, xi)
  idx_lower <- findInterval(q, q_grd, all.inside = TRUE)
  idx_upper <- idx_lower+1L
  int_lower <- p_grd[idx_lower]
  int_upper <- p_grd[idx_upper]
  ps <- mapply(function(.q, .il, .iu) {
    tmp_us <- NULL
    tmp_us <- stats::uniroot(afun, q=.q, alpha=alpha, beta=beta, gamma=gamma, delta=delta, xi=xi,
                             interval=c(.il, .iu), extendInt="no", check.conv=TRUE, tol = tol, maxiter = maxiter)
    if(is.null(tmp_us)) res <- NA else res <- tmp_us$root
  },  q, int_lower, int_upper)
  ps <- pmin(1, pmax(0, ps))

  ps[!is.finite(ps)] <- NA
  ps
}

#' @param x numeric vector of data
#' @rdname wakeby
#' @export
dwakeby <- function(x, alpha, beta, gamma, delta, xi, n_grid=50L, s_grid=5L, tol=1e-15, maxiter=1e3, log=FALSE){
  # Johnson, Kotz, Balakrishnan (1994). Continuous univariate distributions. Vol1 p. 46. ISBN 0-471-58495-9
  stopifnot(wakeby_params_valid(alpha, beta, gamma, delta, xi))
  u <- pwakeby(x, alpha, beta, gamma, delta, xi, n_grid, s_grid, tol, maxiter)
  tt <- (1-u)^(beta+delta)
  res <- ((1-u)^(delta+1))/(alpha*tt+gamma)
  if(log) return(log(res))
  res
}

#' @rdname wakeby
#' @export
is_wakeby_valid <- function(alpha, beta, gamma, delta, xi, n_grid=50L, s_grid=5L){
  grd <- make_pgrid(n_grid, s_grid)
  all(fwakeby(grd, alpha, beta, gamma, delta, xi)>0)
}

#' @rdname wakeby
#' @export
ffwakeby <- function(u, alpha, beta, gamma, delta, xi){
  stopifnot(wakeby_params_valid(alpha, beta, gamma, delta, xi))
  -gamma*(-delta-1)*(1-u)^(-delta-2)-alpha*(beta-1)*(1-u)^(beta-2)
}
