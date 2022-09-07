# internal function for preparing the vector z for sqn.
sqn_prepare_z <- function(q){
  q
}

# internal function for preparing the matrix Y for sqn.
sqn_prepare_Y <- function(p){
  Y <- matrix(NA_real_, 4, 4)
  Y[,1] <- 1
  Y[,2] <- qnorm(p)
  Y[,3] <- p*qnorm(p)
  Y[,4] <- p
  Y
}

#' @title Fitting the sqn functions
#' @description Functions for fitting and sampling from sqn distribution
#' @details
#' `fit_sqn` is for fitting the sqn function to the set of QP values.
#' Number of sqn terms will match the number of QP pairs.
#' `approx_sqn` is for approximating sqn function to the set of data.
#' `is_sqn_valid` is a function for checking if the sqn is valid
#' @param a vector of `a`-coefficient parameters of sqn distribution
#' @param p vector of cumulative probabilities corresponding to quantile values
#' @param q vector of quantile values (data)
#' @param tol tolerance for solve(). Default is .Machine$double.eps^2
#' @rdname fit_sqn
#' @export
#' @examples
#' p <- c(0.1, 0.5, 0.6, 0.9)
#' q <- c(4, 9, 10, 12)
#' fit_sqn(p,q)
fit_sqn <- function(p, q, tol=.Machine$double.eps^2){
  z <- sqn_prepare_z(q)
  Y <- sqn_prepare_Y(p)
  solve(Y, z, tol=tol)
}

#' @param p the vector of probability values the quantiles q correspond to. This would be specified if sqn is fitted to the empirical CDF. Default is NULL.
#' @param thin logical. Should original data be thinned. Default is FALSE.
#' @param n_grid in case data thinning is performed, integer number of quantiles to extract from data, if data vector `q` is longer than this value
#' @param s_grid in case data thinning is performed, probability grid shape parameter passed to `qpd::make_pgrid()`. Default is 10.
#' @param tol tolerance for solve() and qr.solve(), default is .Machine$double.eps^2
#' @rdname fit_sqn
#' @importFrom stats quantile
#' @export
approx_sqn <- function(q, p=NULL, thin=FALSE, n_grid=1e3, s_grid=2L, tol=.Machine$double.eps^2){
  n <- length(q)
  if (is.null(p)){ # do it if p is not specified
    if(thin && (n > n_grid)){
      n <- n_grid
      p <- make_pgrid(n_grid, s=s_grid, trim=TRUE)
      qs <- unname(stats::quantile(q, probs=p))
    } else {
      qs <- sort(q)
      # ecdf assignment trick to avoid 0 and 1
      p <- (seq_along(qs)-0.5)/n
    }
  } else {
    stopifnot("Lengths of q and q vectors do not match"=length(p)==n)
    idx <- p!=0 & p!=1
    p <-p[idx]
    qs <- q[idx]
  }

  z <- matrix(sqn_prepare_z(qs[is.finite(qs)]), ncol=1, byrow=FALSE)
  Y <- sqn_prepare_Y(p[is.finite(qs)])
  M4i <- t(Y) %*% Y
  Mi <- tryCatch({
     #message("Trying to invert the matrix with solve()")
     solve(M4i, tol=tol)
    },
    error = function(e) {
      #message("Solve failed. Trying Cholesky decomposition.")
      res <- tryCatch(
        chol2inv(chol(M4i)),
        error=function(ee){
        #  message("Cholesky also failed. Trying QR decomposition.")
          return(qr.solve(M4i, tol=tol))
        }
      )
      return(res)
    }#,
    #finally = {
    #  message("Solved!")
    #}
    )
  ((Mi %*% t(Y)) %*% z)[,1]
}

#' @title sqn distribution functions
#' @description Functions for sampling from sqn distribution
#' @details
#' `qsqn` is a quantile function.
#' `fsqn` is a quantile density function q(u).
#' `dqsqn` is the reciprocal of it is density quantile function f(Q(p)).
#' `psqn` is an approximation of the cumulative density function.
#' `rsqn` is an RNG.
#' @param a vector of `a`-coefficient parameters of sqn distribution
#' @param p vector of cumulative probabilities corresponding to quantile values
#' @param log should the log density be returned. Default FALSE
#' @rdname sqn
#' @export
#' @examples
#' a <- c(9,  1.8, -1.13, 9)
#' p <- c(0.1, 0.5, 0.9)
#' qsqn(p, a)
qsqn <- function(p, a){
  stopifnot(length(a)==4)
  a[1]+a[2]*qnorm(p)+a[3]*p*qnorm(p)+a[4]*p
}

#' @rdname sqn
#' @export
#' @examples
#' fsqn(p, a)
fsqn <- function(p, a, log=FALSE){
  a[2]*dqnorm(p)+a[3]*qnorm(p)+a[3]*p*dqnorm(p)+a[4]
}

#' @rdname sqn
#' @export
#' @examples
#' dqsqn(p, a)
dqsqn <- function(p, a, log=FALSE){
  res <- fsqn(p, a, log=FALSE)
  if(log) return(ifelse(is.finite(res),-log(res),res))
  1/res
}

#' @param q real vector of values
#' @param n_grid integer size of helper grid to be passed to `make_pgrid`. Default is 50
#' @param s_grid integer beta shape of helper grid to be passed to `make_pgrid`. Default is 2
#' @param tol tolerance value, default is 1e-6
#' @param maxiter maximum number of iterations for approximation, default is 1e6
#' @rdname sqn
#' @export
#'
#' @examples
#' x <- c(5, 9, 14)
#' psqn(x, a)
#' @importFrom stats approx
psqn <- function(q, a, n_grid=50L, s_grid=2L, tol=1e-15, maxiter=1e3){

  afun <- function(u, q, a) {q - qsqn(u, a)}
  p_grd <- sort(c(tol, qpd::make_pgrid(n=n_grid, s=s_grid), 1-tol))
  q_grd <- qsqn(p_grd, a)
  idx_lower <- findInterval(q, q_grd, all.inside = TRUE)
  idx_upper <- idx_lower+1L
  int_lower <- p_grd[idx_lower]
  int_upper <- p_grd[idx_upper]
  ps <- mapply(function(.q, .il, .iu) {
    tmp_us <- NULL
    tmp_us <- stats::uniroot(afun, q=.q, a=a,
                             interval=c(.il, .iu), extendInt="no",
                             check.conv=TRUE, tol = tol, maxiter = maxiter)
    if(is.null(tmp_us)) res <- NA else res <- tmp_us$root
    #tmp_ps
  },  q, int_lower, int_upper)
  ps <- pmin(1, pmax(0, ps))

  ps[!is.finite(ps)] <- NA
  ps
}

#' @param n integer value correponding to number of samples to draw
#' @rdname sqn
#' @export
#'
#' @examples
#' a <- c(2.4, 0.4, -0.08, 9)
#' rsqn(100, a)
#' @importFrom stats runif
rsqn <- function(n, a){
  qsqn(stats::runif(n), a)
}

#' @rdname fit_sqn
#' @examples
#' a <- c(9,  1.8, -1.13, 9)
#' is_sqn_valid(a)
#' @export
is_sqn_valid <- function(a, n_grid=50L, s_grid=2L){
 grd <- make_pgrid(n_grid, s_grid)
 all(fsqn(grd, a)>0)
}
