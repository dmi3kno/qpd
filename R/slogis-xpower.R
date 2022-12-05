# internal function for preparing the vector z for metalog.
slogis_xpower_prepare_z <- function(q){
  q
}

# internal function for preparing the matrix Y for metalog.
slogis_xpower_prepare_Y <- function(p, n, nterms){
  Y <- matrix(NA_real_, n, nterms)
  Y[,1] <- 1
  odd <- TRUE
  if (nterms>1) Y[,2] <- log(p);
  if (nterms>2) Y[,3] <- (-log1p(-p));
  if(nterms>3)
    for(i in 4:nterms)
      if(i%%2==0) Y[,i] <- (p)^((i-2)/2)
            else  Y[,i]  <- (-(1-p)^((i-3)/2))
  Y
}

#' @title Fitting the SLxP functions
#' @description Functions for fitting and sampling from SLxP distribution
#' @details
#' `fit_slogis_xpower` is for fitting the SLxP function to the set of QP values.
#' Number of SLxP terms will match the number of QP pairs.
#' `approx_slogis_xpower` is for approximating SLxP function to the set of data.
#' `is_slogis_xpower_valid` is a function for checking if the SLxP is valid
#' @param a vector of `a`-coefficient parameters of SlxP distribution
#' @param p vector of cumulative probabilities corresponding to quantile values
#' @param q vector of quantile values (data)
#' @param log.p are probabilities provided on log scale. Default FALSE
#' @param tol tolerance for solve(). Default is .Machine$double.eps^2
#' @rdname slogis_xpower
#' @export
#' @examples
#' p <- c(0.1, 0.5, 0.9)
#' q <- c(4, 9, 12)
#' a <-fit_slogis_xpower(p,q)
fit_slogis_xpower <- function(p, q, log.p=FALSE, tol=.Machine$double.eps^2){
  n <- length(q)
  if(length(p)!=n) stop("Length of p and q must be equal")
  if(log.p) p <- exp(p)
  z <- slogis_xpower_prepare_z(q)
  Y <- slogis_xpower_prepare_Y(p, n, n)
  solve(Y, z, tol=tol)
}

#' @param nterms integer number of terms for approximating metalog. Default is 3
#' @param p the grid of probability values the quantiles q correspond to. This would be specified if metalog is fitted to the empirical CDF. Default is NULL.
#' @param thin logical. Should original data be thinned. Default is FALSE.
#' @param n_grid in case data thinning is performed, integer number of quantiles to extract from data, if data vector `q` is longer than this value
#' @param s_grid in case data thinning is performed, probability grid shape parameter passed to `qpd::make_pgrid()`. Default is 10.
#' @param tol tolerance for solve() and qr.solve(), default is .Machine$double.eps^2
#' @rdname slogis_xpower
#' @importFrom stats quantile
#' @export
approx_slogis_xpower <- function(q, nterms=3L, p=NULL, thin=FALSE, n_grid=1e3, s_grid=2L, tol=.Machine$double.eps^2){
  stopifnot("Metalog should have at least 2 terms!"=nterms>1)
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
    stopifnot("Lengths of q and p vectors do not match"=length(p)==n)
    idx <- p!=0 & p!=1
    p <-p[idx]
    qs <- q[idx]
  }

  if(length(nterms)!=1L) stop("Incorrectly specified number of metalog terms")
  z <- matrix(slogis_xpower_prepare_z(qs[is.finite(qs)]), ncol=1, byrow=FALSE)
  Y <- slogis_xpower_prepare_Y(p[is.finite(qs)], n, nterms)
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



#' Quantile function for Skew-Logistic Expanded Power quantile-parameterized distribution
#'
#' Calculate QF, QDF, DQF `qslogis_xpower`, `fslogis_xpower`, `dqslogis_xpower` for the SLxP-QPD
#' @param u numeric probability
#' @param a numeric vector of SLxP-QPD parameters.
#' @param log logical, should the result be returned as a log. Default is FALSE
#'
#' @return returns QF (`qslogis_xpower`), QDF (`fslogis_xpower`), DQF (`dqslogis_xpower`) or random variate of skew-logistic distribution
#' @export
#' @rdname slogis_xpower
#' @examples
#' p_grd <- make_pgrid()
#'
#' fslogis_xpower(p_grd, a)
#' dqslogis_xpower(p_grd, a)
qslogis_xpower <- function(u, a){
  n <- length(a)
  stopifnot(n>=3)
  #stopifnot(all(tail(a,-1)>0))
  res <- a[1] + a[2]*log(u)-a[3]*log1p(-u)
  if(n>3)
    for(i in 4:n)
      if(i%%2==0) res <- res + a[i]*(u)^((i-2)/2)
            else  res <- res - a[i]*(1-u)^((i-3)/2)
  res
}

#' @export
#' @rdname slogis_xpower
fslogis_xpower <- function(u, a, log=FALSE){
  n <- length(a)
  stopifnot(n>=3)
  #stopifnot(all(tail(a,-1)>0))
  res <- a[2]/u + a[3]/(1-u)
  if(n>3)
    for(i in 4:n)
      if(i%%2==0) res <- res + ((i-2)/2)*a[i]*(u)^((i-2)/2-1)
            else  res <- res + (i-3/2)*a[i]*(1-u)^((i-3)/2-1)
  if(log) return(log(res))
  res
}

#' @export
#' @rdname slogis_xpower
dqslogis_xpower <- function(u, a, log=FALSE){
  res <- fslogis_xpower(u, a)
  if(log) return(ifelse(is.finite(res),-log(res),res))
  1/res
}

#' @param n number of samples to draw
#' @export
#' @rdname slogis_xpower
rslogis_xpower <- function(n, a){
  qslogis_xpower(runif(n), a)
}
