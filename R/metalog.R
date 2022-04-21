# internal function for preparing the vector z for metalog.
metalog_prepare_z <- function(q, bl, bu){
  if(is.finite(bl) || is.finite(bu)){
    z <-log((q-bl)/(bu-q))
    if(is.infinite(bu)) z<-log(q-bl) # bl is defined
    if(is.infinite(bl)) z<-(-log(bu-q))# bu is defined
  } else {return(q)}
  z
}

# internal function for preparing the matrix Y for metalog.
metalog_prepare_Y <- function(p, n, nterms){
  log_odds <- log(p)-log1p(-p) #stats::qlogis(p) #log(p/(1-p))
  pmhalf <- p-0.5
  Y <- matrix(NA_real_, n, nterms)
  Y[,1] <- 1
  odd <- TRUE
  if (nterms>1) Y[,2] <- log_odds;
  if (nterms>2) Y[,3] <- pmhalf*log_odds;
  if (nterms>3) Y[,4] <- pmhalf;
  if (nterms>4)
    for (m in 5:nterms){
      if(odd){
        pmhalf <- pmhalf*pmhalf;
        Y[,m] <- pmhalf;
      } else {
        Y[,m] <- pmhalf*log_odds
      }
      odd <- isFALSE(odd)
    }
  Y
}

#' @title Fitting the metalog functions
#' @description Functions for fitting and sampling from metalog distribution
#' @details
#' `fit_metalog` is for fitting the metalog function to the set of QP values.
#' Number of metalog terms will match the number of QP pairs.
#' `approx_metalog` is for approximating metalog function to the set of data.
#' `is_metalog_valid` is a function for checking if the metalog is valid
#' @param a vector of `a`-coefficient parameters of metalog distribution
#' @param p vector of cumulative probabilities corresponding to quantile values
#' @param q vector of quantile values (data)
#' @param bl real value of lower boundary (for bounded metalog). Default -Inf
#' @param bu real value of upper boundary (for bounded metalog). Default Inf
#' @param log.p are probabilities provided on log scale. Default FALSE
#' @rdname fit_metalog
#' @export
#' @examples
#' p <- c(0.1, 0.5, 0.9)
#' q <- c(4, 9, 12)
#' fit_metalog(p,q)
fit_metalog <- function(p, q, bl=-Inf, bu=Inf, log.p=FALSE){
  n <- length(q)
  if(length(p)!=n) stop("Length of p and q must be equal")
  if(log.p) p <- exp(p)
  z <- metalog_prepare_z(q, bl, bu)
  Y <- metalog_prepare_Y(p, n, n)
  solve(Y, z)
}

#' @param nterms integer number of terms for approximating metalog. Default is 3
#' @param p_grid the grid of probability values the quantiles q correspond to. This would be specified if metalog is fitted to the empirical CDF. Default is NULL.
#' @param thin logical. Should original data be thinned. Default is FALSE.
#' @param n_grid in case data thinning is performed, integer number of quantiles to extract from data, if data vector `q` is longer than this value
#' @param s_grid in case data thinning is performed, probability grid shape parameter passed to `qpd::make_pgrid()`. Default is 10.
#' @param tol tolerance for qr.solve()
#' @rdname fit_metalog
#' @importFrom stats quantile
#' @export
approx_metalog <- function(q, nterms=3L, bl=-Inf, bu=Inf, p_grid=NULL, thin=FALSE, n_grid=1e3, s_grid=2L, tol=1e-7){
  stopifnot("Metalog should have at least 2 terms!"=nterms>1)
  n <- length(q)
  if (is.null(p_grid)){ # do it if p_grid is not specified
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
    stopifnot("Lengths of q and p_grid vectors do not match"=length(p_grid)==n)
    idx <- p_grid!=0 & p_grid!=1
    p <-p_grid[idx]
    qs <- q[idx]
  }

  if(length(nterms)!=1L) stop("Incorrectly specified number of metalog terms")
  z <- matrix(metalog_prepare_z(qs[is.finite(qs)], bl, bu), ncol=1, byrow=FALSE)
  Y <- metalog_prepare_Y(p[is.finite(qs)], n, nterms)
  M4i <- t(Y) %*% Y
  Mi <- tryCatch({
     #message("Trying to invert the matrix with solve()")
     solve(M4i)
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

#' @rdname fit_metalog
#' @importFrom stats  median
#' @export
approx_max_metalog <- function(q, bl=-Inf, bu=Inf, p_grid=NULL, thin=FALSE, n_grid=1e3, s_grid=2L, tol=1e-7){
  nterms <- 2
  a <- stats::median(q)
  metalog_valid <- TRUE
  while(metalog_valid){
    tmp_a <- approx_metalog(q, nterms =nterms, bl=bl, bu=bu, p_grid=p_grid, thin=thin, n_grid=n_grid, s_grid=s_grid, tol=tol)
    metalog_valid <- is_metalog_valid(tmp_a, bl=bl, bu=bu)
    if(metalog_valid) a <- tmp_a
    nterms <- nterms+1L
  }
  a
}

#' @title Metalog distribution functions
#' @description Functions for sampling from metalog distribution
#' @details
#' `qmetalog` is a quantile function.
#' `fmetalog` is a quantile density function q(u).
#' `dqmetalog` is the reciprocal of it is density quantile function f(Q(p)).
#' `pmetalog` is an approximation of the cumulative density function.
#' `rmetalog` is an RNG.
#' @param a vector of `a`-coefficient parameters of metalog distribution
#' @param p vector of cumulative probabilities corresponding to quantile values
#' @param bl real value of lower boundary (for bounded metalog). Default -Inf
#' @param bu real value of upper boundary (for bounded metalog). Default Inf
#' @param log.p are probabilities provided on log scale. Default FALSE
#' @param log should the log density be returned. Default FALSE
#' @rdname metalog
#' @export
#' @examples
#' a <- c(9,  1.8, -1.13)
#' p <- c(0.1, 0.5, 0.9)
#' qmetalog(p, a)
qmetalog <- function(p, a, bl=-Inf, bu=Inf, log.p=FALSE){
  n <- length(a)
  if(log.p) p <- exp(p)
  res <- a[1]
  logitp <- log(p/(1-p)) #stats::qlogis(p) #log(p/(1-p))
  pmhalf <- p-0.5
  odd <- TRUE
  #stopifnot(n>0)
  if(n>1) res <- res+a[2]*logitp
  if(n>2) res <- res+a[3]*pmhalf*logitp
  if(n>3) res <- res+a[4]*pmhalf
  if(n>4)
   for(m in 5:n){
    if(odd){
      res = res+a[m]*pmhalf^((m-1)/2)
    }else{
      res = res+a[m]*pmhalf^(m/2-1)*logitp
    }
    odd=isFALSE(odd)
   }
  res <- ifelse(p==0, bl,
                  ifelse(p==1, bu, res))
  if(is.infinite(bl) && is.infinite(bu)) return(res)
  if(is.infinite(bu)) return(ifelse(p==0, bl, bl+exp(res)))# bl is defined
  if(is.infinite(bl)) return(ifelse(p==1, bu, bu-exp(-res)))# bu is defined
  (bl+bu*exp(res))/(1+exp(res))  #logitmetalog case
}

#' @rdname metalog
#' @export
#' @examples
#' fmetalog(p, a)
fmetalog <- function(p, a, bl=-Inf, bu=Inf, log.p=FALSE, log=FALSE){
  n <- length(a)
  if(log.p) p <- exp(p)
  pt1mp <- p*(1-p)
  logitp <- log(p/(1-p)) #stats::qlogis(p) #log(p/(1-p))
  pmhalf <- p-0.5
  odd <- TRUE
  #stopifnot(n>0)
  if(n>1) res <- a[2]/pt1mp
  if(n>2) res <- res+a[3]*(pmhalf/pt1mp+logitp)
  if(n>3) res <- res+a[4]
  if(n>4)
   for(m in 5:n){
    if(odd){
      res <- res+a[m]*(m-1)/2*pmhalf^((m-3)/2)
    }else{
      res<- res+a[m]*(pmhalf^(m/2-1)/pt1mp+(m/2-1)*pmhalf^(m/2-2)*logitp)
    }
    odd=isFALSE(odd)
   }
  res <- ifelse(p==0, bl,
                ifelse(p==1, bu, res))
  if(is.infinite(bl) && is.infinite(bu)){if(log) return(log(res)) else return(res)}
  eQm <- exp(qmetalog(p,a))
  if(is.infinite(bu)) {res <- ifelse(p==0, 0, (res*eQm)); if(log) return(log(res)) else return(res)}# bl is defined
  if(is.infinite(bl)) {res <- ifelse(p==1, 0, (res/eQm)); if(log) return(log(res)) else return(res)}# bu is defined
    #both are defined, logitmetalog case
   res <- res*(bu-bl)*eQm/(1+eQm)^2
  #res <- ifelse(p==0 | p==1, 0, res*(bu-bl)*eQm/(1+eQm)^2)
  if(log) return(log(res))
  res
}

#' @rdname metalog
#' @export
#' @examples
#' dqmetalog(p, a)
dqmetalog <- function(p, a, bl=-Inf, bu=Inf, log.p=FALSE, log=FALSE){
  res <- 1/fmetalog(p,a,bl, bu, log.p, log=FALSE)
  if(log) return(ifelse(is.finite(res),log(res),res))
  res
}

#' @param q real vector of values
#' @param n_grid integer size of helper grid to be passed to `make_pgrid`. Default is 50
#' @param s_grid integer beta shape of helper grid to be passed to `make_pgrid`. Default is 5
#' @param tol tolerance value, default is 1e-6
#' @param maxiter maximum number of iterations for approximation, default is 1e6
#' @param log.p should log probability be returned
#' @rdname metalog
#' @export
#'
#' @examples
#' x <- c(5, 9, 14)
#' pmetalog(x, a)
#' @importFrom stats approx

pmetalog <- function(q, a, bl=-Inf, bu=Inf, n_grid=50L, s_grid=2L, tol=1e-15, maxiter=1e3, log.p=FALSE){

  afun <- function(u, q, a, bl, bu) {q - qmetalog(u, a, bl, bu)}
  p_grd <- sort(c(tol, qpd::make_pgrid(n=n_grid, s=s_grid), 1-tol))
  q_grd <- qmetalog(p_grd, a, bl, bu)
  idx_lower <- findInterval(q, q_grd, all.inside = TRUE)
  idx_upper <- idx_lower+1L
  int_lower <- p_grd[idx_lower]
  int_upper <- p_grd[idx_upper]
  ps <- mapply(function(.q, .il, .iu) {
    tmp_us <- NULL
    tmp_us <- stats::uniroot(afun, q=.q, a=a, bl=bl, bu=bu,
                             interval=c(.il, .iu), extendInt="no", check.conv=TRUE, tol = tol, maxiter = maxiter)
    if(is.null(tmp_us)) res <- NA else res <- tmp_us$root
    #tmp_ps
  },  q, int_lower, int_upper)
  ps <- pmin(1, pmax(0, ps))

  ps[!is.finite(ps)] <- NA
  if(log.p) return(log(ps)) else return(ps)
}

#' @param n integer value correponding to number of samples to draw
#' @rdname metalog
#' @export
#'
#' @examples
#' a <- c(2.4, 0.4, -0.08)
#' rmetalog(100, a)
#' @importFrom stats runif
rmetalog <- function(n, a, bl=-Inf, bu=Inf){
  qmetalog(stats::runif(n), a, bl,bu)
}

#' @rdname fit_metalog
#' @examples
#' a <- c(9,  1.8, -1.13)
#' is_metalog_valid(a)
#' @export
is_metalog_valid <- function(a, bl=-Inf, bu=Inf){
 grd <- make_tgrid()
 all(fmetalog(grd, a, bl, bu)>0)
}
