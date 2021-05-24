# function for preparing the vector z for metalog.
metalog_prepare_z <- function(q, bl, bu){
  if(!is.null(bl) || !is.null(bu)){
    z <-log((q-bl)/(bu-q))
    if(is.null(bu)) z<-log(q-bl) # bl is defined
    if(is.null(bl)) z<-(-log(bu-q))# bu is defined
  } else {return(q)}
  z
}

# function for preparing the matrix Y for metalog.
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

#' @title Metalog distribution functions
#' @description Functions for fitting and sampling from metalog distribution
#' @details
#' `fit_metalog` is for fitting the metalog function to the set of QP values.
#' Number of metalog terms will match the number of QP pairs.
#' `approx_metalog` is for approximating metalog function to the set of data.
#' `qmetalog` is a quantile function.
#' `fmetalog` is a quantile density function q(u). The reciprocal of it is density quantile function f(Q(p)).
#' `pmetalog` is an approximation of the cumulative density function.
#' `rmetalog` is an RNG.
#' `is_metalog_valid` is a function for checking if the metalog is valid
#'
#' @param p vector of cumulative probabilities corresponding to quantile values
#' @param q vector of quantile values (data)
#' @param bl real value of lower boundary (for bounded metalog). Default NULL
#' @param bu real value of upper boundary (for bounded metalog). Default NULL
#' @rdname metalog
#' @export
#' @examples
#' p <- c(0.1, 0.5, 0.9)
#' q <- c(4, 9, 12)
#' fit_metalog(p,q)
fit_metalog <- function(p, q, bl=NULL, bu=NULL){
  n <- length(q)
  if(length(p)!=n) stop("Length of p and q must be equal")
  z <- metalog_prepare_z(q, bl, bu)
  Y <- metalog_prepare_Y(p, n, n)
  solve(Y, z)
}

#' @param nterms integer number of terms for approximating metalog. Default is 3
#' @param thin logical. Should original data be thinned. Default is FALSE.
#' @param thin_to integer number of quantiles to extract from data, if data vector `q` is longer than this value
#' @param s in case data thinning is performed, probability grid shape parameter passed to `qpd::make_pgrid()`. Default is 10.
#' @rdname metalog
#' @importFrom stats quantile
#' @export
approx_metalog <- function(q, nterms=3L, bl=NULL, bu=NULL, thin=FALSE, thin_to=1000L, s=2L){
  n <- length(q)
  if(thin){
    n <- thin_to
    p <- make_pgrid(thin_to, s=s, trim=TRUE)
    qs <- unname(stats::quantile(q, probs=p))
  } else {
    qs <- sort(q)
    # ecdf assignment trick to avoid 0 and 1
    p <- (seq_along(qs)-0.5)/n
  }
  if(length(nterms)!=1L) stop("Incorrectly specified number of metalog terms")
  z <- matrix(metalog_prepare_z(qs[is.finite(qs)], bl, bu), ncol=1, byrow=FALSE)
  Y <- metalog_prepare_Y(p[is.finite(qs)], n, nterms)
  ((solve(t(Y) %*% Y) %*% t(Y)) %*% z)[,1]
}

#' @param a vector of `a`-coefficient parameters of metalog distribution
#'
#' @rdname metalog
#' @export
#' @examples
#' a <- c(9,  1.8, -1.13)
#' qmetalog(p, a)
#' @importFrom stats qlogis
qmetalog <- function(p, a, bl=NULL, bu=NULL){
  n <- length(a)
  res <- a[1]
  logitp <- log(p/(1-p)) #stats::qlogis(p) #log(p/(1-p))
  pmhalf <- p-0.5
  odd <- TRUE
  #stopifnot(n>0)
  if(n>1) res = res+a[2]*logitp
  if(n>2) res = res+a[3]*pmhalf*logitp
  if(n>3) res = res+a[4]*pmhalf
  if(n>4)
   for(m in 5:n){
    if(odd){
      res = res+a[m]*pmhalf^((m-1)/2)
    }else{
      res = res+a[m]*pmhalf^(m/2-1)*logitp
    }
    odd=isFALSE(odd)
   }
  if(is.null(bl) && is.null(bu)) return(res)
  if(is.null(bu)) return(ifelse(p==0, bl, bl+exp(res)))# bl is defined
  if(is.null(bl)) return(ifelse(p==1, bu, bu-exp(-res)))# bu is defined
  ifelse(p==0, bl,
          ifelse(p==1, bu,
                       (bl+bu*exp(res))/(1+exp(res)))
                         ) #logitmetalog case
}

#' @rdname metalog
#' @export
#' @examples
#' fmetalog(p, a)
#' @importFrom stats qlogis
fmetalog <- function(p, a, bl=NULL, bu=NULL){
  n <- length(a)
  pt1mp <- p*(1-p)
  logitp <- log(p/(1-p)) #stats::qlogis(p) #log(p/(1-p))
  pmhalf <- p-0.5
  odd <- TRUE
  #stopifnot(n>0)
  if(n>1) res = a[2]/pt1mp
  if(n>2) res = res+a[3]*(pmhalf/pt1mp+logitp)
  if(n>3) res = res+a[4]
  if(n>4)
   for(m in 5:n){
    if(odd){
      res = res+a[m]*(m-1)/2*pmhalf^((m-3)/2)
    }else{
      res = res+a[m]*(pmhalf^(m/2-1)/pt1mp+(m/2-1)*pmhalf^(m/2-2)*logitp)
    }
    odd=isFALSE(odd)
   }
  if(is.null(bl) && is.null(bu)) return(res)
  eQm <- exp(qmetalog(p,a))
  if(is.null(bu)) return(ifelse(p==0, 0, (res*eQm)))# bl is defined
  if(is.null(bl)) return(ifelse(p==1, 0, (res/eQm)))# bu is defined
    #both are defined, logitmetalog case
  ifelse(p==0 | p==1, 0, res*(bu-bl)*eQm/(1+eQm)^2)
}

#' @param x real vector of values
#' @param tol tolerance value, default is 1e-6
#' @param maxiter maximum number of iterations for approximation, default is 1e6
#' @rdname metalog
#' @export
#'
#' @examples
#' x <- c(5, 9, 14)
#' pmetalog(x, a)
#' @importFrom stats approx
pmetalog <- function(x, a, bl=NULL, bu=NULL, tol=1e-6, maxiter=1e6){
  afun <- function(x, u, a, bl, bu, tol, maxiter){
    i <- 1
    while(abs(x-qmetalog(u,a,bl,bu))>tol && i <= maxiter){
      u=u+(x-qmetalog(u,a,bl,bu))/fmetalog(u,a,bl,bu)
      i<-i+1
    }
    u
  }
  p_grd <- make_pgrid(500, 25, trim=TRUE)
  q_grd <- qmetalog(p_grd, a, bl, bu)
  is_ok <- !is.na(q_grd)
  p_guess <- stats::approx(q_grd[is_ok], p_grd[is_ok], x, ties = "ordered")[["y"]]
  mapply(afun, x, p_guess, MoreArgs = list(a=a, tol=tol, bl=bl, bu=bu, maxiter=maxiter))
}

#' @param n integer value correponding to number of samples to draw
#' @rdname metalog
#' @export
#'
#' @examples
#' a <- c(2.4, 0.4, -0.08)
#' rmetalog(100, a)
#' @importFrom stats runif
rmetalog <- function(n, a, bl=NULL, bu=NULL){
  qmetalog(stats::runif(n), a, bl,bu)
}

#' @rdname metalog
#' @examples
#' is_metalog_valid(a)
#' @export
is_metalog_valid <- function(a, bl=NULL, bu=NULL){
 grd <- make_pgrid(100, s=10, trim = TRUE)
 all(fmetalog(grd, a, bl, bu)>0)
}
