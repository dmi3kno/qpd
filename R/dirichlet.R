#' Generalized Dirichlet (Connor-Mosimann) distribution
#'
#' @param x simplex vector x of length n
#' @param a parameter vector a of length n-1
#' @param b parameter vector b of length n-1
#' @rdname gendir
#'
#' @export
#'
#' @examples
#' x <- c(0.15, 0.40, 0.05, 0.40)
#' a <- c(4.6, 6.7, 1.32)
#' b <- c(39.4, 7.77, 5.10)
#' dgendir(x,a,b)
dgendir <- function(x, a, b){
  if(!inherits(x, "matrix"))
    x <- matrix(x, 1, length(x))
  k <- ncol(x)-1 # length of simplex
  res <- rep.int(1, nrow(x))
  for (i in 1:(k-1))
    res <- res * x[,i]^(a[i]-1)*(1-rowSums(x[, 1:i, drop=FALSE]))^(b[i]-a[i+1]-b[i+1])/beta(a[i], b[i])
  res <- res * x[,k]^(a[k]-1)*(1-rowSums(x[,1:k, drop=FALSE])^(b[k]-1))/beta(a[k],b[k])
  res
}


#' @param n number of observations to draw
#' @rdname gendir
#' @export
#' @importFrom  stats rbeta
#' @examples
#' rgendir(1,a,b)
rgendir <-function(n, a, b){
  k <- length(a)
  s <- stats::rbeta(n, a[1],b[1])
  x <- matrix(0, n, k+1)
  x[,1] <- s
  for (j in 2:k){
    x[,j] <- stats::rbeta(n, a[j], b[j])*(1-s)
    s=s+x[,j]
  }
  x[,k+1] <- 1-rowSums(x[,1:k, drop=FALSE])
  x
}

#' Dirichlet distribution
#'
#' @param n number of observations to draw from the Dirichlet distribution
#' @param k vector of parameters
#'
#' @return returns a matrix of samples from the Dirichlet distribution
#' with number of rows n for `rdir()` or a vector or densities for `ddir()`
#' @rdname dir
#' @export
#' @importFrom stats rgamma
#' @examples
#' k <- c(3.14, 10.6, 2.24, 8.67)
#' rdir(20, k)
rdir <- function(n,k) {
  # Simulations from the Dirichlet Distribution,
  # according to the method of Wikipedia
  lk <- length(k)
  sim <- matrix(0,n,lk)
  gams <- matrix(0,n,lk)
  for (i in 1:lk) {
    gams[,i] <- matrix(stats::rgamma(n,k[i]),n,1)
  }
  gamtotal <- matrix(rowSums(gams),n,lk)
  gams/gamtotal
}


#' Title
#'
#' @param x matrix `n*m` of `n` (rows) simplices of dimension `m` (columns)
#'
#' @rdname dir
#' @export
#'
#' @examples
#' x <- rdir(5, k)
#' ddir(x, k)
ddir <- function(x, k){
  stopifnot(ncol(x)==length(k))
  bk <- prod(gamma(k))/gamma(sum(k))
  apply(x, 1, function(xx) prod(xx^(k-1))/bk)
}


