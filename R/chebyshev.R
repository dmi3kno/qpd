
#' Transform values from any interval to \[-1, 1\] and back
#'
#' @param x,z vector of values on regular and \[-1,1\] scale, respectively
#' @param interval vector of length 2 indicating the interval to which x values belong. Default is c(0,1)
#'
#' @return returns values on the other scale
#' @export
#' @rdname m11
#' @examples
#' m11_to_ab(c(-0.5, 0, 0.5))
#' ab_to_m11(c(0.25, 0.5, 0.75))
m11_to_ab <-function(z, interval=c(0,1)){
  stopifnot(length(interval)==2 && interval[2]>interval[1])
  bma <- diff(interval)
  bpa <- sum(interval)
  z*bma/2+bpa/2
}


#' @export
#' @rdname m11
ab_to_m11 <- function(x, interval=c(0,1)){
  stopifnot(length(interval)==2 && interval[2]>interval[1])
  bma <- diff(interval)
  bpa <- sum(interval)
  (2*x-bpa)/bma
}

#' internal function for making a cosine grid
#' @keywords  internal
make_cos_grid <- function(n, j=1){
  N <- n+1
  k <- c(0:n)
  t(sapply(j, function(j) cos(j*pi*(k+0.5)/N)))
}

#' Approximate a function by Chebyshev polynomials and find proxy roots
#'
#' @param fun function to evaluate
#' @param interval interval for input argument
#' @param n degree of Chebyshev polynomial
#' @param ... other arguments passed to a function
#' @param adjusted should Chebyshev polynomial be adjusted
#'
#' @return a vector of Chebyshev polynomial coefficients
#' @export
#' @rdname chebyshev
#'
#' @examples
#' get_chebyshev_coeffs(qpd::fexp, lambda=0.25, n=5)
get_chebyshev_coeffs <- function(fun, interval=c(0,1), n, ..., adjusted=FALSE){
  # Numerical Recipes 3rd ed, p234, equation 5.8.7
  N <- n+1
  js <- c(0:n)
  zk <- make_cos_grid(n, j=1)
  xk <- m11_to_ab(zk, interval = interval)
  #evaluate a function at these points
  f <- drop(fun(xk, ...)) # 1 x N
  Ijk <- 2/N * make_cos_grid(n,js)
  a <-  Ijk %*% f
  if(adjusted) a[1] <- a[1]/2
  drop(a)
}
#' @param a numeric vector of Chebyshev coefficients
#' @export
#' @rdname chebyshev
find_chebyshev_roots <- function(a, interval=c(0,1)){
  n <- length(a)-1
  vone <-rep(1, n-1)
  # correcting the absolute Chebyshev term
  a[1] <- a[1]/2
  # define Frobenius-Chebyshev companion matrix
  A <- diago(c(vone/2, 1/2), offset=-1L) +
    diago(c(1, vone/2), offset=+1L)
  CF <- A+ rbind(matrix(0, nrow=n-1, ncol=n),
                 (-0.5)*a[1:n]/a[n+1])
  roots <- eigen(CF)$values
  #cat("roots:", Re(roots[Im(roots)==0]), "\n")
  res_m11 <- Re(roots[Im(roots)==0])
  res_ab <- m11_to_ab(res_m11, interval=interval)
  sort(res_ab[res_ab>=min(interval) & res_ab <= max(interval)])
}

#' @param a numeric vector of Chebyshev coefficients
#' @param x numeric vector of values X to evaluate approximating polynomials for
#' @export
#' @rdname chebyshev
eval_chebyshev_poly <- function(a, x, interval = c(0,1)){
  # n is up-to and including
  # evaluate at value x
  stopifnot(is.numeric(n) & n>0 & n%%1==0)
  # remap vector to [-1,1] space
  z <- ab_to_m11(x, interval=interval)
  N <- length(a) # order of polynomial
  n<- N-1L
  T_k <- matrix(0, N,N)
  T_k[1,1] <- 1
  T_k[2,2] <- 1
  if(n>=2)
    for(i in seq.int(3, N, by=1))
      T_k[i,] <- 2*c(0, T_k[i-1, 1:n])-T_k[i-2,]
  T_k <- T_k[, N:1] # reverse the column order to decreasing power
  res <- outer(z[1:length(z)], n:0, "^") %*% t(T_k)
  res %*% a - a[1]/2
}

#' @param a numeric vector of Chebyshev coefficients
#' @export
#' @rdname chebyshev
fasteval_chebyshev_poly <- function(a, x, interval=c(0,1)){
  N <- length(a) # order of polynomial
  # clenshaw-horner recurrence
  z <- ab_to_m11(x, interval = interval)
  b1 <- b2 <- 0
  for (j in 1:N){
    b0 <- 2*z*b1-b2+a[N+1-j]; b3 <- b2; b2 <- b1; b1 <- b0
  }
  (0.5*(b0-b3) + 0.5*a[1]) -0.5*a[1] #adjusted
}
