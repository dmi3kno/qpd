#' @keywords internal
partition_interval <- function(interval, parts){
  grd <- make_pgrid(n=parts+1L, 1, trim = FALSE) # linear partitioning
  bma <- interval[2]-interval[1]
  all_points <- interval[1]+bma*grd
  list(from=head(all_points, -1),
       to=tail(all_points, -1))
}


#' Transform values from any interval to \[-1, 1\] and back
#'
#' @param x,z vector of values on regular and \[-1,1\] scale, respectively
#' @param from,to the interval(s) to which x values belong. Defaults are from=0,to=1
#'
#' @return returns values on the other scale
#' @export
#' @rdname m11
#' @examples
#' m11_to_ab(c(-0.5, 0, 0.5))
#' ab_to_m11(c(0.25, 0.5, 0.75))
m11_to_ab <-function(z, from=0, to=1){
  stopifnot(length(from)==1 && length(to)==1 && to>from)
  bma <- to-from
  bpa <- to+from
  z*bma/2+bpa/2
}


#' @export
#' @rdname m11
ab_to_m11 <- function(x, from=0, to=1){
  stopifnot(length(from)==1 && length(to)==1 && to>from)
  bma <- to-from
  bpa <- to+from
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
#' @param n degree of Chebyshev polynomial
#' @param from,to the interval for the input argument
#' @param ... other arguments passed to a function
#'
#' @return a vector of Chebyshev polynomial coefficients
#' @export
#' @rdname chebyshev
#'
#' @examples
#' get_chebyshev_coeffs(qpd::fexp, lambda=0.25, n=5)
get_chebyshev_coeffs <- function(fun, n, from=0, to=1, ...){
  stopifnot(length(from)==1 && length(to)==1 && to>from)
  # Numerical Recipes 3rd ed, p234, equation 5.8.7
  N <- n+1
  js <- c(0:n)
  zk <- make_cos_grid(n, j=1)
  xk <- m11_to_ab(zk, from=from, to=to)
  #evaluate a function at these points
  f <- drop(fun(xk, ...)) # 1 x N
  Ijk <- 2/N * make_cos_grid(n,js)
  a <-  Ijk %*% f
  a[1] <- a[1]/2 # adjust the a_0 coefficient
  drop(a)
}
#' @param a numeric vector of Chebyshev coefficients
#' @export
#' @rdname chebyshev
find_chebyshev_roots <- function(a, from=0, to=1){
  stopifnot(length(from)==1 && length(to)==1 && to>from)
  N <-length(a)
  n <- N-1
  vone <-rep(1, n-1)
  # Zero-free intervals
  # Boyd JP. 2006. Computing real roots of a polynomial in Chebyshev series form through subdivision. Applied Numerical Mathematics. 56(8):1077–1091. https://doi.org/10/d5wxwt
  # message("This is a zero-free interval!")
  B0 <- sum(abs(a[-1]))
  B1 <- sum(abs((1:n)*a[-1]))
  h <- pi/n
  tj <- pi*(0:n)/h
  NR1 <- (B0 < a[1])
  NR2 <- (B1*h*0.5 < min(abs(eval_chebyshev_poly(a, cos(tj), from = -1, to=1))))
  if(NR1 && NR2) return(numeric(0))

  # define Chebyshev-Frobenius companion matrix
  A <- diago(c(vone/2, 1/2), offset=-1L) +
       diago(c(1, vone/2), offset=+1L)
  if(a[N]==0) a[N] <- .Machine$double.eps
  CF <- A + rbind(matrix(0, nrow=n-1, ncol=n),
                 (-0.5)*a[1:n]/a[N])
  roots <- eigen(CF)$values
  res_m11 <- Re(roots[Im(roots)==0]) # real roots only
  res_ab <- m11_to_ab(res_m11, from=from, to=to)
  sort(res_ab[res_ab>=from & res_ab <= to])
}

#' @param interval interval to check the roots on
#' @export
#' @rdname qdf_check
check_roots <- function(fun, ..., n=13, interval=c(0,1), parts=10){
  stopifnot(length(interval)==2 && interval[2]>interval[1])
  intrvls <- partition_interval(interval, parts=parts)
  res <- lapply(seq_along(intrvls$from), function(i){
    from <- intrvls$from[i]
    to <- intrvls$to[i]
    cC <- get_chebyshev_coeffs(fun, from=from, to=to, n=n, ...)
    roots <- find_chebyshev_roots(cC, from=from, to=to)
    roots
  })
  Reduce(`c`, res)
}

# @param a numeric vector of Chebyshev coefficients
# @param x numeric vector of values X to evaluate approximating polynomials for
# @export
# @rdname chebyshev
#eval_chebyshev_poly <- function(a, x, from=0, to=1){
#  stopifnot(length(from)==1 && length(to)==1 && to>from)
#  len_x <-length(x)
#  # remap vector to [-1,1] space
#  z <- ab_to_m11(x, from=from, to=to)
#  N <- length(a) # order of polynomial
#  n<- N-1L
#
#  za <- acos(z)
#  j <- 0:n
#  Tj <- matrix(nrow=len_x, ncol=N)
#  for(i in seq_along(j))
#    Tj[, i] <- cos(j[i]*za)
#  drop(a %*% t(Tj))
#}

#
#sloweval_chebyshev_poly <- function(a,x, from=0, to=1){
#  stopifnot(length(from)==1 && length(to)==1 && to>from)
#  # remap vector to [-1,1] space
#  z <- ab_to_m11(x, from=from, to=to)
#  N <- length(a) # order of polynomial
#  n<- N-1L
#
#  T_k <- matrix(0, N, N)
#  T_k[1,1] <- 1
#  T_k[2,2] <- 1
#  if(n>=2)
#    for(i in seq.int(3, N, by=1))
#      T_k[i,] <- 2*c(0, T_k[i-1, 1:n])-T_k[i-2,]
#  T_k <- T_k[, N:1] # reverse the column order to decreasing power
#  res <- outer(z[1:length(z)], n:0, "^") %*% t(T_k)
#  drop(res %*% a)
#}

#' @param x numeric vector of values X to evaluate approximating polynomials for
#' @export
#' @rdname chebyshev
eval_chebyshev_poly <- function(a, x, from=0, to=1){
  stopifnot(length(from)==1 && length(to)==1 && to>from)
  N <- length(a) # order of polynomial
  # clenshaw-horner recurrence
  z <- ab_to_m11(x, from=from, to=to)
  b1 <- b2 <- 0
  for (j in 1:N){
    b0 <- 2*z*b1-b2+a[N+1-j]; b3 <- b2; b2 <- b1; b1 <- b0
  }
  (0.5*(b0-b3) + 0.5*a[1]) #-0.5*a[1] # already adjusted
}

#' Validate quantile function by checking the corresponding quantile density function
#'
#' @param fun quantile density function (QDF) to check
#' @param n integer polynomial degree to fit to the function. Default is 13
#' @param parts integer number of parts to split the function range. Default is 10
#' @param ... other parameters passed to QDF function
#' @return logical value indicating if the QDF is valid (i.e. returns only non-negative values)
#' @export
#' @rdname qdf_check
#' @importFrom stats na.omit filter
is_qdf_valid <- function(fun, ..., n=13L, parts=10){
  rs <- check_roots(fun, n=n, parts=parts, interval=c(0,1), ...)
  if(length(rs)==0) return(TRUE)
  if(length(rs)==1) rs <- c(0, rs, 1)
  mid <- stats::na.omit(stats::filter(rs, rep(1, 2), sides=1)/2)
  vals <- fun(mid, ...)
  all(vals>0)
}
