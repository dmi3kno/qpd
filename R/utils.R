#' internal function for calculating the logit
#' @keywords  internal
logit <- function(x){log(x)-log1p(-x)}

#' internal function for calculating the inverse logit
#' @keywords  internal
invlogit <- function(x){1/(1+exp(-x))}

#' internal function for checking if the number is odd
#' @keywords  internal
is.odd <- function(x){round(x,0) %% 2 != 0}

#' internal function for splitting the integer value into a homogenous vector
#' @param x integer value to split
#' @param d integer divisor
#' @return vector of length d of integers
#' @keywords  internal
split_int <- function(x, d){
  rep(x %/% d, d)+
    c(rep(1, x%%d), rep(0, d-x%%d))
}

#' internal function for hyperbolic secant borrowed from `pracma`
#' @param x numeric vector
#' @return vector of hyperbolic secants
#' @keywords  internal
sech <- function(x){
  stopifnot(is.numeric(x) || is.complex(x))
  1 / cosh(x)  # 2 / (exp(x) + exp(-x))
}

#' internal function for hyperbolic cotangent borrowed from `pracma`
#' @param x numeric vector
#' @return vector of hyperbolic secants
#' @keywords  internal
coth <- function(x){
  stopifnot(is.numeric(x) || is.complex(x)) 
  cosh(x)/sinh(x)
}

#' internal function for hyperbolic cosecant borrowed from `pracma`
#' @param x numeric vector
#' @return vector of hyperbolic secants
#' @keywords  internal
csch <- function(x){
  stopifnot(is.numeric(x) || is.complex(x)) 
  1/sinh(x)
}

#' Diagonal matrix with offset
#'
#' @param x a matrix, vector or 1D array, or missing.
#' @param nrow,ncol optional dimensions for the result when x is not a matrix.
#' @param names (when x is a matrix) logical indicating if the resulting vector,
#' the diagonal of x, should inherit names from dimnames(x) if available.
#' @param offset an integer value to offset the diagonal by. Default is 1
#'
#' @return a square matrix
#' @importFrom utils head tail
#' @export
diago <- function(x=1, nrow=NA, ncol=NA, names=TRUE, offset = 1L){
  if(is.na(nrow) && is.na(ncol))
    nrow <- ncol <- length(x)
  m <- diag(x, nrow, ncol, names)
  if(offset==0) return(m)
  nr  <- nrow(m)
  nc <- ncol(m)
  col_idx <- seq(nc)
  ma <- matrix(rep(0, times=nr*abs(offset)), nrow=nr)
  pre <- post<- NULL
  if(offset>0) {
    pre <- ma
    s_col_idx <- utils::head(col_idx,-offset)
  } else {
    post <-ma
    s_col_idx <- utils::tail(col_idx,offset)
  }

  cbind(pre,
        m[, s_col_idx],
        post)
}

#' internal function for monotonic power
#' @references Ahmadabadi MN, Farjami Y, Moghadam MB. 2012. Approximating Distributions by Extended Generalized Lambda Distribution (XGLD). Communications in Statistics - Simulation and Computation. 41(1):1–23. https://doi.org/10.1080/03610911003681503
#' @param x numeric vector of centered depths i.e. x=u-0.5
#' @param pow numeric power to raise the depth into
#' @return monopowered centered depth
#' @keywords  internal
monopower <- function(x, pow=1){
  ifelse(x>=0, x^pow, -abs(x)^pow)
}

