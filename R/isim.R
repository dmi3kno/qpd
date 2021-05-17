#' Convert index simplex to cumulative probability
#'
#' @param s vector or matrix of simplices
#' @param idx integer vector of indices
#'
#' @return matrix of cumulative probabilities
#' @export
#'
#' @examples
#' s <- c(0.1, 0.4, 0.15, 0.35)
#' ss <- matrix(c(0.1, 0.4, 0.15, 0.35,
#'                0.13, 0.26, 0.05, 0.56),
#'                ncol=4, byrow=TRUE)
#' idx <- c(1L, 2L, 4L, 3L)
#' isim_to_cprob(s, idx)
#' isim_to_cprob(ss, idx)
isim_to_cprob <- function(s, idx){
  UseMethod("isim_to_cprob", s)
}

#' @export
isim_to_cprob.numeric <- function(s, idx){
  lidx <- length(idx)
  ss <- matrix(s, ncol=lidx, byrow = TRUE)
  isim_to_cprob(ss, idx)
}

#' @export
isim_to_cprob.matrix <- function(s, idx){
  lidx <- length(idx)
  t(apply(s[,idx, drop=FALSE],1,cumsum)[-lidx,,drop=FALSE])
}
