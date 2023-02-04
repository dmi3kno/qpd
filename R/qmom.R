#' Robust moments
#'
#' Unnormalized robust moments.
#' \deqn{\sigma_r=Q(3/4)-Q(1/4)}
#' \deqn{s_r=\frac{Q(3/4)+Q(1/4)-2Q(1/2)}{\sigma_r}}
#' \deqn{\kappa_r=\frac{Q(7/8)-Q(5/8)+Q(3/8)-Q(1/8)}{\sigma_r}}
#' These are implemented as `qdv_qf`, `qsk_qf`, and `qkr_qf`, respectively. Note that the robust measure of location is the median \eqn{\mu_r=Q(1/2)}.
#'
#' @param qf quantile function
#' @param ... parameters, passed to the of quantile function
#' @param lambda probability < 0.5 corresponding to the tail in the robust deviation. Default is 0.25
#' @param zeta probability `zeta<lambda` corresponding to the tail in robust kurtosis. Default is `lambda/2`
#' @return numeric value of robust moment
#' @export
#' @rdname qmom
#'
#' @examples
#' qdv_qf(qnorm, 0, 1)
#' qsk_qf(qnorm, 0, 1)
#' qkr_qf(qnorm, 0, 1)
qdv_qf <- function(qf,...,lambda=0.25){
  stopifnot(lambda<0.5)
  num <- (qf(1-lambda,...)-qf(lambda,...))
  num
}

#' @export
#' @rdname qmom
qsk_qf <- function(qf,...,lambda=0.25){
  stopifnot(lambda<0.5)
  num <- (qf(1-lambda,...)+qf(lambda,...)-2*qf(0.5,...))
  denom <-  (qf(1-lambda,...)-qf(lambda,...))
  num/denom
}

#' @export
#' @rdname qmom
qkr_qf <- function(qf, ..., lambda=0.25, zeta=lambda/2){
  stopifnot(lambda<0.5, zeta<lambda)
  num <- (qf(1-lambda+zeta,...)-qf(1-lambda-zeta,...)) +
    (qf(lambda+zeta,...)-qf(lambda-zeta,...))
  denom <-  (qf(1-lambda,...)-qf(lambda,...))
  num/denom
}
