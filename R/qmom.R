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
#' @rdname qmomqf
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
#' @rdname qmomqf
qsk_qf <- function(qf,...,lambda=0.25){
  stopifnot(lambda<0.5)
  num <- (qf(1-lambda,...)+qf(lambda,...)-2*qf(0.5,...))
  denom <-  (qf(1-lambda,...)-qf(lambda,...))
  num/denom
}

#' @export
#' @rdname qmomqf
qkr_qf <- function(qf, ..., lambda=0.25, zeta=lambda/2){
  stopifnot(lambda<0.5, zeta<lambda)
  num <- (qf(1-lambda+zeta,...)-qf(1-lambda-zeta,...)) +
    (qf(lambda+zeta,...)-qf(lambda-zeta,...))
  denom <-  (qf(1-lambda,...)-qf(lambda,...))
  num/denom
}

#' Unnormalized empirical robust moments.
#' \deqn{\sigma_r=Q(3/4)-Q(1/4)}
#' \deqn{s_r=\frac{Q(3/4)+Q(1/4)-2Q(1/2)}{\sigma_r}}
#' \deqn{\kappa_r=\frac{Q(7/8)-Q(5/8)+Q(3/8)-Q(1/8)}{\sigma_r}}
#' These are implemented as `qdv_eqf`, `qsk_eqf`, and `qkr_eqf`, respectively. Note that the robust measure of location is the median \eqn{\mu_r=Q(1/2)}.
#'
#' @param x numerical sample to compute the quantiles from
#' @param type parameters, passed to `quantile` function
#' @param lambda probability < 0.5 corresponding to the tail in the robust deviation. Default is 0.25
#' @param zeta probability `zeta<lambda` corresponding to the tail in robust kurtosis. Default is `lambda/2`
#' @return numeric value of robust moment
#' @export
#' @rdname qmomeqf
#'
#' @examples
#'
#' qdv_eqf(1:100) #49.5
#' qsk_eqf(1:100)
#' qkr_eqf(1:100)

qdv_eqf <- function(x, type=5,lambda=0.25){
  stopifnot(lambda<0.5)
  qs <- quantile(x, probs=c(lambda, 1-lambda), type=type, names=FALSE)
  num <- (qs[2]-qs[1])
  num
}

#' @export
#' @rdname qmomeqf
qsk_eqf <- function(x, type=5, lambda=0.25){
  stopifnot(lambda<0.5)
  qs <- quantile(x, probs=c(lambda, 0.5, 1-lambda), type=type, names=FALSE)
  num <- (qs[1]+qs[3]-2*qs[2])
  denom <- (qs[3]-qs[1])
  num/denom
}

#' @export
#' @rdname qmomeqf
qkr_eqf <- function(x, type=5, lambda=0.25, zeta=lambda/2){
  stopifnot(lambda<0.5, zeta<lambda)
  qs <- quantile(x, probs=c(lambda-zeta, lambda, lambda+zeta,
                          1-lambda-zeta, 1-lambda, 1-lambda+zeta),
                 type=type, names=FALSE)
  num <- (qs[6]-qs[4]) + (qs[3]-qs[1])
  denom <- (qs[5]-qs[2])
  num/denom
}
