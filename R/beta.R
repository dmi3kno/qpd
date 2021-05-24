#' Fit beta distribution to symmetrical quantile triplet of probabilities
#'
#' @param lq,mq,uq  values of lower, median and upper quantiles. All values should be between 0 and 1.
#' @param alpha share of density in each tail. Default is 0.25, which means that `lq` corresponds to 25th quantile and `uq` to 75th quantile.
#' @param ... other arguments passed to `stats::optim`
#'
#' @return data.frame with one row and two columns: `alpha` and `beta`
#' @export
#'
#' @examples
#' fit_beta(0.6, 0.65, 0.71)
fit_beta <- function(lq, mq, uq, alpha=0.25, ...){

  stopifnot(all(c(lq, mq, uq)<=1 &
                c(lq, mq, uq)>=0) )

  Q <- function(ab, lq, mq, uq, alpha){
    (stats::pbeta(lq, exp(ab[[1]]), exp(ab[[2]]))-alpha)^2+
      (stats::pbeta(mq, exp(ab[[1]]), exp(ab[[2]]))-0.5)^2+
      (stats::pbeta(uq, exp(ab[[1]]), exp(ab[[2]]))-(1-alpha))^2
  }

  iqr_norm <- stats::qnorm(1-alpha)-stats::qnorm(alpha)

  cc <- (iqr_norm^2)/4*
    ((lq*(1-mq))^0.5-
       (mq*(1-lq ))^0.5+
       (mq*(1-uq ))^0.5-
       (uq*(1-mq ))^0.5
    )^(-2)

  init_ab <- list(a=log(cc*mq+1/4),
                  b=log(cc*(1-mq)+1/4))

  res <- stats::optim(par=init_ab, fn=Q, ...,
               lq=lq, mq=mq, uq=uq, alpha=alpha)

  data.frame(alpha=exp(res$par[1]), beta=exp(res$par[2]), stringsAsFactors = FALSE, row.names = NULL)
}
