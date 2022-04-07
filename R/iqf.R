#' Function factory for creating an inverse quantile function
#'
#' @param QFUN bare name of the quantile function to be inverted. Assumes that the first argument in the quantile function is `p` for probability.
#'
#' @return function. Approximated CDF for a given quantile function
#' @rdname iqf
#' @export
#'
#' @examples
#' qmyfun <- function(p, b, k){b*(log(p/(1-p))+k*p)}
#' pmyfun <- iqf(qmyfun)
#' x <- qmyfun(1:9/10, 2, 3)
#' pmyfun(x, 2, 3)
iqf <- function(QFUN){
  function(q, ..., lower=0, upper=1, tol=1e-6, silent=TRUE){
   afun <- function(x, p) {x - QFUN(p, ...)}
   ps <- sapply(q, function(.q) {
    tmp_ps <- NULL
    tmp_ps <- try(stats::uniroot(afun, lower=lower, upper = upper, x=.q, tol = tol), silent=silent)
    ifelse(is.null(tmp_ps) || inherits(tmp_ps, "try-error"), NA, tmp_ps$root)
    })

  ps[ps < 0] <- 0
  ps[ps > 1] <- 1
  ps[!is.finite(ps)] <- NA
  ps
  }
}
