#' Function factory for creating an inverse quantile function
#'
#' @param QFUN bare name of the quantile function to be inverted. Assumes that the first argument in the quantile function is `p` for probability.
#'
#' @return function. Approximated CDF for a given quantile function.
#' The inverse function will have the following additional arguments, passed to stats::uniroot()
#' lower=0, upper=1, tol=1e-6, silent=TRUE.
#' The resulting inverse quantile function is fully vectorized with regards to
#' all of its arguments (shorter vectors are recycled).
#'
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
   dots <- list(...)
   dots_rec <- lapply(dots, rep, length.out=length(q)) # crude recycling
   dots_t <- do.call(Map, c(f = list, dots_rec))
   afun <- function(x, p, .fun, .arglst) {x - do.call(.fun, c(p, .arglst))}

   ps <- sapply(seq_along(q), function(i) {
    tmp_ps <- NULL
    tmp_ps <- try(stats::uniroot(afun, .fun=QFUN, .arglst=dots_t[[i]],
                                 lower=lower, upper = upper, x=q[i], tol = tol),
                  silent=silent)
    ifelse(is.null(tmp_ps) || inherits(tmp_ps, "try-error"), NA, tmp_ps$root)
    })

  ps[ps < 0] <- 0
  ps[ps > 1] <- 1
  ps[!is.finite(ps)] <- NA
  ps
  }
}
