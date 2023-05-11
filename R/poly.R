#' Evaluate a polynomial in reverse order
#'
#' Evaluate the polynomial using the Horner scheme
#' @param x numeric vector of values to evaluate a polynomial at
#' @param b numeric vector of coefficients (highest degree first, degree 0 last)
#' @param deriv logical. Should the value of the derivative be returned instead?. Default is FALSE.
#'
#' @return quantiles, QDF, DQF, random samples or probabilities of GLD (CSW parameterization)
#' @rdname poly
#' @export
#'
#' @examples
#' eval_poly(1:3, c(0.1, 0.2, 0.6, 1.23, 4.43))
eval_poly <- function(x, b, deriv=FALSE){
  sapply(x, function(x){
    n <- length(b)
    res <- b[n]; der <- 0
    for(i in (n-1):1){
      der <- res + x * der
      res <- b[i]+x*res
    }
    if(deriv) return(der)
    return(res)
  })
}
