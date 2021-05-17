
#' Make probability grid using beta-distribution method
#'
#' @param n integer length of the grid
#' @param s beta distribution shape parameter. For uniform grid choose 1. Default 10.
#' @param trim logical, should the 0 and 1 (tails of the grid) be trimmed. Default is `TRUE`
#'
#' @return probability grid vector of length n
#' @export
#'
#' @examples
#' make_pgrid(100)
#' make_pgrid(100, 1, FALSE) #uniform grid including 0 and 1
#' @importFrom utils head tail
#' @importFrom stats pbeta
make_pgrid <- function(n, s=10, trim=TRUE){
  if(trim) n <- n + 2L
  x<- seq(0, 1, length.out=n)
  res <- stats::pbeta(x, s, s)
  if(trim) return(utils::tail(utils::head(res, -1),-1))
  res
}
