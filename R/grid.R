#' internal function for making a scaled equispaced, non-bounded grid
#' @param f from value
#' @param t to value
#' @param n number of equispaced values
#' @keywords  internal
seq_grid <- function(f, t, n){
  f+(t-f)*(seq_len(n)-0.5)/n
}

#' Compute ECDF for a sample
#'
#' Function for creating an empirical CDF from a numerical sample
#' @param x numeric sample
#' @return data frame with 2 columns p - vector of empirical probabilities, q - vector of sorted original sample values
#'
#' @rdname ecdf
#' @export
#'
#' @examples
#' make_ecdf_df(runif(10))
#'
make_ecdf_df <- function(x){
  # ecdf assignment trick to avoid 0 and 1
  data.frame(p=(rank(x, na.last = "keep")-0.5)/length(x),
             q=x, stringsAsFactors = FALSE)
}

#' Make probability grid
#' @description Functions for creating a probability grid using different method
#' `make_pgrid` uses beta-distribution method
#' `make_tgrid` uses tiered linear method
#'
#' @param n integer length of the grid. Default is 50
#' @param s beta distribution shape parameter, passed to both `shape1` and `shape2` of `pbeta()`. For uniform grid choose 1. Default is 2.
#' @param trim logical, should the 0 and 1 (tails of the grid) be trimmed. Default is `TRUE`
#'
#' @return probability grid vector of length n
#' @rdname make_grid
#' @export
#'
#' @examples
#' make_pgrid(100)
#' make_pgrid(100, 1, FALSE) #uniform grid including 0 and 1
#' @importFrom utils head tail
#' @importFrom stats pbeta
make_pgrid <- function(n=50L, s=2L, trim=TRUE){
  if(trim) n <- n + 2L
  x<- seq(0, 1, length.out=n)
  res <- stats::pbeta(x, s, s)
  if(trim) return(utils::tail(utils::head(res, -1),-1))
  res
}

#' @param tier integer number of tiers in the linear grid. Each tier contains `tail` share of the previous tier. Default is 3.
#' @param tail real number representing share of grid in each tier. Default is 0.25
#'
#' @rdname make_grid
#' @export
#'
#' @examples
#' make_tgrid(100,3,0.1)
#' @importFrom utils head tail
#' @importFrom stats pbeta
make_tgrid <- function(n=50L, tier=3L, tail=0.25){
  from <- 0
  nn <- split_int(n%/%2, tier)
  res <- vector(mode="list", length = tier)
  for(i in rev(seq_len(tier))){
    a <- seq_grid(from, min(tail^(i-1),0.5), nn[i])
    from <- tail^i
    res[[i]] <- c(a, 1-a)
  }
  if(is.odd(n)) {res <- c(0.5, res) }
  sort(unlist(res))
}


