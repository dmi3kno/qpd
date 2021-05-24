
# internal function for normalizing probabilities in the QDirichlet frame
# returns QDirichlet frame
#' @keywords internal
qdframe_normalize_probs <- function(df, id, spt, prb){
  df <- as.data.frame(df)
  orig_cols <- names(df)
  spts <- unique(df[[spt]])
  stopifnot(length(spts)==3L)
  mid_spt <- spts[2]
  # lag=c(NA, head(x,-1))
  lag_mid_prb <- c(0, head(df[order(df[[id]]), prb][df[[spt]]==mid_spt], -1))
  norm_c_df <- data.frame(".norm_c"=1-cumsum(lag_mid_prb))
  norm_c_df[[id]] <- unique(df[order(df[[id]]), id])
  merged_df <- merge(df, norm_c_df, by=id, all.x = TRUE)
  merged_df[[prb]] <- merged_df[[prb]]/merged_df[[".norm_c"]]
  merged_df[,orig_cols]
}

# internal function for fitting beta distribution to every category in QDirichlet frame
# returns ab frame
#' @keywords internal
qdframe_fit_beta <- function(df, id, spt, prb){
  ab_df_lst <- lapply(unique(df[[id]]),
                      function(.x) {v <- df[df[[id]]==.x, prb]
                      ab_df <- do.call(qpd::fit_beta,as.list(v))
                      ab_df[[id]] <- .x
                      ab_df
                      })
  res_df <- do.call(rbind, ab_df_lst)
  # expecting that qpd::fit_beta() returns a dataframe with c("alpha", "beta") as column names
  # also expecting that recently added id column is last
  stopifnot(identical(names(res_df), c("alpha", "beta", id)))

  seq_to_maxid <- seq_len(max(res_df[[id]]))
  if(identical(sort(res_df[[id]]),  seq_to_maxid))
    mssng_id <-  max(res_df[[id]])+1L
  else
    mssng_id <- setdiff(seq_to_maxid, res_df[[id]])
  # order rows and columns
  res <- res_df[order(res_df[[id]]), c(id, setdiff(names(res_df), id))]
  N <- nrow(res)
  mssng_ab <- data.frame(alpha=res[["beta"]][N],
                         beta=NA_real_,
                         stringsAsFactors = FALSE, row.names = NULL)
  mssng_ab[[id]] <- mssng_id
  # reorder columns in mssng_ab dataset
  mssng_ab[ , c(id, setdiff(names(mssng_ab), id))]
  rbind(res, mssng_ab)
}

# internal function for converting abframe into parameters of generalized dirichlet
#' @keywords internal
#' @importFrom stats qbeta setNames
abframe_dir<- function(df, id, alpha, beta, spt){
  a <- df[[alpha]]
  b <- df[[beta]]
  i <- df[[id]]
  n <- nrow(df)
  #cat("a ",a, "\n")
  ac <- rev(cumsum(rev(a)))
  # lead=c(tail(x,-1),NA)
  beta_corr <- c(tail(ac, -1), NA_real_)
  a_cum <- cumsum(a)
  a_sum <- sum(a_cum)+sum(b, na.rm=TRUE)
  N <- a_sum / n
  mu_rec <- ifelse(is.na(b), a/a_cum, a/(a_cum+b))
  mu_norm <- mu_rec/sum(mu_rec)
  a_norm <- mu_norm * N
  b_norm <- N - a_norm
  m <- sapply(spt, stats::qbeta, shape1 = a_norm, shape2 = b_norm)

  dir_res <- stats::setNames(data.frame(a=a_norm, .id=i), c("a", id))
  return(dir_res[,c(id, "a")])
}

# internal function for converting abframe into parameters of generalized dirichlet
#' @keywords internal
#' @importFrom stats qbeta setNames
abframe_gendir<- function(df, id, alpha, beta, spt){
  a <- df[[alpha]]
  b <- df[[beta]]
  i <- df[[id]]
  n <- nrow(df)
  # normalizing for generalized dirichlet
  gd_a_norm1 <- ifelse(is.na(b), 1, a/(a+b))
  gd_a_norm2 <- ifelse(is.na(b), 1, a*(a+1)/((a+b)*(a+b+1)))
  gd_b_norm1 <- cumprod(b/(a+b))
  gd_b_norm2 <- cumprod(b*(b+1)/((a+b)*(a+b+1)))
  gd_S <- gd_a_norm1*c(1, head(gd_b_norm1,-1))
  gd_T <- gd_a_norm2*c(1, head(gd_b_norm2,-1))
  gd_C <- gd_S*(gd_S-gd_T)/(gd_T-gd_S^2)
  gd_D <- (1-gd_S)*(gd_S-gd_T)/(gd_T-gd_S^2)
  gm <- sapply(spt, stats::qbeta, shape1 = gd_C, shape2 = gd_D)

  gendir_res <- stats::setNames(data.frame(a=a, b=b, .id=i), c("a", "b", id))
  return(gendir_res[seq_len(n-1L),c(id, "a", "b")])
}

#' Fit Dirichlet or Generalized Dirichlet distribution to the dataframe of elicited conditional QPPs
#'
#' @param df data frame of conditional QPPs with 3 columns: category id, SPTs (probability triplets per category), values (probability value for Qdirichlet)
#' @param id_col character name of the column containing category id
#' @param spt_col character name of the column containing pribability triplets per category
#' @param prb_col  character name of the column containing values (probabilities for QDirichlet)
#'
#' @rdname elicit
#' @return data frame with 2(3) columns containing category id and the parameters of Dirichlet or
#' Generalized Dirichlet distribution
#' @export
#'
#' @examples
#' df <- data.frame(id=rep(1:3, each=3),
#'                  spt=rep(0.25*1:3, 3),
#'                  values=c(7, 9, 12, 34, 41, 50, 5, 7.5, 15)/100)
#'  fit_dir(df, "id", "spt", "values")
#'  fit_gendir(df, "id", "spt", "values")
fit_dir <- function(df, id_col, spt_col, prb_col){
  spt_v <- unique(df[[spt_col]])
  qf_norm <- qdframe_normalize_probs(df, id=id_col, spt=spt_col, prb=prb_col)
  abframe <- qdframe_fit_beta(qf_norm, id=id_col, spt=spt_col, prb=prb_col)
  abframe_dir(abframe, id=id_col, alpha="alpha", beta="beta", spt=spt_v)
}


#' @rdname elicit
#' @return
#' @export
fit_gendir <- function(df, id_col, spt_col, prb_col){
  spt_v <- unique(df[[spt_col]])
  qf_norm <- qdframe_normalize_probs(df, id=id_col, spt=spt_col, prb=prb_col)
  abframe <- qdframe_fit_beta(qf_norm, id=id_col, spt=spt_col, prb=prb_col)
  abframe_gendir(abframe, id=id_col, alpha="alpha", beta="beta", spt=spt_v)
}
