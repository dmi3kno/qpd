
library(tidyverse)
library(qpd)

## categories. These are conditional probabilities
### Normal - Overweight - Obese - Underweight

elicrec <- tibble::tribble(
	~cat, ~lq, ~mq, ~uq,
	"norm", 0.6012526096033404, 0.6555323590814196, 0.7056367432150313, 
	"over", 0.17954070981210857, 0.20250521920668063, 0.24634655532359084,
	"obese", 0.06889352818371608, 0.08768267223382042, 0.11273486430062629) # a=1.69605752864251 b=1.296


#fitted_ab <- tibble::tribble(~cat, ~alpha, ~beta,
#																													"norm", 24.79463794306823, 13.213042458017421,
#																													"over", 6.479672591559014, 4.3052760039826055,
#																													"obese", 2.7253997552201112, 1.6882363331590782)

# L3+U4=U3+L4=1-0.65-0.20=0.15
#lower = 0.04 = 0.15 - 0.11
#mid = 0.06 = 1-0.65-0.20-0.09
#upper = 0.08 = 0.15 - 0.07

lastrec <- list(
	cat="under",
	lq=1-sum(elicrec$mq[1:2])-elicrec$uq[3],
	mq=1-sum(elicrec$mq[1:3]),
	uq=1-sum(elicrec$mq[1:2])-elicrec$lq[3])

#df_qs <- dplyr::bind_rows(elicrec, lastrec)

normalize_cond_quantilies <- function(lq, mq, uq, 
																		.names=c("lq_norm", "mq_norm", "uq_norm")){
	d_mq=1-lag(cumsum(mq), default = 0)
	setNames(tibble::tibble(col1=lq/d_mq,	col2=mq/d_mq,	col3=uq/d_mq),
										.names)
}


df_qs <- 
	elicrec %>% 
	mutate(normalize_cond_quantilies(lq, mq, uq)) %>% 
	#left_join(fitted_ab, by="cat") %>% 
	mutate(ab=pmap(list(lq_norm, mq_norm, uq_norm), qpd::fit_beta)) %>% 
	unnest(ab) %>%
	bind_rows(tibble(cat="under", alpha=.$beta[nrow(.)], beta=NA_real_)) %>% 
	mutate(beta_corr=lead(rev(cumsum(rev(alpha))),1),
								alf_cum=cumsum(alpha),
								asum=sum(alf_cum)+sum(beta, na.rm=TRUE),
								N=asum/n(),
								mu_rec=ifelse(is.na(beta), alpha/(alf_cum), alpha/(alf_cum + beta)),
								mu_norm=mu_rec/sum(mu_rec),
								alpha_norm=mu_norm*N,
								beta_norm=N-alpha_norm,
								m_lq=qbeta(0.25, shape1 =alpha_norm , shape2 = beta_norm),
								m_mq=qbeta(0.50, shape1 =alpha_norm , shape2 = beta_norm),
								m_uq=qbeta(0.75, shape1 =alpha_norm , shape2 = beta_norm),
								gd_anorm1 =  ifelse(is.na(beta), 1, alpha/(alpha+beta)), # normalizing for Generalized Dirichlet
								gd_anorm2 = ifelse(is.na(beta), 1, alpha*(alpha+1)/((alpha+beta)*(alpha+beta+1))),# normalizing for Generalized Dirichlet
								gd_bnorm1 =  cumprod(beta/(alpha+beta)), # normalizing for Generalized Dirichlet
								gd_bnorm2 = cumprod(beta*(beta+1)/((alpha+beta)*(alpha+beta+1))),# normalizing for Generalized Dirichlet
								gd_S=gd_anorm1*lag(gd_bnorm1,1, default = 1),
								gd_T= gd_anorm2*lag(gd_bnorm2,1, default = 1),
								gd_C=gd_S*(gd_S-gd_T)/(gd_T-gd_S^2),
								gd_D=(1-gd_S)*(gd_S-gd_T)/(gd_T-gd_S^2),
								gm_lq=qbeta(0.25, shape1 = gd_C, shape2 = gd_D),
								gm_mq=qbeta(0.50, shape1 = gd_C, shape2 = gd_D),
								gm_uq=qbeta(0.75, shape1 = gd_C, shape2 = gd_D)
								) # mu reconciled 

df_qs

# N reconciled

#N_rec <- 27.08

#sum(df_qs$mu_rec)

#ar <- double(nrow(df_qs))
#ar[1] <- N_rec * df_qs[["mu_norm"]][1]
#for(i in 2:(nrow(df_qs)-1)){
#	ar[i] <- (N_rec - sum( ar[seq(i-1)] ) ) * df_qs[["mu_norm"]][i]
#}
#ar[nrow(df_qs)] <- N_rec-sum(ar[seq(nrow(df_qs)-1)])
#ar

# feedback for the diriechlet
#br <- N_rec-ar
#
#marg <- tibble::tibble(
#	cat = c("norm", "over", "obese", "under"),
#	lq = qbeta(0.25, ar, br),
#	mq = qbeta(0.5, ar, br),
#	uq = qbeta(0.75, ar, br)
#)
#### Feedback for Dirichlet
# Say expert wants to modify marginals for one category
# quartiles of the first category should be
#marg$lq[1] <- 0.60
#marg$mq[1] <- 0.65
#marg$uq[1] <- 0.72

#upd_marg <- marg %>% 
# mutate(ab=pmap(list(lq, mq, uq), fit_beta)) %>% 
#	unnest(ab) %>% 
#	mutate(N_upd = sum(alpha+beta)/n(),
#		      mu_norm=alpha/(alpha+beta),
#								mu_norm=mu_norm/sum(mu_norm),
#								a_upd=N_upd*mu_norm,
#								b_upd=N_upd-a_upd)

#### Feedback for Connor-Mosimann
# Feedback for Connor-Mosimann distribution
# marginals are not beta, but can be approximated as such



