library(qpd)
library(winference)
library(gridExtra)

metalog_reparametrize <- function(thetas, idx, qntls){
  if(!inherits(thetas, "matrix") )
    thetas <- matrix(thetas, nrow=1, byrow = TRUE)
  N <- nrow(thetas)
  ps <- isim_to_cprob(thetas, idx)
  q_m <- matrix(rep.int(qntls, N), nrow=N, byrow=TRUE)
  res <- sapply(seq.int(N), function(i) fit_metalog(ps[i,], q_m[i,]))
  t(res)
}

rprior <- function(N, parameters){
  #particles <- matrix(N, ncol=4)
  particles <- qpd::rdir(N, parameters$k)
  particles[,-3]
}

logdprior <- function(thetaparticles, parameters){
  thetas <- cbind(thetaparticles[,1], thetaparticles[,2], 1-rowSums(thetaparticles), thetaparticles[,3])
  logdensities <- qpd::ddir(thetas, parameters$k)
  log(logdensities)
}

xdqfmetalog <- function(xs, theta, idx, q_star){
  cat("thetas: ", theta, "\n")
  theta <- c(theta, 1-sum(theta))[idx]
  p <- isim_to_cprob(theta, idx)
  cat("ps: ", p, "\n")
  a <- fit_metalog(p, q_star)
  if(!is_metalog_valid(a)) return(NA_real_)
  ys <- pmetalog(xs, a)
  res <- qpd::fmetalog(ys, a)
  log(1/res)
}

loglikelihood <- function(thetas, xs, parameters, ...){
  idx <- parameters$idx
  q_star <- parameters$q_star
  sapply(seq.int(nrow(thetas)), function(i) sum(xdqfmetalog(xs, thetas[i,], idx, q_star), na.rm = TRUE))
}

simulate <- function(theta, parameters, ...){
  idx <- parameters$idx
  q_star <- parameters$q_star
  a <- metalog_reparametrize(theta, idx, q_star)[1,]
  qmetalog(runif(nobs), a)
}

nobs <- 100
#obs <- qpd::rmetalog(nobs, a_star)
idx <- c(1L, 2L, 4L, 3L)
# hyperparameter vector for Qdirichlet
k <- c(3.14, 10.6, 2.24, 8.67)
# true values of quantiles
q_star <- c(5,9,12)
# true p (and, as a result, theta) for data generation
p_star <- c(0.1, 0.4, 0.8)
theta_star <- diff(c(0, p_star, 1))[idx]
a_star <- qpd::fit_metalog(p_star, q_star) # not needed
# check
all.equal(cumsum(theta_star[idx])[-length(theta_star)], p_star)

target <- list(rprior=rprior, dprior=logdprior, simulate=simulate, loglikelihood=loglikelihood,
               parameter_names=c("k", "idx", "q_star"), parameters=list(k=k, idx=idx, q_star=q_star),
               thetadim=3, ydim=1)

obs <- target$simulate(theta_star, target$parameters)

#curve(sapply(x, FUN=function(u) qmetalog(u, a_star)), from = 1e-10, to=1-1e-10, n=500, ylab="quantile")

#test: sample from distribution
#samples <- qmetalog(runif(1e4), a_star)
#hist(samples, prob=TRUE, nclass=1000)

#theta <- rprior(100, parameters)
tuning_params <- list(niterations=10000, nchains=1, cov_proposal=diag(1e-4, 3, 3),
                      adaptation=2500, init_chains=matrix(theta_star[-3], nrow=1))
mh <- metropolishastings(obs, target, tuning_params)
chains_df <- mhchainlist_to_dataframe(mh$chains)
burnin <- 2500
g1 <- ggplot(chains_df %>% filter(iteration>burnin, iteration%%10==1),
            aes(x=iteration, y=X.1, group=ichain))+geom_line()+
  geom_hline(yintercept = theta_star[1], col="red", linetype=2)
g1

g2 <- ggplot(chains_df %>% filter(iteration>burnin, iteration%%10==1),
             aes(x=iteration, y=X.2, group=ichain))+geom_line()+
  geom_hline(yintercept = theta_star[2], col="red", linetype=2)
g2

g3 <- ggplot(chains_df %>% filter(iteration>burnin, iteration%%10==1),
             aes(x=iteration, y=X.3, group=ichain))+geom_line()+
  geom_hline(yintercept = theta_star[4], col="red", linetype=2)
g3


grid.arrange(g1, g2, g3, ncol=2)

d1 <- ggplot(chains_df %>% filter(iteration>burnin),
       aes(x=X.1, fill=as.factor(ichain), group=as.factor(ichain)))+geom_density(alpha=0.1)+
  geom_vline(xintercept = theta_star[1], col="red", linetype=2)
d2 <- ggplot(chains_df %>% filter(iteration>burnin),
             aes(x=X.2, fill=as.factor(ichain), group=as.factor(ichain)))+geom_density(alpha=0.1)+
  geom_vline(xintercept = theta_star[2], col="red", linetype=2)
d3 <- ggplot(chains_df %>% filter(iteration>burnin),
             aes(x=X.3, fill=as.factor(ichain), group=as.factor(ichain)))+geom_density(alpha=0.1)+
  geom_vline(xintercept = theta_star[4], col="red", linetype=2)

grid.arrange(d1, d2, d3, ncol=2)

obs_sorted<-sort(obs)
compute_d <-  function(y_fake){
  y_fake_sorted <- sort(y_fake)
  return(mean(abs(obs_sorted-y_fake_sorted)))
}
param_algo <- list(nthetas=1024, nmoves=1, proposal=mixture_rmixmod(),
                   minimum_diversity=0.5, R=2, maxtrials=1e5)

results <- wsmc(compute_d, target, param_algo, maxtime=3*60)
results$thetas_history
