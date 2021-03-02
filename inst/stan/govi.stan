//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//
functions {
  // for lifetime model theta is set to zero
  // y is value, returns probability
  real govindarajulu_qf_lcdf(vector p, real sg, real a){
  int N = rows(p);
  return(-N*log1p(sg)- a * sum(log1p(p)) - N*log1p(a) - sum(log1p(1+1/a-p)));
  }

  // y is value, returns probability
  real govindarajulu_qdf_lpdf(vector p, real sg, real a){
  int N = rows(p);
  return(-N*log1p(sg)-N*log1p(a)-N*log1p(a+1)-a*sum(log1p(p))-sum(log1p(1-p)) );
  }

  real genexp_lcdf(vector y, real theta){
    // F(y)=(1-exp(-y))^theta
    return sum(theta*log1m_exp(-y));
  }

  real genexp_qf_lcdf(vector p, real theta){
    int N = rows(p);
    real inv_theta = inv(theta);
    vector[N] p_inv_theta = exp(inv_theta*log(p));
    return -sum(log(-log1m(p_inv_theta)));
  }

  real genexp_qdf_lpdf(vector p, real theta){
  // generalized exponential CDF: F(alpha)=(1-exp(-alpha))^theta
  // generalized exponential QF: Q(p)=-log(1-p^(1/theta))
  //generalized exponential QDF by differentiation:
  //d/dp(-log(1 - p^(1/theta))) = p^(1/theta - 1)/(theta - theta*p^(1/theta))
  // we need to implement negative log of QDF
  //
  int N = rows(p);
  real inv_theta=inv(theta);
  vector[N] p_inv_theta;
   for(i in 1:N){
     p_inv_theta[i]=pow(p[i], inv_theta);
   }
  return -sum(log(p)*theta/(theta-1)+log(theta)+log1m(p_inv_theta));
 }

}

data{
  int<lower=0> N;
  vector[N] y;
  real <lower=0> sg;
  real theta;
}

parameters{
  real <lower=0> a;
}

model{
  a ~ genexp_qdf(theta);
  for(i in 1:N){
  y ~ govindarajulu_qdf(sg, a[i]);
  }
}
