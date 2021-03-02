functions{ // begin the functions block
vector fit_metalog(vector p, vector quantiles) {
    // code by Ben Goodrich <benjamin.goodrich@columbia.edu> 
    int n = rows(p);
    matrix[n, n] Y;
    vector[n] log_odds = n > 1 ? logit(p) : rep_vector(not_a_number(), n);
    vector[n] pmhalf =   n > 2 ? p - 0.5  : rep_vector(not_a_number(), n);
    int odd = 1;
    if (n == 0) reject("p cannot be of size zero");
    if (rows(quantiles) != n) reject("p and quantiles must be of the same size");
  
    Y[ , 1] = rep_vector(1, n);
    if (n > 1) Y[ , 2] = log_odds;
    if (n > 2) Y[ , 3] = pmhalf .* log_odds;
    if (n > 3) Y[ , 4] = pmhalf;
    for (m in 5:n) {
      if (odd) {
          pmhalf .*= pmhalf;
          Y[ , m]  = pmhalf;
        } else Y[ , m] = pmhalf .* log_odds;
        odd = odd == 0;
    }
    return Y \ quantiles;
  }
real metalog_qf_s_cdf(real y, vector a){
    // code by Dmytro Perepolkin <dperepolkin@gmail.com>
    int n = rows(a);
    real res=a[1];
    real logity = logit(y);
    real ymhalf = y-0.5;
    int odd = 1;
    if (n == 0) reject("a cannot be of size zero");
    if (n > 1) res += a[2]*logity;
    if (n > 2) res += a[3]*ymhalf*logity;
    if (n > 3) res += a[4]*ymhalf;
    for (m in 5:n) {
      if (odd) {
          res += a[m]*pow(ymhalf, (m-1)/2.0);
        } else res += a[m]*pow(ymhalf, m/2.0-1)*logity;
        odd = odd == 0;
    }
    return res;
}

real metalog_lqf_s_lcdf(real y, vector a){
  return log(metalog_qf_s_cdf(y, a));
}

real metalog_rng(vector a) {
    return metalog_qf_s_cdf(uniform_rng(0, 1), a);
}
real metalog_qdf_s_pdf(real y, vector a){
    // code by Dmytro Perepolkin <dperepolkin@gmail.com>
    int n = rows(a);
    real res;
    real yt1my = y*(1-y);
    real logity = logit(y);
    real ymhalf = y-0.5;
    int odd = 1;
    if (n == 0) reject("a cannot be of size zero");
    if (n > 1) res = a[2]/yt1my;
    if (n > 2) res += a[3]*(ymhalf/yt1my+logity);
    if (n > 3) res += a[4];
    for (m in 5:n) {
      if (odd) {
          res += a[m]*(m-1)/2.0*pow(ymhalf, (m-3)/2.0);
        } else res += a[m]*(pow(ymhalf, m/2.0-1)/yt1my+(m/2.0-1)*pow(ymhalf, m/2.0-2)*logity);
        odd = odd == 0;
    }
    return res;
 }

 real metalog_lqdf_s_lpdf(real y, vector a){
    // code by Dmytro Perepolkin <dperepolkin@gmail.com>
    return log(metalog_qdf_s_pdf(y, a));
 }

 real metalog_rqdf_s_pdf(real y, vector a){
    // code by Dmytro Perepolkin <dperepolkin@gmail.com>
    return inv(metalog_qdf_s_pdf(y, a));
 }

 real metalog_lrqdf_s_lpdf(real y, vector a){
    // code by Dmytro Perepolkin <dperepolkin@gmail.com>
    return log(inv(metalog_qdf_s_pdf(y, a)));
 }


 vector metalog_approx_u(vector x, vector a, real tol){
  // goal-seeking function for passing the parameters as data
   int N=rows(x);
   vector[N] u_guess = cumulative_sum(rep_vector(1, N))/(N+1); // 1:N/(N+1)
   vector[N] u;
   real u0;
   for(i in 1:N){
     u0 = u_guess[i];
     while(fabs(metalog_qf_s_cdf(u0, a)- x[i]) > tol){
       u0 += (x[i]-metalog_qf_s_cdf(u0, a))/metalog_qdf_s_pdf(u0, a);
       u0 = fmin(fmax(u0, tol), 1-tol);
     }
     u[i] = u0;
   }
   return u;
  }
} // finish the functions block
data{
  int N; // number of observations
  vector[N] x; //data
  vector[4] alpha; // dirichlet parameters
  int idx[4];  // how dirichlet was built
  vector[3] p_est; // estimated probabilities
  vector[3] qntls; // estimated quantiles
}

transformed data{
  vector[N] xs = sort_asc(x); // this is to make sure data comes in sorted
  vector[3] a0 = fit_metalog(p_est, qntls);
  vector[N] u;
  // fill the data vector
  u = metalog_approx_u(xs, a0, 1e-06);
}

parameters{
  simplex[4] theta;
}

transformed parameters {
  vector[4] cumtheta;
  vector[3] a;
  cumtheta = cumulative_sum(theta[idx]);
  a = fit_metalog(cumtheta[1:3], qntls);
  //a =fit_metalog(p_est, qntls); // for test without prior
}
model{
  target += dirichlet_lpdf(theta | alpha);
  for (i in 1:N){
    target += metalog_lrqdf_s_lpdf(u[i] | a);
  }
}

generated quantities{
  simplex[4] pp_theta = dirichlet_rng(alpha);
  vector[4] pp_cumtheta = cumulative_sum(pp_theta[idx]);
  vector[3] pp_a = fit_metalog(pp_cumtheta[1:3], qntls);
  real y_sim = metalog_rng(pp_a);

}
