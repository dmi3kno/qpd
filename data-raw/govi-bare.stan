functions{
  real govindarajulu_qf_s_cdf(real p, real alpha, real sigma) {
    real inv_alpha = inv(alpha);
    // govindarajulu quantile function
    if (alpha<0)
      reject("alpha<0; found alpha =", alpha);
    if (sigma<0)
      reject("sigma<0; found sigma =", sigma);
  return sigma*alpha*pow(p, alpha)*(1+inv_alpha-p);
  }

  real govindarajulu_qf_v_cdf(vector p, real alpha, real sigma) {
    int N = rows(p);
    real inv_alpha=inv(alpha);
    // govindarajulu quantile function
    if (alpha<0)
      reject("alpha<0; found alpha =", alpha);
    if (sigma<0)
      reject("sigma<0; found sigma =", sigma);
   return exp(N*log(sigma)+N*log(alpha)+alpha*sum(log(p))+sum(log1p(inv_alpha-p)));
  }
  real govindarajulu_qf_s_lcdf(real p, real alpha, real sigma) {
    real inv_alpha=inv(alpha);
    // govindarajulu quantile function
    if (alpha<0)
      reject("alpha<0; found alpha =", alpha);
    if (sigma<0)
      reject("sigma<0; found sigma =", sigma);
  return log(sigma)+log(alpha)+alpha*log(p)+log(1+inv_alpha-p);
  }

  real govindarajulu_qf_v_lcdf(vector p, real alpha, real sigma) {
    int N = rows(p);
    real inv_alpha=inv(alpha);
    // govindarajulu quantile function
    if (alpha<0)
      reject("alpha<0; found alpha =", alpha);
    if (sigma<0)
      reject("sigma<0; found sigma =", sigma);
   return N*log(sigma)+N*log(alpha)+alpha*sum(log(p))+sum(log1p(inv_alpha-p));
  }
  real govindarajulu_qdf_s_pdf(real p, real alpha, real sigma) {
    // govindarajulu quantile function
    if (alpha<0)
      reject("alpha<0; found alpha =", alpha);
    if (sigma<0)
      reject("sigma<0; found sigma =", sigma);
  return sigma*alpha*(alpha-1)+pow(p, alpha)*(1-p);
  }

  real govindarajulu_qdf_v_pdf(vector p, real alpha, real sigma) {
    int N = rows(p);
    // govindarajulu quantile function
    if (alpha<0)
      reject("alpha<0; found alpha =", alpha);
    if (sigma<0)
      reject("sigma<0; found sigma =", sigma);
   return exp(N*log(sigma)+N*log(alpha)+N*log(alpha-1)+alpha*sum(log(p))+sum(log1m(p)));
  }
  real govindarajulu_qdf_s_lpdf(real p, real alpha, real sigma) {
    // govindarajulu quantile function
    if (alpha<0)
      reject("alpha<0; found alpha =", alpha);
    if (sigma<0)
      reject("sigma<0; found sigma =", sigma);
  return log(sigma)+log(alpha)+log(alpha-1)+alpha*log(p)+log(1-p);
  }
  real govindarajulu_qdf_v_lpdf(vector p, real alpha, real sigma) {
    int N = rows(p);
    real inv_alpha=inv(alpha);
    vector[N] pow_p_alpha;
    // govindarajulu quantile function
    if (alpha<0)
      reject("alpha<0; found alpha =", alpha);
    if (sigma<0)
      reject("sigma<0; found sigma =", sigma);
   return N*log(sigma)+N*log(alpha)+N*log(alpha-1)+alpha*sum(log(p))+sum(log1m(p));
  }
//////// generalized exponential CDF ////////////
  real genexp_s_cdf(real x, real alpha, real lambda) {
    // generalized exponential cdf
    if (alpha<=0)
      reject("alpha<=0; found alpha =", alpha);
    if (lambda<=0)
      reject("lambda<=0; found lambda =", lambda);
  return (1-exp(-lambda*x))^alpha;
  }
  real genexp_v_cdf(vector x, real alpha, real lambda)  {
    int N = rows(x);
    real inv_alpha=inv(alpha);
    // govindarajulu quantile function
    if (alpha<=0)
      reject("alpha<=0; found alpha =", alpha);
    if (lambda<=0)
      reject("lambda<=0; found lambda =", lambda);
   return exp(alpha*sum(log1m(exp(-lambda*x))));
  }
//////// generalized exponential QF function ////////////
 real genexp_qf_s_cdf(real p, real alpha, real lambda) {
    real inv_lambda = inv(lambda);
    real inv_alpha =inv(alpha);
    // generalized exponential quantile function
    if (alpha<=0)
      reject("alpha<=0; found alpha =", alpha);
    if (lambda<=0)
      reject("lambda<=0; found lambda =", lambda);
  return -inv_lambda*log1m(pow(p,inv_alpha));
  }
  real genexp_qf_v_cdf(vector p, real alpha, real lambda)  {
    int N = rows(p);
    real inv_lambda = inv(lambda);
    real inv_alpha =inv(alpha);
    // generalized exponential quantile function
    if (alpha<=0)
      reject("alpha<=0; found alpha =", alpha);
    if (lambda<=0)
      reject("lambda<=0; found lambda =", lambda);
  return -inv_lambda*log1m_exp(inv_alpha*sum(log(p)));
  }
//////// generalized exponential LOG QF function ////////////
  real genexp_qf_s_lcdf(real p, real alpha, real lambda) {
    real inv_lambda = inv(lambda);
    real inv_alpha =inv(alpha);
    // generalized exponential quantile function
    if (alpha<=0)
      reject("alpha<=0; found alpha =", alpha);
    if (lambda<=0)
      reject("lambda<=0; found lambda =", lambda);
  return log(-inv_lambda)+log(log1m(pow(p,inv_alpha)));
  }
  real genexp_qf_v_lcdf(vector p, real alpha, real lambda)  {
    int N = rows(p);
    real inv_lambda = inv(lambda);
    real inv_alpha =inv(alpha);
    // generalized exponential quantile function
    if (alpha<=0)
      reject("alpha<=0; found alpha =", alpha);
    if (lambda<=0)
      reject("lambda<=0; found lambda =", lambda);
  return log(-inv_lambda)*log(log1m_exp(inv_alpha*sum(log(p))));
  }
//////// generalized exponential QDF function /
  real genexp_qdf_s_pdf(real p, real alpha, real lambda) {
    real inv_alpha =inv(alpha);
    // generalized exponential quantile density function
    if (alpha<=0)
      reject("alpha<=0; found alpha =", alpha);
    if (lambda<=0)
      reject("lambda<=0; found lambda =", lambda);
  return pow(p,inv_alpha-1)/(alpha*lambda*(1-pow(p, inv_alpha)));
  }
  real genexp_qdf_v_pdf(vector p, real alpha, real lambda)  {
    int N = rows(p);
    real inv_alpha =inv(alpha);
    // generalized exponential quantile density function
    if (alpha<=0)
      reject("alpha<=0; found alpha =", alpha);
    if (lambda<=0)
      reject("lambda<=0; found lambda =", lambda);
  return exp((inv_alpha-1)*sum(log(p))) - log(alpha)-log(lambda)-log1m(inv_alpha*sum(log(p)));
  }
//////// generalized exponential LOG QDF function ////////////
  real genexp_qdf_s_lpdf(real p, real alpha, real lambda) {
    real inv_alpha =inv(alpha);
    // generalized exponential quantile density function
    if (alpha<=0)
      reject("alpha<=0; found alpha =", alpha);
    if (lambda<=0)
      reject("lambda<=0; found lambda =", lambda);
  return pow(p,inv_alpha-1)/(alpha*lambda*(1-pow(p, inv_alpha)));
  }
  real genexp_qdf_v_lpdf(vector p, real alpha, real lambda)  {
    int N = rows(p);
    real inv_alpha =inv(alpha);
    // generalized exponential quantile density function
    if (alpha<=0)
      reject("alpha<=0; found alpha =", alpha);
    if (lambda<=0)
      reject("lambda<=0; found lambda =", lambda);
  return exp((inv_alpha-1)*sum(log(p))) - log(alpha)-log(lambda)-sum(log1m(inv_alpha*log(p)));
  }

  vector govindarajulu_root_d(vector u0, vector theta, real[] x_r, int[] x_i){
  // goal-seeking function for passing the parameters as data
    vector[1] res;
    real x = x_r[1];
    real sg = x_r[2];
    real a = x_r[3];
    real tmp = u0[1] + (x - govindarajulu_qf_s_cdf(u0[1], a, sg))/govindarajulu_qdf_s_pdf(u0[1], a, sg);
    res[1] = govindarajulu_qf_s_cdf(tmp, a, sg)-x;
    return(res);
  }

}
