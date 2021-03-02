#include govi-bare.stan

data{
  int N;
  vector[N] xs;
  real sg;
  real a;
}

transformed data{
  vector[1] y_guess = [0.01]';
  vector[N] u;
  vector[1] tmp;
  vector[0] theta; // empty parameter vector
  int x_i[0]; // empty integer data vector
  real x_r[3];
  x_r[2] = sg;
  x_r[3] = a;

  for (i in 1:N){
    x_r[1] = xs[i];
    tmp = algebra_solver(govindarajulu_root_d, y_guess, theta, x_r, x_i);
    u[i] = tmp[1];
  }
}

parameters{
  real r;
}

model{
  //target += genexp_qdf_s_lpdf(a | theta, 1);
  for(i in 1:N){
    print(u[i]);
  }
  //  target += govindarajulu_qdf_s_lpdf(u[i] | a, sigma);
  r ~ normal(0,1);
}
