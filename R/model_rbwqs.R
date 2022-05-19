model_rbqws_regression_cov =  "
data {
  int<lower=0> N;
  int<lower=0> J;
  int<lower=0> C;
  int<lower=0> K;
  vector[N] y;
  matrix[N,K] X;
  matrix[N,C] Chem;
  int cohort[N];
  vector[C] Dalp;
}
parameters {
  simplex[C] W;
  real<lower=0> sigma;
  real<lower=0> sigma_a;
  real<lower=0> sigma_b;
  vector[J] a;
  vector[J] b;
  vector[K] delta;
  real mu_a;
  real mu_b;
}
transformed parameters {
  vector[N] BWQS;
  vector[N] dX;
  vector[N] mu;
  BWQS = Chem*W;
  dX = X*delta;
  mu = a[cohort] + b[cohort].*BWQS + dX;
}
model {
  mu_a ~ normal(0, 100);
  mu_b ~ normal(0, 100);
  sigma_a ~ inv_gamma(0.01,0.01);
  sigma_b ~ inv_gamma(0.01,0.01);
  W ~ dirichlet(Dalp);  
  sigma ~ inv_gamma(0.01,0.01);
  delta ~ normal(0,100); 

  a ~ normal(mu_a, sigma_a);
  b ~ normal(mu_b, sigma_b);
  y ~ normal(mu, sigma);
}
generated quantities {
vector[N] log_lik;
for (nn in 1:N)
log_lik[nn] = normal_lpdf(y[nn]| mu[nn], sigma);

}
"

model_rbqws_regression = "
data {
  int<lower=0> N;
  int<lower=0> J;
  int<lower=0> C;
  vector[N] y;
  matrix[N,C] Chem;
  int cohort[N];
  vector[C] Dalp;
}
parameters {
  simplex[C] W;
  real<lower=0> sigma;
  real<lower=0> sigma_a;
  real<lower=0> sigma_b;
  vector[J] a;
  vector[J] b;
  real mu_a;
  real mu_b;
}
transformed parameters {
  vector[N] BWQS;
  vector[N] mu;
  BWQS = Chem*W;
  mu = a[cohort] + b[cohort].*BWQS;
}
model {
  mu_a ~ normal(0, 100);
  mu_b ~ normal(0, 100);
  sigma_a ~ inv_gamma(0.01,0.01);
  sigma_b ~ inv_gamma(0.01,0.01);
  W ~ dirichlet(Dalp);  
  sigma ~ inv_gamma(0.01,0.01);

  a ~ normal(mu_a, sigma_a);
  b ~ normal(mu_b, sigma_b);
  y ~ normal(mu, sigma);
}
generated quantities {
vector[N] log_lik;
for (nn in 1:N)
log_lik[nn] = normal_lpdf(y[nn]| mu[nn], sigma);

}
"