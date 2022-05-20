model_bwqs_regression <- "data {
int<lower=0> N;          // number of individual
int<lower=0> C;          // number of chemicals
matrix[N,C] X;		       // matrix of indipendent variable
vector[C] Dalp;          // vector of the Dirichlet coefficients
real y[N];               // outcome continuos variable
}
parameters {
real beta0;              // intercepts
real beta1;              // overall effect
simplex[C] W;            // weights
real<lower=0> sigma;     // standard deviation of the model
}
transformed parameters {
vector[N] Xb;
Xb = beta0 + beta1*(X*W);
}
model {
beta0 ~ normal(0, 100);
beta1 ~ normal(0, 100);
W ~ dirichlet(Dalp);
sigma ~ inv_gamma(0.01, 0.01);

y ~ normal(Xb, sigma);
}
generated quantities {
vector[N] log_lik;
for (j in 1:N){
log_lik[j] = normal_lpdf(y[j]| Xb[j], sigma);
}
}
"

model_bwqs_regression_positive <- "data {
int<lower=0> N;          // number of individual
int<lower=0> C;          // number of chemicals
matrix[N,C] X;		       // matrix of indipendent variable
vector[C] Dalp;          // vector of the Dirichlet coefficients
real y[N];               // outcome continuos variable
}
parameters {
real beta0;              // intercepts
real beta1;              // overall effect
simplex[C] W;            // weights
real<lower=0> sigma;     // standard deviation of the model
}
transformed parameters {
vector[N] Xb;
Xb = beta0 + beta1*(X*W);
}
model {
beta0 ~ normal(0, 100);
beta1 ~ normal(0, 100) T[0,];
W ~ dirichlet(Dalp);
sigma ~ inv_gamma(0.01, 0.01);

y ~ normal(Xb, sigma);

}
generated quantities {
vector[N] log_lik;
for (j in 1:N){
log_lik[j] = normal_lpdf(y[j]| Xb[j], sigma);
}
}
"

model_bwqs_regression_negative <- "data {
int<lower=0> N;          // number of individual
int<lower=0> C;          // number of chemicals
matrix[N,C] X;		       // matrix of indipendent variable
vector[C] Dalp;          // vector of the Dirichlet coefficients
real y[N];               // outcome continuos variable
}
parameters {
real beta0;              // intercepts
real beta1;              // overall effect
simplex[C] W;            // weights
real<lower=0> sigma;     // standard deviation of the model
}
transformed parameters {
vector[N] Xb;
Xb = beta0 + beta1*(X*W);
}
model {
beta0 ~ normal(0, 100);
beta1 ~ normal(0, 100) T[,0];
W ~ dirichlet(Dalp);
sigma ~ inv_gamma(0.01, 0.01);

y ~ normal(Xb, sigma);
}
generated quantities {
vector[N] log_lik;
for (j in 1:N){
log_lik[j] = normal_lpdf(y[j]| Xb[j], sigma);
}
}
"

model_bwqs_regression_cov <- "data {
int<lower=0> N;          // number of individual
int<lower=0> C;          // number of chemicals
int<lower=0> K;          // number of covariates
matrix[N,C] X;		       // matrix of indipendent variable
matrix[N,K] KV;		       // matrix of covariates
vector[C] Dalp;          // vector of the Dirichlet coefficients
real y[N];               // outcome continuos variable
}
parameters {
real beta0;              // intercepts
real beta1;              // overall effect
vector[K] delta;         // covariates coefficients
simplex[C] W;            // weights
real<lower=0> sigma;     // standard deviation of the model
}
transformed parameters {
vector[N] Xb;
Xb = beta0 + beta1*(X*W) + KV*delta;
}
model {
beta0 ~ normal(0, 100);
beta1 ~ normal(0, 100);
delta ~ normal(0, 100);
W ~ dirichlet(Dalp);
sigma ~ inv_gamma(0.01, 0.01);

y ~ normal(Xb, sigma);
}
generated quantities {
vector[N] log_lik;
for (nn in 1:N)
log_lik[nn] = normal_lpdf(y[nn]| Xb[nn], sigma);

}
"

model_bwqs_regression_cov_positive <- "data {
int<lower=0> N;          // number of individual
int<lower=0> C;          // number of chemicals
int<lower=0> K;          // number of covariates
matrix[N,C] X;		       // matrix of indipendent variable
matrix[N,K] KV;		       // matrix of covariates
vector[C] Dalp;          // vector of the Dirichlet coefficients
real y[N];               // outcome continuos variable
}
parameters {
real beta0;              // intercepts
real beta1;              // overall effect
vector[K] delta;         // covariates coefficients
simplex[C] W;            // weights
real<lower=0> sigma;     // standard deviation of the model
}
transformed parameters {
vector[N] Xb;
Xb = beta0 + beta1*(X*W) + KV*delta;
}
model {
beta0 ~ normal(0, 100);
beta1 ~ normal(0, 100) T[0,];
delta ~ normal(0, 100);
W ~ dirichlet(Dalp);
sigma ~ inv_gamma(0.01, 0.01);

y ~ normal(Xb, sigma);
}
generated quantities {
vector[N] log_lik;
for (nn in 1:N)
log_lik[nn] = normal_lpdf(y[nn]| Xb[nn], sigma);

}
"

model_bwqs_regression_cov_negative <- "data {
int<lower=0> N;          // number of individual
int<lower=0> C;          // number of chemicals
int<lower=0> K;          // number of covariates
matrix[N,C] X;		       // matrix of indipendent variable
matrix[N,K] KV;		       // matrix of covariates
vector[C] Dalp;          // vector of the Dirichlet coefficients
real y[N];               // outcome continuos variable
}
parameters {
real beta0;              // intercepts
real beta1;              // overall effect
vector[K] delta;         // covariates coefficients
simplex[C] W;            // weights
real<lower=0> sigma;     // standard deviation of the model
}
transformed parameters {
vector[N] Xb;
Xb = beta0 + beta1*(X*W) + KV*delta;
}
model {
beta0 ~ normal(0, 100);
beta1 ~ normal(0, 100) T[,0];
delta ~ normal(0, 100);
W ~ dirichlet(Dalp);
sigma ~ inv_gamma(0.01, 0.01);

y ~ normal(Xb, sigma);
}
generated quantities {
vector[N] log_lik;
for (nn in 1:N)
log_lik[nn] = normal_lpdf(y[nn]| Xb[nn], sigma);

}
"

model_bwqs_logit_cov <- "data {
int<lower=0> N;                    // number of individual
int<lower=0> C;                    // number of chemicals
int<lower=0> K;                    // number of covariates
matrix[N,C] X;	                   // matrix of indipendent variable at time 1
matrix[N,K] KV;	                   // matrix of covariates
vector[C] Dalp;                    // vector of the Dirichlet coefficients
int<lower=0, upper=1> y[N];        // outcome binomial variable
}
parameters {
real beta0;              // intercepts
real beta1;              // overall effect
vector[K] delta;         // covariates coefficients
simplex[C] W;            // weights
}
transformed parameters {
vector[N] Xb;
Xb = beta0 + beta1*(X*W) + KV*delta;
}
model {
beta0 ~ normal(0, 100);
beta1 ~ normal(0, 100);
for(j in 1:K) delta[j] ~ normal(0,100);
W ~ dirichlet(Dalp);

y ~ bernoulli_logit(Xb);
}
generated quantities {
vector[N] log_lik;
for (nn in 1:N)
log_lik[nn] = bernoulli_logit_lpmf(y[nn] | Xb[nn]);

}
"

model_bwqs_logit_cov_positive <- "data {
int<lower=0> N;                    // number of individual
int<lower=0> C;                    // number of chemicals
int<lower=0> K;                    // number of covariates
matrix[N,C] X;	                   // matrix of indipendent variable at time 1
matrix[N,K] KV;	                   // matrix of covariates
vector[C] Dalp;                    // vector of the Dirichlet coefficients
int<lower=0, upper=1> y[N];        // outcome binomial variable
}
parameters {
real beta0;              // intercepts
real beta1;              // overall effect
vector[K] delta;         // covariates coefficients
simplex[C] W;            // weights
}
transformed parameters {
vector[N] Xb;
Xb = beta0 + beta1*(X*W) + KV*delta;
}
model {
beta0 ~ normal(0, 100);
beta1 ~ normal(0, 100) T[0,];
for(j in 1:K) delta[j] ~ normal(0,100);
W ~ dirichlet(Dalp);

y ~ bernoulli_logit(Xb);
}
generated quantities {
vector[N] log_lik;
for (nn in 1:N)
log_lik[nn] = bernoulli_logit_lpmf(y[nn] | Xb[nn]);

}
"

model_bwqs_logit_cov_negative <- "data {
int<lower=0> N;                    // number of individual
int<lower=0> C;                    // number of chemicals
int<lower=0> K;                    // number of covariates
matrix[N,C] X;	                   // matrix of indipendent variable at time 1
matrix[N,K] KV;	                   // matrix of covariates
vector[C] Dalp;                    // vector of the Dirichlet coefficients
int<lower=0, upper=1> y[N];        // outcome binomial variable
}
parameters {
real beta0;              // intercepts
real beta1;              // overall effect
vector[K] delta;         // covariates coefficients
simplex[C] W;            // weights
}
transformed parameters {
vector[N] Xb;
Xb = beta0 + beta1*(X*W) + KV*delta;
}
model {
beta0 ~ normal(0, 100);
beta1 ~ normal(0, 100) T[,0];
for(j in 1:K) delta[j] ~ normal(0,100);
W ~ dirichlet(Dalp);

y ~ bernoulli_logit(Xb);
}
generated quantities {
vector[N] log_lik;
for (nn in 1:N)
log_lik[nn] = bernoulli_logit_lpmf(y[nn] | Xb[nn]);

}
"

model_bwqs_logit <- "data {
int<lower=0> N;                    // number of individual
int<lower=0> C;                    // number of chemicals
matrix[N,C] X;	                   // matrix of indipendent variable at time 1
vector[C] Dalp;                    // vector of the Dirichlet coefficients
int<lower=0, upper=1> y[N];        // outcome binomial variable
}
parameters {
real beta0;              // intercepts
real beta1;              // overall effect
simplex[C] W;            // weights
}
transformed parameters {
vector[N] Xb;
Xb = beta0 + beta1*(X*W);
}
model {
beta0 ~ normal(0, 100);
beta1 ~ normal(0, 100);
W ~ dirichlet(Dalp);

y ~ bernoulli_logit(Xb);
}
generated quantities {
vector[N] log_lik;
for (nn in 1:N)
log_lik[nn] = bernoulli_logit_lpmf(y[nn] | Xb[nn]);

}
"

model_bwqs_logit_positive <- "data {
int<lower=0> N;                    // number of individual
int<lower=0> C;                    // number of chemicals
matrix[N,C] X;	                   // matrix of indipendent variable at time 1
vector[C] Dalp;                    // vector of the Dirichlet coefficients
int<lower=0, upper=1> y[N];        // outcome binomial variable
}
parameters {
real beta0;              // intercepts
real beta1;              // overall effect
simplex[C] W;            // weights
}
transformed parameters {
vector[N] Xb;
Xb = beta0 + beta1*(X*W);
}
model {
beta0 ~ normal(0, 100);
beta1 ~ normal(0, 100) T[0,];
W ~ dirichlet(Dalp);

y ~ bernoulli_logit(Xb);
}
generated quantities {
vector[N] log_lik;
for (nn in 1:N)
log_lik[nn] = bernoulli_logit_lpmf(y[nn] | Xb[nn]);

}
"


model_bwqs_logit_negative <- "data {
int<lower=0> N;                    // number of individual
int<lower=0> C;                    // number of chemicals
matrix[N,C] X;	                   // matrix of indipendent variable at time 1
vector[C] Dalp;                    // vector of the Dirichlet coefficients
int<lower=0, upper=1> y[N];        // outcome binomial variable
}
parameters {
real beta0;              // intercepts
real beta1;              // overall effect
simplex[C] W;            // weights
}
transformed parameters {
vector[N] Xb;
Xb = beta0 + beta1*(X*W);
}
model {
beta0 ~ normal(0, 100);
beta1 ~ normal(0, 100) T[,0];
W ~ dirichlet(Dalp);

y ~ bernoulli_logit(Xb);
}
generated quantities {
vector[N] log_lik;
for (nn in 1:N)
log_lik[nn] = bernoulli_logit_lpmf(y[nn] | Xb[nn]);

}
"

model_bwqs_poisson_cov <- "data {
int<lower=0> N;                // number of individual
int<lower=0> C;                // number of chemicals
int<lower=0> K;                // number of covariates
matrix[N,C] X;                 // matrix of indipendent variable
matrix[N,K] KV;                // matrix of covariates
vector[C] Dalp;                // vector of the Dirichlet coefficients
int<lower=0>  y[N];                     // outcome categorical variable
}
parameters {
real beta0;              // intercepts
real beta1;              // overall effect time 1
vector[K] delta;         // overall effect time 2
simplex[C] W;            // weights
}
transformed parameters {
vector[N] Xb;
vector[N] mu;
Xb = beta0 + beta1*(X*W) + KV*delta;
for(j in 1:N) mu[j] = exp(Xb[j]);
}
model {
beta0 ~ normal(0, 100);
beta1 ~ normal(0, 100);
delta ~ normal(0, 100);
W ~ dirichlet(Dalp);

y ~ poisson(mu);
}
generated quantities {
vector[N] log_lik;
for (j in 1:N){
log_lik[j] = poisson_lpmf(y[j] | mu[j]);
}
}
"

model_bwqs_poisson_cov_positive <- "data {
int<lower=0> N;                // number of individual
int<lower=0> C;                // number of chemicals
int<lower=0> K;                // number of covariates
matrix[N,C] X;                 // matrix of indipendent variable
matrix[N,K] KV;                // matrix of covariates
vector[C] Dalp;                // vector of the Dirichlet coefficients
int<lower=0>  y[N];                     // outcome categorical variable
}
parameters {
real beta0;              // intercepts
real beta1;              // overall effect time 1
vector[K] delta;         // overall effect time 2
simplex[C] W;            // weights
}
transformed parameters {
vector[N] Xb;
vector[N] mu;
Xb = beta0 + beta1*(X*W) + KV*delta;
for(j in 1:N) mu[j] = exp(Xb[j]);
}
model {
beta0 ~ normal(0, 100);
beta1 ~ normal(0, 100) T[0,];
delta ~ normal(0, 100);
W ~ dirichlet(Dalp);

y ~ poisson(mu);
}
generated quantities {
vector[N] log_lik;
for (j in 1:N){
log_lik[j] = poisson_lpmf(y[j] | mu[j]);
}
}
"

model_bwqs_poisson_cov_negative <- "data {
int<lower=0> N;                // number of individual
int<lower=0> C;                // number of chemicals
int<lower=0> K;                // number of covariates
matrix[N,C] X;                 // matrix of indipendent variable
matrix[N,K] KV;                // matrix of covariates
vector[C] Dalp;                // vector of the Dirichlet coefficients
int<lower=0>  y[N];                     // outcome categorical variable
}
parameters {
real beta0;              // intercepts
real beta1;              // overall effect time 1
vector[K] delta;         // overall effect time 2
simplex[C] W;            // weights
}
transformed parameters {
vector[N] Xb;
vector[N] mu;
Xb = beta0 + beta1*(X*W) + KV*delta;
for(j in 1:N) mu[j] = exp(Xb[j]);
}
model {
beta0 ~ normal(0, 100);
beta1 ~ normal(0, 100) T[,0];
delta ~ normal(0, 100);
W ~ dirichlet(Dalp);

y ~ poisson(mu);
}
generated quantities {
vector[N] log_lik;
for (j in 1:N){
log_lik[j] = poisson_lpmf(y[j] | mu[j]);
}
}
"

model_bwqs_poisson <- "data {
int<lower=0> N;                // number of individual
int<lower=0> C;                // number of chemicals
matrix[N,C] X;                 // matrix of indipendent variable
vector[C] Dalp;                // vector of the Dirichlet coefficients
int<lower=0>  y[N];                     // outcome categorical variable
}
parameters {
real beta0;              // intercepts
real beta1;              // overall effect time 1
simplex[C] W;            // weights
}
transformed parameters {
vector[N] Xb;
vector[N] mu;
Xb = beta0 + beta1*(X*W);
for(j in 1:N) mu[j] = exp(Xb[j]);
}
model {
beta0 ~ normal(0, 100);
beta1 ~ normal(0, 100);
W ~ dirichlet(Dalp);

y ~ poisson(mu);
}
generated quantities {
vector[N] log_lik;
for (j in 1:N){
log_lik[j] = poisson_lpmf(y[j] | mu[j]);
}
}
"

model_bwqs_poisson_positive <- "data {
int<lower=0> N;                // number of individual
int<lower=0> C;                // number of chemicals
matrix[N,C] X;                 // matrix of indipendent variable
vector[C] Dalp;                // vector of the Dirichlet coefficients
int<lower=0>  y[N];                     // outcome categorical variable
}
parameters {
real beta0;              // intercepts
real beta1;              // overall effect time 1
simplex[C] W;            // weights
}
transformed parameters {
vector[N] Xb;
vector[N] mu;
Xb = beta0 + beta1*(X*W);
for(j in 1:N) mu[j] = exp(Xb[j]);
}
model {
beta0 ~ normal(0, 100);
beta1 ~ normal(0, 100) T[0,];
W ~ dirichlet(Dalp);

y ~ poisson(mu);
}
generated quantities {
vector[N] log_lik;
for (j in 1:N){
log_lik[j] = poisson_lpmf(y[j] | mu[j]);
}
}
"

model_bwqs_poisson_negative <- "data {
int<lower=0> N;                // number of individual
int<lower=0> C;                // number of chemicals
matrix[N,C] X;                 // matrix of indipendent variable
vector[C] Dalp;                // vector of the Dirichlet coefficients
int<lower=0>  y[N];                     // outcome categorical variable
}
parameters {
real beta0;              // intercepts
real beta1;              // overall effect time 1
simplex[C] W;            // weights
}
transformed parameters {
vector[N] Xb;
vector[N] mu;
Xb = beta0 + beta1*(X*W);
for(j in 1:N) mu[j] = exp(Xb[j]);
}
model {
beta0 ~ normal(0, 100);
beta1 ~ normal(0, 100) T[,0];
W ~ dirichlet(Dalp);

y ~ poisson(mu);
}
generated quantities {
vector[N] log_lik;
for (j in 1:N){
log_lik[j] = poisson_lpmf(y[j] | mu[j]);
}
}
"


model_rbqws_regression_cov <-  "data {
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

model_rbqws_regression <- "data {
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
