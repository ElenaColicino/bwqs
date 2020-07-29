model_2bwqs_regression <- "data {
int<lower=0> N;           // number of individual
int<lower=0> C1;          // number of element of first mix
int<lower=0> C2;          // number of element of second mix
matrix[N,C1] XC1;	        // matrix of first mix
matrix[N,C2] XC2;	        // matrix of second mix
vector[C1] DalpC1;        // vector of the Dirichlet coefficients for first mix
vector[C2] DalpC2;        // vector of the Dirichlet coefficients for second mix
real y[N];                // outcome continuos variable
}
parameters {
real beta0;               // intercepts
real beta1;               // overall effect of first mix
real beta2;               // overall effect of second mix
simplex[C1] WC1;          // weights of first mix
simplex[C2] WC2;          // weights of second mix
real<lower=0> sigma;      // standard deviation of the model
}
transformed parameters {
vector[N] Xb;
Xb = beta0 + beta1*(XC1*WC1) + beta2*(XC2*WC2);
}
model {
beta0 ~ normal(0, 100);
beta1 ~ normal(0, 100);
beta2 ~ normal(0, 100);
WC1 ~ dirichlet(DalpC1);
WC2 ~ dirichlet(DalpC2);
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

model_2bwqs_regression_cov <- "data {
int<lower=0> N;          // number of individual
int<lower=0> C1;         // number of elements in the first mix
int<lower=0> C2;         // number of elements in the second mix
int<lower=0> K;          // number of covariates
matrix[N,C1] XC1;	       // matrix of elements in the first mix
matrix[N,C2] XC2;	       // matrix of elements in the second mix
matrix[N,K] KV;	         // matrix of covariates
vector[C1] DalpC1;       // vector of the Dirichlet coefficients for first mix
vector[C2] DalpC2;       // vector of the Dirichlet coefficients for second mix
real y[N];               // outcome continuos variable
}
parameters {
real beta0;              // intercepts
real beta1;              // overall effect of first mix
real beta2;              // overall effect of second mix
vector[K] delta;         // covariates coefficients
simplex[C1] WC1;         // weights of first mix
simplex[C2] WC2;         // weights of second mix
real<lower=0> sigma;     // standard deviation of the model
}
transformed parameters {
vector[N] Xb;
Xb = beta0 + beta1*(XC1*WC1) + beta2*(XC2*WC2) + KV*delta;
}
model {
beta0 ~ normal(0, 100);
beta1 ~ normal(0, 100);
beta2 ~ normal(0, 100);
for(j in 1:K) delta[j] ~ normal(0,100);
WC1 ~ dirichlet(DalpC1);
WC2 ~ dirichlet(DalpC2);
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

model_2bwqs_logit <- "data {
int<lower=0> N;              // number of individual
int<lower=0> C1;             // number of element of first mix
int<lower=0> C2;             // number of element of second mix
matrix[N,C1] XC1;	           // matrix of first mix
matrix[N,C2] XC2;	           // matrix of second mix
vector[C1] DalpC1;           // vector of the Dirichlet coefficients for first mix
vector[C2] DalpC2;           // vector of the Dirichlet coefficients for second mix
int<lower=0, upper=1> y[N];  // outcome binomial variable
}
parameters {
real beta0;               // intercepts
real beta1;               // overall effect of first mix
real beta2;               // overall effect of second mix
simplex[C1] WC1;          // weights of first mix
simplex[C2] WC2;          // weights of second mix
}
transformed parameters {
vector[N] Xb;
Xb = beta0 + beta1*(XC1*WC1) + beta2*(XC2*WC2);
}
model {
beta0 ~ normal(0, 100);
beta1 ~ normal(0, 100);
beta2 ~ normal(0, 100);
WC1 ~ dirichlet(DalpC1);
WC2 ~ dirichlet(DalpC2);

y ~ bernoulli_logit(Xb);
}
generated quantities {
vector[N] log_lik;
for (nn in 1:N)
log_lik[nn] = bernoulli_logit_lpmf(y[nn] | Xb[nn]);
}
"

model_2bwqs_logit_cov <- "data {
int<lower=0> N;              // number of individual
int<lower=0> C1;             // number of elements in the first mix
int<lower=0> C2;             // number of elements in the second mix
int<lower=0> K;              // number of covariates
matrix[N,C1] XC1;	           // matrix of elements in the first mix
matrix[N,C2] XC2;	           // matrix of elements in the second mix
matrix[N,K] KV;	             // matrix of covariates
vector[C1] DalpC1;           // vector of the Dirichlet coefficients for first mix
vector[C2] DalpC2;           // vector of the Dirichlet coefficients for second mix
int<lower=0, upper=1> y[N];  // outcome continuos variable
}
parameters {
real beta0;                  // intercepts
real beta1;                  // overall effect of first mix
real beta2;                  // overall effect of second mix
vector[K] delta;             // covariates coefficients
simplex[C1] WC1;             // weights of first mix
simplex[C2] WC2;             // weights of second mix
}
transformed parameters {
vector[N] Xb;
Xb = beta0 + beta1*(XC1*WC1) + beta2*(XC2*WC2) + KV*delta;
}
model {
beta0 ~ normal(0, 100);
beta1 ~ normal(0, 100);
beta2 ~ normal(0, 100);
for(j in 1:K) delta[j] ~ normal(0,100);
WC1 ~ dirichlet(DalpC1);
WC2 ~ dirichlet(DalpC2);

y ~ bernoulli_logit(Xb);
}
generated quantities {
vector[N] log_lik;
for (nn in 1:N)
log_lik[nn] = bernoulli_logit_lpmf(y[nn] | Xb[nn]);
}
"
