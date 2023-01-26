## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(BWQS)
library(rstan)
library(MASS)
library(knitr)
library(clusterGeneration)
library(kableExtra)

## ----example_bwqs, echo=FALSE-------------------------------------------------
# fix the seed
set.seed(1990)

# Sample sizes
N = 1000

# Mean & SD of variables
mu = c(0,0,0,1,3)
sd = c(1,1,1,3,1)
sd2 = sd %*% t(sd)

# Correlation Matrix
rho = 0.65
corMat = cbind(c(1,rho,rho^2,rho^2,rho^2),
               c(rho,1,rho^2, rho^2, rho^2),
               c(rho^2,rho^2,1,rho^2,rho^2),
               c(rho^2,rho^2,rho^2,1, rho),
               c(rho^2,rho^2,rho^2,rho,1))

# Covariance Matrix
Sigma = sd2*corMat

# Simulating five correlated exposure variables (Metals)
X = as.data.frame(mvrnorm(N, mu = mu, Sigma = Sigma, empirical = TRUE))
colnames(X) = c('Al','As','Cd','Pb','V')

# Quantile trasformation
Xq = as.matrix(quantile_split(X, mix_name = colnames(X), q=4))

# Intercept coefficient
beta0 = -2

# Overall effect
beta1 = -0.8

# Weights
W = c(0.5,0.20,0.2,0.05,0.05)

# sigma of the model
sigma = 1

# Simulate covariates
D = cbind(rnorm(N,1.3,1), rnorm(N,0.5,0.5))
colnames(D) = c("cov1","cov2")

# coefficient of covariates
delta = c(0.5, 1.1)

# Outcome simulation (continuos)
y = rnorm(N, beta0 + beta1*(Xq %*% W) + D %*% delta, sd = sigma)

# Aggregate data in a data.frame
Data = as.data.frame(cbind(y,X,D))

head(round(Data,3))

## -----------------------------------------------------------------------------
# we run the model ans save results in "fit_bwqs" variable
fit_bwqs = bwqs(y ~ cov1 + cov2, mix_name = c('Al','As','Cd','Pb','V'),
                data = Data, q = 4, family = "gaussian")

## ----tables-------------------------------------------------------------------
fit_bwqs$summary_fit %>%
  kbl() %>%
  kable_classic(full_width = F, html_font = "Cambria")

## ----plots--------------------------------------------------------------------
bwqs_plot(fit_bwqs, parms = "W", size = 2)

## ----plots1-------------------------------------------------------------------
bwqs_plot(fit_bwqs, parms = c("beta0","beta1",'delta',"sigma"), size = 2)

## ----waic---------------------------------------------------------------------
bwqs_waic(fit_bwqs$fit)

## ----hbwqs_dataset, echo=FALSE------------------------------------------------
set.seed(1990)
N = c(250,150,100,100,200)
C = 4
names_ch = paste0("Ch_",1:C)
# Mean and variance of each chemical
mu <- rnorm(C)
sd <- runif(C,0,0.2)
sd2 <- sd %*% t(sd)

# Correlation Matrix
corMat = rcorrmatrix(C,alphad=1)

# Simulating three correlated exposure varMatrix
Sigma <- sd2*corMat
X <- as.data.frame(mvrnorm(sum(N), mu=mu, Sigma = Sigma, empirical=TRUE))
colnames(X)<-names_ch

# Quantile transformation
Xq = matrix(NA,sum(N),C)
for(i in 1:C) Xq[,i] = ecdf(X[[names_ch[i]]])(X[[names_ch[i]]])*4

beta0_1 = 0.65           # intercept 1st cohort
beta1_1 = (-1.7)         # slope 1st cohort
beta0_2 = -0.38          # intercept 2nd cohort
beta1_2 = -0.94          # slope 2nd cohort
beta0_3 = 0.15           # intercept 3rd cohort
beta1_3 = -0.11          # slope 3rd cohort
beta0_4 = (-0.98)        # intercept 4th cohort
beta1_4 = -1.45          # slope 4th cohort
beta0_5 = 1.13           # intercept 5th cohort
beta1_5 = -1.88          # slope 5th cohort
W = c(0.1,0.4,0.3,0.2)   # weights
sigma = 1                # sigma

K = cbind(runif(sum(N),20,40),
        rbinom(sum(N),1,prob = 0.5),
        runif(sum(N),11,18))
colnames(K) = c("m_age","sex","education")
d = c(-0.4,0.5,0.2)

Nh = cumsum(N)
y1 = rnorm(N[1], beta0_1 + beta1_1 * (Xq %*% W)[1:Nh[1]] + (K %*% d)[1:Nh[1]],sigma)
y2 = rnorm(N[2], beta0_2 + beta1_2 * (Xq %*% W)[(Nh[1]+1):Nh[2]] + (K %*% d)[(Nh[1]+1):Nh[2]],sigma)
y3 = rnorm(N[3], beta0_3 + beta1_3 * (Xq %*% W)[(Nh[2]+1):Nh[3]] + (K %*% d)[(Nh[2]+1):Nh[3]],sigma)
y4 = rnorm(N[4], beta0_4 + beta1_4 * (Xq %*% W)[(Nh[3]+1):Nh[4]] + (K %*% d)[(Nh[3]+1):Nh[4]],sigma)
y5 = rnorm(N[5], beta0_5 + beta1_5 * (Xq %*% W)[(Nh[4]+1):Nh[5]] + (K %*% d)[(Nh[4]+1):Nh[5]],sigma)

dt = data.frame(y = c(y1,y2,y3,y4,y5), K, X,
                cohort = rep(c(1,2,3,4,5),N))
head(round(dt,3))

## -----------------------------------------------------------------------------
# we run the model ans save results in "fit_bwqs" variable
fit_hbwqs = bwqs_r(y ~ m_age + sex + education, 
                   mix_name = c("Ch_1","Ch_2","Ch_3","Ch_4"),
                   cluster_name = "cohort", Dalp = rep(1,4),
                   iter = 10000, thin = 3, data = dt, 
                   q = 4, family = "gaussian")

## ----hbwqs_results------------------------------------------------------------
fit_hbwqs$summary_fit %>%
  kbl() %>%
  kable_classic(full_width = F, html_font = "Cambria")

