## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(BWQS)
library(MASS)
library(knitr)

## ----example_bwqs-------------------------------------------------------------
# fix the seed
set.seed(1990)

# Sample sizes
N <- 1000

# Mean & SD of variables
mu <- c(0,0,0,1,3)
sd <- c(1,1,1,3,1)
sd2 <- sd %*% t(sd)

# Correlation Matrix
rho <- 0.65
corMat <- cbind(c(1,rho,rho^2,rho^2,rho^2),
                c(rho,1,rho^2, rho^2, rho^2),
                c(rho^2,rho^2,1,rho^2,rho^2),
                c(rho^2,rho^2,rho^2,1, rho),
                c(rho^2,rho^2,rho^2,rho,1))

# Covariance Matrix
Sigma <- sd2*corMat

# Simulating five correlated exposure variables (Metals)
X <- as.data.frame(mvrnorm(N, mu=mu, Sigma = Sigma, empirical=TRUE))
colnames(X)<-c('Al','As','Cd','Pb','V')

# Quantile trasformation
Xq <- as.matrix(quantile_split(X, mix_name = colnames(X), q=4))

# Intercept coefficient
beta0 <- -2

# Overall effect
beta1 <- -0.8

# Weights
W <- c(0.5,0.20,0.2,0.05,0.05)

# sigma of the model
sigma <- 1

# Simulate covariates
D = cbind(rnorm(N,1.3,1), rnorm(N,0.5,0.5))
colnames(D) = c("cov1","cov2")

# coefficient of covariates
delta = c(0.5, 1.1)

# Outcome simulation (continuos)
y <- rnorm(N, beta0 + beta1*(Xq %*% W) + D %*% delta, sd = sigma)

# Aggregate data in a data.frame
Data <-as.data.frame(cbind(y,X,D))

# we run the model ans save results in "fit_bwqs" variable
fit_bwqs <- bwqs(y ~ cov1 + cov2, mix_name = c('Al','As','Cd','Pb','V'),
                 data = Data, q=4, prior = "negative", family = "gaussian")

## ----tables-------------------------------------------------------------------
fit_bwqs$summary_fit

## ----plots--------------------------------------------------------------------
a = bwqs_plot(fit_bwqs, parms = "W", size = 2)
a

## ----plots1-------------------------------------------------------------------
b = bwqs_plot(fit_bwqs, parms = c("beta0","beta1",'delta',"sigma"), size = 2)
b

## ----waic---------------------------------------------------------------------
bwqs_waic(fit_bwqs$fit)

