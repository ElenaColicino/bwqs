#' Fitting Bayesian Weighted Quantile Sum regression models with two mixture
#'
#' Fits Bayesian Weighted Quantile Sum (BWQS) regressions for continuous and binomial outcomes with
#' two different mixture. This model allows two cosider the effect of two different mixture. This
#' approach is useful when different mixture can lead opposite effect on the outcome variable.
#'
#' @param formula Object of class \code{formula} specifying the relationship between the outcome and the
#' covariates of the model not involved in the mixture variable. If the model has no covariates specify
#' \code{y ~ NULL}.
#' @param mix_name_1 A character vector listing the variables contributing to first mixture effect.
#' @param mix_name_2 A character vector listing the variables contributing to second mixture effect.
#' @param data The \code{data.frame} containing the variables (covariates and elements of both of the mixtures)
#' to be included in the model.
#' @param q An \code{integer} to specify how mixture variables will be ranked, e.g. in quartiles
#' (\code{q = 4}), deciles (\code{q = 10}), or percentiles (\code{q = 100}). If \code{q = NULL} then
#' the values of the mixture variables are taken (these must be standardized).
#' @param DalpC1 A \code{vector} containing the parameters of the Dirichlet distribution for the weights
#' of the first mixture. For further details see https://en.wikipedia.org/wiki/Dirichlet_distribution.
#' @param DalpC2 A \code{vector} containing the parameters of the Dirichlet distribution for the weights
#' of the second mixture.
#' @param chains An \code{integer} to specify the number of chain in Hamiltonian Monte Carlo algorithm
#' (1 by default).
#' @param iter An \code{integer} to specify the lenght of chain in Hamiltonian Monte Carlo algorithm
#' (10000 by default).
#' @param thin An \code{integer} to specify the thinning parameter in Hamiltonian Monte Carlo algorithm.
#' @param seed An \code{integer} value to fix the seed. If \code{seed = NULL} the seed are randomly choosen.
#' @param start_value A \code{vector} containing the initial value of the prior distribution,
#' if it is equal to NULL random values are chosen.
#' @param c_int A \code{vector} of two elements to specify the credible intervals for parameters, for 95% credible
#' interval \code{c_int = c(0.025,0.975)} (default).
#' @param family A \code{string} to specify the type of outcome. Possible values are "gaussian" (default) and
#' "binomial".
#'
#'@details
#'The function \code{bwqs_2} uses the package \code{rstan} which allows the connection with STAN,
#'a specific software, written in C++ for bayesian inference, for further information see https://mc-stan.org/.
#'
#' @return
#'  \code{bwqs} returns a list with two argument:
#'  \item{fit}{An \code{S4} object with all details of the Hamiltonian Monte Carlo, all the extractions
#'  from the posterior distribution and all values of the parameters}
#'  \item{summary_fit}{Table with the statistics of the parameters: mean, standard error of the mean,
#'    standard deviation, lower and upper values for the credible interval (with credible level specified
#'    by \code{c_int}), n_eff and Rhat. For further details see https://cran.r-project.org/web/packages/rstan/rstan.pdf}
#'
#' @author
#' Nicolo Foppa Pedretti, Elena Colicino
#'
#' @import rstan
#' @import Rcpp
#' @import methods
#' @import Matrix
#'
#' @examples
#' #load libraries
#' library(MASS)
#' library(BWQS)
#' library(Matrix)
#'
#' # fix the seed
#' set.seed(456)
#'
#' # Sample sizes
#' N <- 800
#'
#' # Mean & SD of variables
#' mu <- c(0,0,0,1,3,8,2,1,3,1)
#' sd <- c(1,1,1,3,1,2,1,4,5,3)
#' sd2 <- sd %*% t(sd)
#'
#' # Correlation Matrix
#' #rho <- 0.65
#' corMat = matrix(runif(100), 10, 10)
#' corMat = (corMat * lower.tri(corMat)) + t(corMat * lower.tri(corMat))
#' diag(corMat) <- 1
#' corMat = nearPD(corMat, posd.tol=1.e-04)$mat
#'
#' # Covariance Matrix
#' Sigma <- sd2*corMat
#'
#' # Simulated dataset
#' X <- as.data.frame(mvrnorm(N, mu=mu, Sigma = Sigma, empirical=TRUE))
#' colnames(X)<-c("X1","X2","X3","X4","X5","X6","X7","X8","X9","X10")
#' #corrplot.mixed(cor(Xq), lower.col = "black", number.cex = .8)
#' # Intercept values
#' beta0<- 2
#'
#' # Overall effect mix 1
#' beta1 <- 1
#'
#' # Overall effect mix 2
#' beta2 <- (-2.8)
#'
#' # Weights
#' W1 <- c(0.5,0.20,0.2,0.05,0.05)
#' W2 <- c(0.05,0.15,0.7,0.05,0.05)
#'
#' # sigma
#' sigma <- 1#sd(beta1*(Xq %*% W))
#'
#' # Quantile extraction
#' Xq <- as.matrix(quantile_split(X, mix_name = colnames(X), q=4))
#'
#' # Outcome simulation (continuos)
#' y <- rnorm(N, beta0 + beta1*(Xq[,1:5] %*% W1) + beta2*(Xq[,6:10] %*% W2), sigma)
#
#' # Aggregate data in a data.frame
#' Data <-as.data.frame(cbind(y,X))
#'
#' fit_2bwqs <- bwqs_2(y ~ NULL,
#'                     mix_name_1 = c("X1","X2","X3","X4","X5"),
#'                     mix_name_2 = c("X6","X7","X8","X9","X10"),
#'                     data = Data, q=4,
#'                     family = "gaussian",
#'                     iter = 10000, thin = 5)
#'
#' fit_2bwqs$summary_fit
#'
#' @export

bwqs_2 <- function(formula, mix_name_1, mix_name_2, data, q,
                   DalpC1 = NULL, DalpC2 = NULL,
                   chains = 1, iter = 10000, thin = 3,
                   seed=2019, start_value=NULL,
                   c_int=c(0.025,0.975), family="gaussian"){

  formula = as.formula(formula)
  y_name  <- all.vars(formula)[1]
  KV_name <- all.vars(formula)[-1]
  X1_name  <- mix_name_1
  X2_name  <- mix_name_2
  seed <- ifelse(is.null(seed), ceiling(runif(1,0,100000)), seed)

  check_input_2(formula, mix_name_1, mix_name_2, data, q, DalpC1,
                DalpC2, chains, iter, thin, seed,
                start_value, c_int, family)

  if(length(KV_name)==0){
    data = as.data.frame(data[,c(y_name,X1_name,X2_name)])
  } else{
    data = as.data.frame(data[,c(y_name,KV_name,X1_name,X2_name)])
  }

  data <- na.omit(data)
  if(nrow(data)==0) stop("No dataset available")

  if(is.null(q)){
    X1 = data[,mix_name_1]
    X2 = data[,mix_name_2]
  } else{
    X1 = quantile_split(data=data,mix_name=mix_name_1,q)[,mix_name_1]
    X2 = quantile_split(data=data,mix_name=mix_name_2,q)[,mix_name_2]
  }

  if(length(KV_name)==0){
    KV = NULL
  } else{
    KV = data[,KV_name]
  }

  if(is.null(DalpC1)){
    DalpC1 = rep(1, length(mix_name_1))
  } else{
    DalpC1 = DalpC1
  }

  if(is.null(DalpC2)){
    DalpC2 = rep(1, length(mix_name_2))
  } else{
    DalpC2 = DalpC2
  }

  if(is.null(start_value)){
    start_value = "random"
  }
  
  suppressWarnings({
  switch (family,
          gaussian = {if(!is.null(KV)){
            data_reg <- list(
              N   = nrow(data),
              C1  = length(mix_name_1),
              C2  = length(mix_name_2),
              K   = length(KV_name),
              XC1 = cbind(X1),
              XC2 = cbind(X2),
              KV  = cbind(KV),
              DalpC1 = DalpC1,
              DalpC2 = DalpC2,
              y = as.vector(data[,y_name])
            )

            fit <- stan(model_code = model_2bwqs_regression_cov,
                        data = data_reg,
                        init = start_value,
                        chains = 1,
                        warmup = iter/2,
                        iter = iter,
                        cores = 1,
                        thin = thin,
                        refresh = 0,
                        algorithm = 'NUTS',
                        seed = seed,
                        control=list(max_treedepth = 20,
                                     adapt_delta = 0.999999999999999))

          } else{
            data_reg <- list(
              N   = nrow(data),
              C1  = length(mix_name_1),
              C2  = length(mix_name_2),
              XC1 = cbind(X1),
              XC2 = cbind(X2),
              DalpC1 = DalpC1,
              DalpC2 = DalpC2,
              y = as.vector(data[,y_name])
            )

            fit <- stan(model_code = model_2bwqs_regression,
                        data = data_reg,
                        init = start_value,
                        chains = 1,
                        warmup = iter/2,
                        iter = iter,
                        cores = 1,
                        thin = thin,
                        refresh = 0,
                        algorithm = 'NUTS',
                        seed = seed,
                        control=list(max_treedepth = 20,
                                     adapt_delta = 0.999999999999999))

          }},

          binomial = {if(!is.null(KV)){
            data_reg <- list(
              N   = nrow(data),
              C1  = length(mix_name_1),
              C2  = length(mix_name_2),
              K   = length(KV_name),
              XC1 = cbind(X1),
              XC2 = cbind(X2),
              KV  = cbind(KV),
              DalpC1 = DalpC1,
              DalpC2 = DalpC2,
              y = as.vector(data[,y_name])
            )

            fit <- stan(model_code = model_2bwqs_logit_cov,
                        data = data_reg,
                        init = start_value,
                        chains = 1,
                        warmup = iter/2,
                        iter = iter,
                        cores = 1,
                        thin = thin,
                        refresh = 0,
                        algorithm = 'NUTS',
                        seed = seed,
                        control=list(max_treedepth = 20,
                                     adapt_delta = 0.999999999999999))

          } else{
            data_reg <- list(
              N   = nrow(data),
              C1  = length(mix_name_1),
              C2  = length(mix_name_2),
              XC1 = cbind(X1),
              XC2 = cbind(X2),
              DalpC1 = DalpC1,
              DalpC2 = DalpC2,
              y = as.vector(data[,y_name])
            )

            fit <- stan(model_code = model_2bwqs_logit,
                        data = data_reg,
                        init = start_value,
                        chains = 1,
                        warmup = iter/2,
                        iter = iter,
                        cores = 1,
                        thin = thin,
                        refresh = 0,
                        algorithm = 'NUTS',
                        seed = seed,
                        control=list(max_treedepth = 20,
                                     adapt_delta = 0.999999999999999))

          }}
  )
  })

  if(length(KV_name)==0){
    if(family=="gaussian"){
      parameters <- c('beta0','beta1','beta2',
                      paste0('W1_',mix_name_1),
                      paste0('W2_',mix_name_2),
                      'sigma')
    } else{
      parameters <- c('beta0','beta1','beta2',
                      paste0('W1_',mix_name_1),
                      paste0('W2_',mix_name_2))
    }

  }else{
    if(family=="gaussian"){
      parameters <- c('beta0','beta1','beta2',paste0('C_',KV_name),
                      paste0('W1_',mix_name_1),paste0('W2_',mix_name_2),
                      'sigma')
    } else{
      parameters <- c('beta0','beta1','beta2',paste0('C_',KV_name),
                      paste0('W1_',mix_name_1),paste0('W2_',mix_name_2))
    }
  }

  names(fit)[1:length(parameters)] <- parameters
  sum_fit <- round(summary(fit,pars = parameters,
                           probs = c_int)$summary,5)
  return(list(fit=fit, summary_fit = sum_fit))
}


