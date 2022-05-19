#' Fitting Bayesian Weighted Quantile Sum regression models
#'
#' Fits Bayesian Weighted Quantile Sum (BWQS) regressions for continuous and binomial outcomes. This model
#' provides estimation for the mixture composition and overall effect of the mixture on the outcomes using
#' bayesian framework.
#'
#' @param formula Object of class \code{formula} specifying the relationship between the outcome and the
#' covariates of the model not involved in the mixture variable. If the model has no covariates specify
#' \code{y ~ NULL}.
#' @param mix_name A character vector listing the variables contributing to a mixture effect.
#' @param data The \code{data.frame} containing the variables (covariates and elements of the mixture)
#' to be included in the model.
#' @param q An \code{integer} to specify how mixture variables will be ranked, e.g. in quartiles
#' (\code{q = 4}), deciles (\code{q = 10}), or percentiles (\code{q = 100}). If \code{q = NULL} then
#' the values of the mixture variables are taken (these must be standardized or the domain must be the same).
#' @param Dalp A \code{vector} containing the parameters of the Dirichlet distribution of the weights, the number
#' of the elements of the vector has to be equal to the number of chemicals. If \code{Dalp = NULL}
#' the domain is explored uniformly.
#' @param chains An \code{integer} to specify the number of chain in Hamiltonian Monte Carlo algorithm.
#' Default value \code{chains = 1}.
#' @param iter An \code{integer} to specify the lenght of chain in Hamiltonian Monte Carlo algorithm.
#' Default value \code{iter = 10000}.
#' @param thin An \code{integer} to specify the thinning parameter in Hamiltonian Monte Carlo algorithm.
#' @param seed An \code{integer} value to fix the seed. If \code{seed = NULL} the seed are randomly choosen.
#' @param start_value A \code{vector} containing the initial value of the prior distribution,
#' if it is equal to NULL random values are chosen.
#' @param c_int A \code{vector} of two elements to specify the credible intervals for parameters, for 95% credible
#' interval \code{c_int = c(0.025,0.975)} (default).
#' @param family A \code{string} to specify the type of outcome. Possible values are "gaussian" (default),
#' "binomial" and "poisson".
#' @param prior A \code{string} to specify the direction of prior distribution. Possible values are
#' "None"(default), "positive" and "negative".
#'
#' @details
#' The function \code{bwqs} uses the package \code{rstan} which allows the connection with STAN,
#' a specific software, written in C++ for bayesian inference, for further information see https://mc-stan.org/.
#'
#' @return
#' \code{bwqs} returns a list with two argument:
#' \item{fit}{An \code{S4} object with all details of the Hamiltonian Monte Carlo, all the extractions
#' from the posterior distribution and all values of the parameters}
#' \item{summary_fit}{Table with the statistics of the parameters: mean, standard error of the mean,
#' standard deviation, lower and upper values for the credible interval (with credible level specified
#' by \code{c_int}), n_eff and Rhat. For further details see https://cran.r-project.org/web/packages/rstan/rstan.pdf}
#'
#' @author
#' Nicolo Foppa Pedretti, Elena Colicino
#'
#' @import rstan
#' @import Rcpp
#' @import methods
#'
#' @export

bwqs <- function(formula, mix_name, data, q, Dalp = NULL,
                 chains = 1, iter = 10000, thin = 3,
                 seed=2019, start_value=NULL,
                 c_int=c(0.025,0.975), family="gaussian", prior="None"){

  formula = as.formula(formula)
  y_name  <- all.vars(formula)[1]
  KV_name <- all.vars(formula)[-1]
  X_name  <- mix_name
  seed <- ifelse(is.null(seed), ceiling(runif(1,0,100000)), seed)

  check_input(formula, mix_name, data, q, Dalp,
              chains, iter, thin, seed,
              start_value, c_int, family, prior)

  if(length(KV_name)==0){
    data = as.data.frame(data[,c(y_name,X_name)])
  } else{
    data = as.data.frame(data[,c(y_name,KV_name,X_name)])
  }

  data <- na.omit(data)
  if(nrow(data)==0) stop("No dataset available")

  if(is.null(q)){
    X = data[,mix_name]
  } else{
    X = quantile_split(data=data,mix_name=mix_name,q)[,mix_name]
  }

  if(length(KV_name)==0){
    KV = NULL
  } else{
    KV = data[,KV_name]
  }

  if(is.null(Dalp)){
    Dalp = rep(1, length(mix_name))
  } else{
    Dalp = Dalp
  }

  if(is.null(start_value)){
    start_value = "random"
  }

  suppressWarnings({
  switch (family,
          gaussian = {if(!is.null(KV)){

            data_reg <- list(
              N = nrow(data),
              C = length(mix_name),
              K = length(KV_name),
              X = cbind(X),
              KV = cbind(KV),
              Dalp = Dalp,
              y = as.vector(data[,y_name])
            )

            if(prior == "None"){

            fit <- stan(model_code = model_bwqs_regression_cov,
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

              if(prior == "positive"){

                fit <- stan(model_code = model_bwqs_regression_cov_positive,
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

                fit <- stan(model_code = model_bwqs_regression_cov_negative,
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

              }
            }

          } else{
            data_reg <- list(
              N = nrow(data),
              C = length(mix_name),
              X = X,
              Dalp = Dalp,
              y = as.vector(data[,y_name])
            )

            if(prior == "None"){

              fit <- stan(model_code = model_bwqs_regression,
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

              if(prior == "positive"){

                fit <- stan(model_code = model_bwqs_regression_positive,
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

                fit <- stan(model_code = model_bwqs_regression_negative,
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
              }
            }

          }},

          binomial = {if(!is.null(KV)){
            data_reg <- list(
              N = nrow(data),
              C = length(mix_name),
              K = length(KV_name),
              X = cbind(X),
              KV = cbind(KV),
              Dalp = Dalp,
              y = as.vector(data[,y_name])
            )

            if(prior == "None"){

              fit <- stan(model_code = model_bwqs_logit_cov,
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

              if(prior == "positive"){

                fit <- stan(model_code = model_bwqs_logit_cov_positive,
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

                fit <- stan(model_code = model_bwqs_logit_cov_negative,
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

              }
            }


          } else{
            data_reg <- list(
              N = nrow(data),
              C = length(mix_name),
              X = X,
              Dalp = Dalp,
              y = as.vector(data[,y_name])
            )

            if(prior == "None"){

              fit <- stan(model_code = model_bwqs_logit,
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

              if(prior == "positive"){

                fit <- stan(model_code = model_bwqs_logit_positive,
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

                fit <- stan(model_code = model_bwqs_logit_negative,
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
              }
            }


          }},

          poisson = {if(!is.null(KV)){
            data_reg <- list(
              N = nrow(data),
              C = length(mix_name),
              K = length(KV_name),
              X = cbind(X),
              KV = cbind(KV),
              Dalp = Dalp,
              y = as.vector(data[,y_name])
            )

            if(prior == "None"){

              fit <- stan(model_code = model_bwqs_poisson_cov,
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

              if(prior == "positive"){

                fit <- stan(model_code = model_bwqs_poisson_cov_positive,
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

                fit <- stan(model_code = model_bwqs_poisson_cov_negative,
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

              }
            }


          } else{
            data_reg <- list(
              N = nrow(data),
              C = length(mix_name),
              X = X,
              Dalp = Dalp,
              y = as.vector(data[,y_name])
            )

            if(prior == "None"){

              fit <- stan(model_code = model_bwqs_poisson,
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

              if(prior == "positive"){

                fit <- stan(model_code = model_bwqs_poisson_positive,
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

                fit <- stan(model_code = model_bwqs_poisson_negative,
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
              }
            }


          }}

  )})

  if(length(KV_name)==0){
    if(family=="gaussian"){
      parameters <- c('beta0','beta1',
                      paste0('W_',mix_name),'sigma')
    } else{
      parameters <- c('beta0','beta1',
                      paste0('W_',mix_name))
    }

  }else{
    if(family=="gaussian"){
      parameters <- c('beta0','beta1',paste0('C_',KV_name),
                      paste0('W_',mix_name),'sigma')
    } else{
      parameters <- c('beta0','beta1',paste0('C_',KV_name),
                      paste0('W_',mix_name))
    }
  }

  names(fit)[1:length(parameters)] <- parameters
  sum_fit <- round(summary(fit,pars = parameters,
                           probs = c_int)$summary,5)
  return(list(fit=fit, summary_fit = sum_fit))
}


#' Compute of Watanabe-Akaike Information Criterion
#'
#' This function compute the WAIC and LOO given a posterior distribution.
#' The values of the logarithmic likelihood are extracted from Stan object
#' of the bwqs function.
#'
#' @param stanfit Object \code{fit} is returned from the \code{bwqs} function.
#'
#' @return The function returns a list with:
#' \item{waic}{value of total WAIC for the model}
#' \item{elpd_waic}{expected log pointwise predictive density for a new dataset}
#' \item{p_waic}{the estimated effective number of parameters}
#' \item{elpd_loo}{approximation to the actual out-of-sample prediction error under the model}
#' \item{p_loo}{effective number of parameters as
#' the bias adjustment corresponding to the overfitting inherent}
#'
#' @details This function is implemented based on the formulas given
#' in http://www.stat.columbia.edu/~gelman/research/unpublished/waic_stan.
#' For examples of this function see example of \code{bwqs} function
#'
#' @export
#'

bwqs_waic <- function(stanfit){

  log_lik <- extract(stanfit, "log_lik")$log_lik

  dim(log_lik) <- if (length(dim(log_lik))==1) c(length(log_lik),1) else
    c(dim(log_lik)[1], prod(dim(log_lik)[2:length(dim(log_lik))]))

  S <- nrow(log_lik)
  n <- ncol(log_lik)
  lpd <- log(colMeans(exp(log_lik)))
  p_waic <- colVars(log_lik)
  elpd_waic <- lpd - p_waic
  waic <- -2*elpd_waic

  loo_weights_raw <- 1/exp(log_lik-max(log_lik))
  loo_weights_normalized <- loo_weights_raw/
    matrix(colMeans(loo_weights_raw),nrow=S,ncol=n,byrow=TRUE)
  loo_weights_regularized <- pmin (loo_weights_normalized, sqrt(S))
  elpd_loo <- log(colMeans(exp(log_lik)*loo_weights_regularized)/
                    colMeans(loo_weights_regularized))
  p_loo <- lpd - elpd_loo
  pointwise <- cbind(waic,lpd,p_waic,elpd_waic,p_loo,elpd_loo)
  total <- colSums(pointwise)
  se <- sqrt(n*colVars(pointwise))
  return(list(waic=total["waic"], elpd_waic=total["elpd_waic"],
              p_waic=total["p_waic"], elpd_loo=total["elpd_loo"], p_loo=total["p_loo"]))
}


#' Plot function for Bayesian Weighted Quantile Sum regression models
#'
#' This function provide a plot of the parameters with a credible interval. The paramenters of
#' interest could be manually selected and also the level of credible interval.
#'
#' @param stanfit the object coming from \code{bwqs()} function.
#' @param parms a vector with the name of parameters that we want to plot possible value are
#' \code{'beta0'}, \code{'beta1'}, \code{'W'}, \code{'delta'} (components of covariates) and
#' \code{'sigma'} (if \code{'gaussian'} family is selected).
#' @param cri_level An \code{integer} containing the credible level for the interval. This number must
#' be between 0 and 1.
#' @param point_color color of point that indicates the value for the fitted parameters (default \code{black}).
#' @param line_color color of line that indicates the credible interval (default \code{black}).
#' @param shape symbols to indicate the value of parameters.
#' @param size dimension of \code{shape}.
#'
#' @import rstan
#' @importFrom stats as.formula na.omit quantile runif
#' @importFrom ggplot2 geom_point aes_string
#' @details For examples of this function see example of \code{bwqs} function
#' @export

bwqs_plot = function(stanfit, parms = NULL, cri_level = 0.95, point_color="black", line_color="black", shape=21, size=5){
  if(is.null(parms)) stop("bwqs_plot: Parameters not decleared")
  summary_r = summary(stanfit$fit,pars=parms)$summary
  pointdata = data.frame(xname = summary_r[rownames(summary_r),1],
                         ypos = seq(length(rownames(summary_r)),1))
  suppressMessages({
  stan_plot(stanfit$fit, pars = parms, point_est = "mean",
            show_density = F, show_outer_line = FALSE, ci_level = cri_level,
            fill_color = line_color) + geom_point(data = pointdata,
                                                  aes_string("xname", "ypos"),color = point_color,
                                                  shape = shape, size = size, fill = point_color)
  })
}


#' Quantile function for BWQS
#'
#' This function allows quantile splitting for BWQS model.
#' @param data A \code{data.frame} where the columns have the name of
#' the mixture components.
#' @param mix_name A character list of components that have been
#' considered in the mixture.
#' @param q An \code{integer} to specify how mixture variables will be standardize, e.g. in quartiles
#' (\code{q = 4}), deciles (\code{q = 10}), or percentiles (\code{q = 100}).
#' @param shift If \code{SHITF = TRUE} the quantiles are counted from 0, otherwise 1.
#' @return A \code{data.frame} with quantiles for element of the list only.
#' @details For examples of this function see example of \code{bwqs} function
#' @export
#'

quantile_split <- function(data, mix_name = mix_name, q, shift=TRUE){

  if(shift){
    for (i in 1:length(mix_name)){
      dat_num = as.numeric(unlist(data[, mix_name[i]]))
      data[[mix_name[i]]] = cut(dat_num,
                                breaks = unique(quantile(dat_num, probs = seq(0, 1, by = 1/q), na.rm = TRUE)),
                                labels = FALSE,
                                include.lowest = TRUE)-1
    }

  } else{
    for (i in 1:length(mix_name)){
      dat_num = as.numeric(unlist(data[, mix_name[i]]))
      data[[mix_name[i]]] = cut(dat_num, breaks = unique(quantile(dat_num, probs = seq(0, 1, by = 1/q), na.rm = TRUE)),
                                labels = FALSE, include.lowest = TRUE)
    }

  }

  return(data)
}

#' Quantile function for BWQS with continuous data
#'
#' This function allows quantile splitting (empirical cumulative function)
#' for BWQS model.
#' @param data A \code{data.frame} where the columns have the name of
#' the mixture components.
#' @param mix_name A character list of components that have been
#' considered in the mixture.
#' @param q An \code{integer} to specify how mixture variables will be standardize, e.g. in quartiles
#' (\code{q = 4}), deciles (\code{q = 10}), or percentiles (\code{q = 100}).
#' @return A \code{data.frame} with specified quantiles for element of the list only.
#' @details For examples of this function see example of \code{bwqs} function
#' @export
#'

quantile_split_2 <- function(data, mix_name = mix_name, q){
  for(i in mix_name) data[,i] = ecdf(data[,i])(data[,i])*q
  return(data)
}

#' Fitting Bayesian Weighted Quantile Sum regression models
#'
#' Fits Random Bayesian Weighted Quantile Sum (BWQS) regressions for continuous outcomes. This model
#' provides estimation for the mixture composition and overall effect of the mixture across different
#' groups on the outcomes using bayesian framework.
#'
#' @param formula Object of class \code{formula} specifying the relationship between the outcome and the
#' covariates of the model not involved in the mixture variable. If the model has no covariates specify
#' \code{y ~ NULL}.
#' @param mix_name A character vector listing the variables contributing to a mixture effect.
#' @param cluster_name A character string that specifiy which is the column of the dataset which
#' contains the group number. Note that the \code{cluster_name} should be numeric, strings and factors
#' are not allowed
#' @param data The \code{data.frame} containing the variables (covariates and elements of the mixture)
#' to be included in the model.
#' @param q An \code{integer} to specify how mixture variables will be ranked, e.g. in quartiles
#' (\code{q = 4}), deciles (\code{q = 10}), or percentiles (\code{q = 100}). If \code{q = NULL} then
#' the values of the mixture variables are taken (these must be standardized or the domain must be the same).
#' @param Dalp A \code{vector} containing the parameters of the Dirichlet distribution of the weights, the number
#' of the elements of the vector has to be equal to the number of chemicals. If \code{Dalp = NULL}
#' the domain is explored uniformly.
#' @param chains An \code{integer} to specify the number of chain in Hamiltonian Monte Carlo algorithm.
#' Default value \code{chains = 1}.
#' @param iter An \code{integer} to specify the lenght of chain in Hamiltonian Monte Carlo algorithm.
#' Default value \code{iter = 10000}.
#' @param thin An \code{integer} to specify the thinning parameter in Hamiltonian Monte Carlo algorithm.
#' @param seed An \code{integer} value to fix the seed. If \code{seed = NULL} the seed are randomly choosen.
#' @param start_value A \code{vector} containing the initial value of the prior distribution,
#' if it is equal to NULL random values are chosen.
#' @param c_int A \code{vector} of two elements to specify the credible intervals for parameters, for 95% credible
#' interval \code{c_int = c(0.025,0.975)} (default).
#' @param family A \code{string} to specify the type of outcome. With the current implementation  the
#' possible values are only continuous - "gaussian" (default).
#'
#' @details
#' The function \code{bwqs} uses the package \code{rstan} which allows the connection with STAN,
#' a specific software, written in C++ for bayesian inference, for further information see https://mc-stan.org/.
#'
#' @return
#' \code{bwqs} returns a list with two argument:
#' \item{fit}{An \code{S4} object with all details of the Hamiltonian Monte Carlo, all the extractions
#' from the posterior distribution and all values of the parameters}
#' \item{summary_fit}{Table with the statistics of the parameters: mean, standard error of the mean,
#' standard deviation, lower and upper values for the credible interval (with credible level specified
#' by \code{c_int}), n_eff and Rhat. For further details see https://cran.r-project.org/web/packages/rstan/rstan.pdf}
#'
#' @author
#' Nicolo Foppa Pedretti, Elena Colicino
#'
#' @import rstan
#' @import Rcpp
#' @import methods
#'
#' @export

bwqs_r <- function(formula, mix_name, cluster_name, data, q, Dalp = NULL,
                   chains = 1, iter = 1000, thin = 3, seed=2019, start_value=NULL,
                   c_int=c(0.025,0.975), family="gaussian"){

  formula = as.formula(formula)
  y_name  <- all.vars(formula)[1]
  KV_name <- all.vars(formula)[-1]
  X_name  <- mix_name

  check_input_r(formula, mix_name, cluster_name, data, q, Dalp,
                chains, iter, thin, seed,
                start_value, c_int, family)

  if(length(KV_name)==0){
    data = as.data.frame(data[,c(y_name,X_name,cluster_name)])
  } else{
    data = as.data.frame(data[,c(y_name,KV_name,X_name,cluster_name)])
  }

  data <- na.omit(data)
  if(nrow(data)==0) stop("No dataset available")

  if(is.null(q)){
    Chem = data[,mix_name]
  } else{
    Chem = quantile_split_2(data=data,mix_name=mix_name,q)[,mix_name]
  }

  if(length(KV_name)==0){
    KV = NULL
  } else{
    KV = data[,KV_name]
  }

  if(is.null(Dalp)){
    Dalp = rep(1, length(mix_name))
  } else{
    Dalp = Dalp
  }

  if(is.null(start_value)){
    start_value = "random"
  }

  if(!(cluster_name %in% colnames(data))){
    stop("Cluster specification not found")
  }

  switch (family,
          gaussian = {if(!is.null(KV)){
            data_reg <- list(
              N      = nrow(data),
              J      = length(unique(data[,cluster_name])),
              C      = length(mix_name),
              K      = length(KV_name),
              cohort = as.vector(data[,cluster_name]),
              Chem   = cbind(X),
              X      = cbind(KV),
              Dalp   = Dalp,
              y      = as.vector(data[,y_name])
            )

            fit <- stan(model_code = model_rbwqs_regression_cov,
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
              N      = nrow(data),
              J      = length(unique(data[,cluster_name])),
              C      = length(mix_name),
              cohort = as.vector(data[,cluster_name]),
              Chem   = cbind(X),
              Dalp   = Dalp,
              y      = as.vector(data[,y_name])
            )

            fit <- stan(model_code = model_rbwqs_regression,
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

  if(length(KV_name)!=0){
      parameters <- c(paste0('W_',mix_name),
                      "sigma",
                      "sigma_b0",
                      "sigma_b1",
                      paste0('beta0_',seq(1,length(unique(data[,cluster_name])))),
                      paste0('beta1_',seq(1,length(unique(data[,cluster_name])))),
                      paste0('Cov_',KV_name))
    }else{
      parameters <- c(paste0('W_',mix_name),
                      "sigma",
                      "sigma_b0",
                      "sigma_b1",
                      paste0('beta0_',seq(1,length(unique(data[,cluster_name])))),
                      paste0('beta1_',seq(1,length(unique(data[,cluster_name])))))
    }

  names(fit)[1:length(parameters)] <- parameters
  sum_fit <- round(summary(fit,pars = parameters,
                           probs = c_int)$summary,5)
  return(list(fit=fit, summary_fit = sum_fit))
}








