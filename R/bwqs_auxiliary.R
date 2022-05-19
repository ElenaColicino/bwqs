check_input <- function(formula, mix_name, data, q, Dalp,
                        chains, iter, thin, seed,
                        start_value, c_int, family, prior){

  form = ifelse(class(formula) != "formula", TRUE, FALSE)
  if(form) stop("formula must be of class formula")

  chars = ifelse(!is.character(mix_name), TRUE, FALSE)
  if(chars) stop("mix_name must be a character vector")

  df = ifelse(!is.data.frame(data), TRUE, FALSE)
  if(df) stop("data must be a data.frame")

  quant = ifelse((!is.numeric(q) | length(q) != 1) & !is.null(q), TRUE, FALSE)
  if(quant) stop("q must be numeric of length 1, or a NULL")

  dalp = ifelse((sum(is.na(Dalp))!=0 & length(Dalp)>1), TRUE, FALSE)
  if(dalp) stop("Dalp must be numeric vector of lenght(mix_names) , or NA")

  chns = ifelse((!is.numeric(chains) | length(chains)!= 1), TRUE, FALSE)
  if(chns) stop("chains must be numeric of length 1")

  it = ifelse((!is.numeric(iter) | length(iter)!= 1), TRUE, FALSE)
  if(it) stop("iter must be numeric of length 1")

  thn = ifelse((!is.numeric(thin) | length(thin)!= 1), TRUE, FALSE)
  if(thn) stop("thin must be numeric of length 1")

  sd = ifelse((!is.numeric(seed) | length(seed)!= 1), TRUE, FALSE)
  if(sd) stop("seed must be numeric of length 1")

  st_v = ifelse((!is.null(start_value) & length(start_value)!= length(c(all.vars(formula),mix_name))), TRUE, FALSE)
  if(st_v) stop("start_value must be numeric vector of length model parameters")

  conf_int = ifelse((!is.null(c_int) & !is.vector(c_int)), TRUE, FALSE)
  if(conf_int) stop("c_int must be numeric vector")

  fam = ifelse(family!="gaussian" & family!="binomial" & family!="poisson", TRUE, FALSE)
  if(fam) stop("family not available")
  
  pr = ifelse(prior!="None" & prior!="positive" & prior!="negative", TRUE, FALSE)
  if(pr) stop("prior not available")

}


colVars <- function(a){
  n <- dim(a)[[1]]
  c <- dim(a)[[2]]
  return(.colMeans(((a - matrix(.colMeans(a, n, c), nrow = n, ncol = c, byrow = TRUE)) ^ 2), n, c) * n / (n - 1))
}

check_input_r <- function(formula, mix_name, cluster_name, data, q, Dalp,
                          chains, iter, thin, seed, 
                          start_value, c_int, family, rw, SW){
  
  form = ifelse(class(formula) != "formula", TRUE, FALSE)
  if(form) stop("formula must be of class formula")
  
  chars1 = ifelse(!is.character(mix_name), TRUE, FALSE)
  if(chars1) stop("mix_name must be a character vector")
  
  df = ifelse(!is.data.frame(data), TRUE, FALSE)
  if(df) stop("data must be a data.frame")
  
  quant = ifelse((!is.numeric(q) | length(q) != 1) & !is.null(q), TRUE, FALSE)
  if(quant) stop("q must be numeric of length 1, or a NULL")
  
  dalp = ifelse((sum(is.na(Dalp))!=0 & length(Dalp)>1), TRUE, FALSE)
  if(dalp) stop("Dalp must be numeric vector of lenght(mix_names) , or NA")
  
  clus = ifelse((!is.character(cluster_name)), TRUE, FALSE)
  if(clus) stop("cluster_name must be character vector")
  
  chns = ifelse((!is.numeric(chains) | length(chains)!= 1), TRUE, FALSE)
  if(chns) stop("chains must be numeric of length 1")
  
  it = ifelse((!is.numeric(iter) | length(iter)!= 1), TRUE, FALSE)
  if(it) stop("iter must be numeric of length 1")
  
  thn = ifelse((!is.numeric(thin) | length(thin)!= 1), TRUE, FALSE)
  if(thn) stop("thin must be numeric of length 1")
  
  sd = ifelse((!is.numeric(seed) | length(seed)!= 1), TRUE, FALSE)
  if(sd) stop("seed must be numeric of length 1")
  
  st_v = ifelse((!is.null(start_value) & length(start_value)!= length(c(all.vars(formula),mix_name))), TRUE, FALSE)
  if(st_v) stop("start_value must be numeric vector of length model parameters")
  
  conf_int = ifelse((!is.null(c_int) & !is.vector(c_int)), TRUE, FALSE)
  if(conf_int) stop("c_int must be numeric vector")
  
  fam = ifelse(family!="gaussian", TRUE, FALSE)
  if(fam) stop("family not available")
}
