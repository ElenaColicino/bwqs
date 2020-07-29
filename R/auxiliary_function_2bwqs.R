check_input_2 <- function(formula, mix_name_1, mix_name_2, data, q, DalpC1, DalpC2,
                          chains, iter, thin, seed,
                          start_value, c_int, family){

  form = ifelse(class(formula) != "formula", TRUE, FALSE)
  if(form) stop("formula must be of class formula")

  chars1 = ifelse(!is.character(mix_name_1), TRUE, FALSE)
  if(chars1) stop("mix_name_1 must be a character vector")

  chars2 = ifelse(!is.character(mix_name_2), TRUE, FALSE)
  if(chars2) stop("mix_name_2 must be a character vector")

  df = ifelse(!is.data.frame(data), TRUE, FALSE)
  if(df) stop("data must be a data.frame")

  quant = ifelse((!is.numeric(q) | length(q) != 1) & !is.null(q), TRUE, FALSE)
  if(quant) stop("q must be numeric of length 1, or a NULL")

  dalp1 = ifelse((sum(is.na(DalpC1))!=0 & length(DalpC1)>1), TRUE, FALSE)
  if(dalp1) stop("DalpC1 must be numeric vector of lenght(mix_names) , or NA")

  dalp2 = ifelse((sum(is.na(DalpC2))!=0 & length(DalpC2)>1), TRUE, FALSE)
  if(dalp2) stop("DalpC2 must be numeric vector of lenght(mix_names) , or NA")

  chns = ifelse((!is.numeric(chains) | length(chains)!= 1), TRUE, FALSE)
  if(chns) stop("chains must be numeric of length 1")

  it = ifelse((!is.numeric(iter) | length(iter)!= 1), TRUE, FALSE)
  if(it) stop("iter must be numeric of length 1")

  thn = ifelse((!is.numeric(thin) | length(thin)!= 1), TRUE, FALSE)
  if(thn) stop("thin must be numeric of length 1")

  sd = ifelse((!is.numeric(seed) | length(seed)!= 1), TRUE, FALSE)
  if(sd) stop("seed must be numeric of length 1")

  st_v = ifelse((!is.null(start_value) & length(start_value)!= length(c(all.vars(formula),
                                                                        mix_name_1, mix_name_2))), TRUE, FALSE)
  if(st_v) stop("start_value must be numeric vector of length model parameters")

  conf_int = ifelse((!is.null(c_int) & !is.vector(c_int)), TRUE, FALSE)
  if(conf_int) stop("c_int must be numeric vector")

  fam = ifelse(family!="gaussian" & family!="binomial" & family!="poisson", TRUE, FALSE)
  if(fam) stop("family not available")

}
