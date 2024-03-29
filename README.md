# Bayesian Weighted Quantile Sum

This package provides the code for running BWQS regressions. 

In order to run the model, install Rtools (see the version of Rtools for your specific R version at this link: https://cran.r-project.org/bin/windows/Rtools/) 
on your machine and run the following command on your console: 

```
devtools::install_github("ElenaColicino/bwqs", build_vignettes = TRUE)
library(BWQS)
```

The additional argument `build_vignettes = TRUE` allows the package to build the vignette with all the examples. In order to read the vignette in your browser type:

```
browseVignettes("BWQS")
```

## Resources

For additional information on the BWQS regression please read this paper:
https://pubmed.ncbi.nlm.nih.gov/32613152/.

Update: the function `bwqs_r` for the *Hierarchical Bayesian Bayesian Weighted Quantile Sum* has been added to the package.  
