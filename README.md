# fastconstrained
An R package for quickly fitting the constrained lasso under one linear constraint using successive subsets from the unconstrained lasso.

# Starting out
Before using, please install and attach the required RcppArmadillo functions from https://github.com/lamttran/ADMMRcppArma with the following code (the devtools package is required):
- Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS=TRUE)
- install_github("lamttran/ADMMRcppArma", ref = "main") 
- library(ADMMRcppArma)

Do the same with the intasymm package:
- install_github("lamttran/fastconstrained", subdir="pkg", ref = "main") 
- library(fastconstrained)
