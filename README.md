# ECLasso
An R package for quickly fitting the constrained lasso under one linear constraint using successive subsets from the unconstrained lasso.

# Starting out
Before using, please install and attach the required RcppArmadillo functions from https://github.com/lamttran/ADMMRcppArma with the following code (the devtools package is required):
- Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS=TRUE)
- install_github("lamttran/ADMMRcppArma", ref = "main") 
- library(ADMMRcppArma)

Do the same with the intasymm package:
- install_github("lamttran/ECLasso", subdir="pkg", ref = "main") 
- library(ECLasso)

# Workflow
The main function in the package is fast_constrained, which takes your data (x, y, and potentially a censoring vector) and the following:
- The model type, for now one of logistic or Cox models.
- The constraint vector (C) which linearly combines the model regression coefficients to equal the scalar d.
- A penalty.factor, either a scalar 1 or a vector equal to the number of parameters. This multiples the lasso penalty such that elements equal to 1 are lasso-penalized and elements equal to 0 are unpenalized.
- Depth, which is how far along the unconstrained penalty sequence to solve for. Generally anywhere from 10-20 is sufficient, as there may not be enough subsets to solve for deeper depth.
- Tolerance, which determines termination criteria for the internal ADMM algorithm.

# Try it out with simulated data
You can perform a simple sum-to-zero constrained lasso example with logistic and Cox models using simulated data with the script ADMMsimulations.R
