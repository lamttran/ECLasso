# ECLasso
An R package for quickly fitting the constrained lasso under one or more linear constraints using candidate subsets from the unconstrained lasso. This package is an implentation of the approach described in the manuscript "A fast solution to the lasso problem with equality constraints."

# Starting out
Before using, please install and attach the required RcppArmadillo functions from https://github.com/lamttran/ADMMRcppArma with the following code (the devtools package is required):
- Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS=TRUE)
- install_github("lamttran/ADMMRcppArma", ref = "main") 
- library(ADMMRcppArma)

Do the same with the intasymm package:
- install_github("lamttran/ECLasso", subdir="pkg", ref = "main") 
- library(ECLasso)

# Workflow
The main function in the package is ECLasso, which takes your data (x, y, and potentially a censoring vector) and the following:
- The model type, for now one of logistic or Cox models.
- The constraint matrix (C) which linearly combines the model regression coefficients to equal the vector d.
- A penalty.factor, either a scalar 1 or a vector equal to the number of parameters. This multiples the lasso penalty such that elements equal to 1 are lasso-penalized and elements equal to 0 are unpenalized.
- Depth, which is how far along the unconstrained penalty sequence to solve for. Generally anywhere from 10-20 is sufficient, as there may not be enough subsets to solve for deeper depth.
- Tolerance, which determines termination criteria for the internal ADMM algorithm.
- Method, one of "combined_subset", "admm_subset", "naive". The first two implement our candidate subset approach, the former using a very efficient combined ADMM/Newton-Raphson algorithm and the latter using ADMM only. The last solves the constrained lasso using ADMM only on the whole predictor space.

# Try it out some simulated data
A script to quickly get our method running on simulated data is included in the repository, entited "newuser_simulations.R". A more involved script to replicate the simulations in the paper is "paper_simulations.R".
