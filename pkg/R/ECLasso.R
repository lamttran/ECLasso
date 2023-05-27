#' Function for solving constrained lasso.
#'
#' This function solves the constrained lasso using our candidate subset approach.
#' It is a wrapper for three fitting methods: combined candidate subset, ADMM-only candidate subset
#' and naive whole-predictor. 
#' @param x The design matrix
#' @param y The response vector
#' @param censor If using a Cox model, the censoring vector
#' @param C The constraints that linearly combine the regression coefficients. 
#' @param d The vector of values that C times the regression coefficients should equal
#' @param family The type of model. Note for now only logistic and Cox models are supported.
#' @param penalty.factor The penalization vector that multiplies the lasso penalty lambda e.g. elements = 1
#' are lasso penalized and elements = 0 are un-penalized
#' @param depth How many lambdas into the unconstrained sequence to solve for
#' @param tolerance Tolerance of ADMM, influences the primal and dual criteria for termination
#' @param method One of "combined_subset", "admm_subset", or "naive". The first two
#' implement our candidate subset approach, the former using a combined ADMM/Newton-Raphson algorithm
#' and the latter using ADMM only. The last solves the constrained lasso using ADMM only on the whole predictor space.
#' @param trace If false, returns solution path. If true, also returns run-time, subsets, and penalty sequence
#' @import glmnet
#' @import survival
#' @import MASS
#' @import spatstat.utils
#' @import Matrix
#' @keywords ECLasso
#' @export

ECLasso <- function(x, y, censor = NULL, C = NULL, d = NULL,
                    family = c("cox", "logistic"), penalty.factor = 1, depth = 15, 
                    tolerance = 1e-7, method = c("combined_subset", "admm_subset", "naive"), 
                    intercept = F, trace = F){
  
  family = match.arg(family)

  if(method == "combined_subset"){
    constrained_fit = fast_constrained(x = x, y = y, censor = censor, C = C, d = d, intercept = intercept,
                                       family = family, penalty.factor = penalty.factor, 
                                       depth = depth, tolerance = tolerance, trace = trace)
  } else if(method == "admm_subset"){
    constrained_fit = fast_no_nr(x = x, y = y, censor = censor, C = C, d = d, intercept = intercept,
                                 family = family, penalty.factor = penalty.factor, 
                                 depth = depth, tolerance = tolerance, trace = trace)
  } else if(method == "naive"){
    constrained_fit = naive_fit(x = x, y = y, censor = censor, C = C, d = d, intercept = intercept,
                                family = family, penalty.factor = penalty.factor, 
                                depth = depth, tolerance = tolerance, trace = trace)
  }
  
  if(trace == F){
    return(list(solution = constrained_fit$solution))
  } else if(trace == T){
    return(list(solution = constrained_fit$solution, time = constrained_fit$time, 
                lambda_seq = constrained_fit$lambda_seq))
  }
}
