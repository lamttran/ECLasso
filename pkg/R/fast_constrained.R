#' Function for solving constrained lasso
#' 
#' This function solves the constrained lasso using our candidate subset approach
#' @param x The design matrix
#' @param y The response vector
#' @param censor If using a Cox model, the censoring vector
#' @param C The constraint that linearly combines the regression coefficients. Should be in a 1xp matrix
#' @param d The scalar value that C times the regrssion coefficients should equal
#' @param family The type of model. Note for now only logistic and Cox models are supported.
#' @param penalty.factor The penalization vector that multiplies the lasso penalty lambda e.g. elements = 1
#' are lasso penalized and elements = 0 are unpenalized
#' @param depth How many lambdas into the unconstrained sequence to solve for
#' @param tolerance Tolerance of ADMM, influences the primal and dual criteria for termination
#' @param trace If false, returns regression coefficients. If true, also returns number of ADMM iterations required
#' @import glmnet
#' @import survival
#' @import MASS
#' @import spatstat.utils
#' @keywords fast_constrained
#' @export

fast_constrained <- function(x, y, censor = NULL, C = NULL, d = NULL,
                             family = c("cox", "logistic"),
                             penalty.factor = 1, depth = 15, tolerance = 1e-7,
                             trace = F){
  
  n = nrow(x)
  p = ncol(x)
  if (length(penalty.factor) == 1) penalty.factor = rep(penalty.factor, p)
  if (is.null(d)) {d = 0}
  if (is.null(C)) {#C is not specified, set it to a row of 1.
    C = matrix(1, nrow = 1, ncol = p)
  }
  
  
  if(family == "logistic"){
    censor = NULL
    fit1 = glmnet(x = x, y = y, family = "binomial", intercept = F, standardize = F,
                  penalty.factor = penalty.factor)
  } else if(family == "cox"){
    fit1 = glmnet(x = x, y = Surv(y, censor), family = "cox", standardize = F, 
                  penalty.factor = penalty.factor)
  }
  
  lambda_seq = n * fit1$lambda
  
  par_keep = list()
  par_keep[[1]] = which(coef(fit1)[, 1]!= 0)
  for(i in 2:length(fit1$lambda)){
    par_keep[[i]] = sort(union(par_keep[[i - 1]],
                               which(coef(fit1)[, i]!= 0)))
  }
  
  numparam = lengths(par_keep)
  unique_index = c()
  
  for(i in 1:length(unique(numparam))){
    unique_index[i] = min(which(numparam == unique(numparam)[i]))
  }
  
  nlambda = length(unique_index)
  
  solution = c()
  time.iter = c() #tracking the time to fit for a lambda
  i.glmnet.iter = c() #for a given admm lambda, which glmnet lambda does it correspond to?
  
  i.max.glmnet = 0
  i.glmnet = 2
  eta_list = list(iterate_1 = list())
  
  for(i in 1:depth){
    stop = F
    print(i)
    ptm <- proc.time()
    lambda = lambda_seq[i] 
    
    if(family == "logistic"){
      index_i = par_keep[[unique_index[min(i.glmnet, nlambda)]]] - 1
    } else if(family == "cox"){
      index_i = par_keep[[unique_index[min(i.glmnet, nlambda)]]]
    }
    
    while(length(index_i) < ifelse(family == "logistic", 1, 2)){
      i.glmnet = i.glmnet + 1
      index_i = par_keep[[unique_index[min(i.glmnet, nlambda)]]]
    }
    
    xdata.new = as.matrix(x[, index_i], ncol = length(index_i))
    
    if(family == "logistic"){
      fit.admm1 = onestepADMM(x = xdata.new, y = y,
                              family = "logistic", tol = tolerance, 
                              maxit = 100, lambda = lambda,
                              intercept = F, C = matrix(C[index_i], nrow = 1), d = d,
                              penalty.factor = penalty.factor[index_i])
    } else if(family == "cox"){
      fit.admm1 = onestepADMM(x = xdata.new, y = y, censor = censor, 
                              family = "cox", tol = tolerance, 
                              maxit = 100, lambda = lambda,
                              intercept = F, C = matrix(C[index_i], nrow = 1), d = d,
                              penalty.factor = penalty.factor[index_i])
    }
    
    new_guess = rep(0, p)
    new_guess[index_i] = inner.solution(admm.fit = fit.admm1, family = family, x = xdata.new, 
                                        y = y, intercept = F, censor = censor, 
                                        lambda = lambda, C = C[index_i], 
                                        penalty.factor = penalty.factor[index_i], d = d)
    
    
    kkt.outer = KKTtestfct(par_guess = new_guess, family = family, lambda = lambda, 
                           x = x, y = y, C = C, censor = censor)
    
    while(kkt.outer == F){
      i.glmnet = i.glmnet + 1
      if(i.glmnet > nlambda){
        stop == TRUE
        break
      }
      if(family == "logistic"){
        index_i = par_keep[[unique_index[min(i.glmnet, nlambda)]]] - 1
      } else if(family == "cox"){
        index_i = par_keep[[unique_index[min(i.glmnet, nlambda)]]]
      }
      
      #xdata.new = x[, index_i]
      xdata.new = as.matrix(x[, index_i], ncol = length(index_i))
      
      fit.admm1 = onestepADMM(x = xdata.new, y = y, censor = censor, 
                              family = family, tol = tolerance, 
                              maxit = 200, lambda = lambda,
                              intercept = F, C = matrix(C[index_i], nrow = 1), d = d,
                              penalty.factor = penalty.factor[index_i])
      
      new_guess = rep(0, p)
      new_guess[index_i] = inner.solution(admm.fit = fit.admm1, family = family, x = xdata.new, 
                                          y = y, censor = censor, intercept = F,
                                          lambda = lambda, C = C[index_i], 
                                          penalty.factor = penalty.factor[index_i], d = d)
      
      kkt.outer = KKTtestfct(par_guess = new_guess, family = family, lambda = lambda, 
                             x = x, y = y.logistic, C = C, censor = censor)
    }
    if(stop == T){break}
    solution = cbind(solution, new_guess)
    time.iter = c(time.iter, (proc.time() - ptm)[3])
    i.glmnet.iter = c(i.glmnet.iter, i.glmnet)
  }
  
  if(trace == F){
    return(list(solution = solution))
  } else if(trace == T){
    return(list(solution = solution, time = time.iter, subsets = i.glmnet.iter))
  }
}