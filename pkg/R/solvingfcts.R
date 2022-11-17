#' Solvers of the constrained lasso, called by the main ECLasso function
#'
#' @param x The design matrix
#' @param y The response vector
#' @param censor If using a Cox model, the censoring vector
#' @param C The constraints that linearly combine the regression coefficients. 
#' @param d The vector of values that C times the regrssion coefficients should equal
#' @param family The type of model. Note for now only logistic and Cox models are supported.
#' @param penalty.factor The penalization vector that multiplies the lasso penalty lambda e.g. elements = 1
#' are lasso penalized and elements = 0 are unpenalized
#' @param depth How many lambdas into the unconstrained sequence to solve for
#' @param tolerance Tolerance of ADMM, influences the primal and dual criteria for termination
#' @param trace If false, returns solution path. If true, also returns run-time, subsets, and penalty sequence
#' @name solvingfcts
NULL
#> NULL

#' @rdname solvingfcts
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
  i.glmnet = 3
  eta_list = list(iterate_1 = list())
  kkt.failed.indices = NULL
  
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
    
    index_i = sort(union(index_i, kkt.failed.indices))
    
    xdata.new = as.matrix(x[, index_i], ncol = length(index_i))
    
    if(family == "logistic"){
      fit.admm1 = onestepADMM(x = xdata.new, y = y,
                              family = "logistic", tol = tolerance, 
                              maxit = 100, lambda = lambda,
                              intercept = F, C = matrix(C[,index_i], nrow = nrow(C)), d = d,
                              penalty.factor = penalty.factor[index_i])
    } else if(family == "cox"){
      fit.admm1 = onestepADMM(x = xdata.new, y = y, censor = censor, 
                              family = "cox", tol = tolerance, 
                              maxit = 100, lambda = lambda,
                              intercept = F, C = matrix(C[,index_i], nrow = nrow(C)), d = d,
                              penalty.factor = penalty.factor[index_i])
    }
    
    #if the inner solution is trivial, so is the whole solution
    if(fit.admm1$iteration[2] != 200 & max(abs(fit.admm1$z.out) == 0)){
      new_guess = rep(0, p)
      new_guess[index_i] = fit.admm1$z.out
      solution = cbind(solution, new_guess)
      time.iter = c(time.iter, (proc.time() - ptm)[3])
      i.glmnet.iter = c(i.glmnet.iter, i.glmnet)
      next 
    }
    
    #if the subsetted constraint matrix is not full rank, use the next unconstrained subset
    fullrank = (rankMatrix(C[,index_i, drop = F]) == nrow(C))
    if(fullrank == T){
      new_guess = rep(0, p)
      new_guess[index_i] = inner.solution(admm.fit = fit.admm1, family = family, x = xdata.new, 
                                          y = y, intercept = F, censor = censor, 
                                          lambda = lambda, C = C[,index_i, drop = F], 
                                          penalty.factor = penalty.factor[index_i], d = d)
      
      
      fullset.test = KKTtestfct(par_guess = new_guess, family = family, lambda = lambda, 
                                x = x, y = y, C = C, censor = censor, trace = T)
      kkt.outer = fullset.test$kkt.logical
      
      if(kkt.outer == F){
        kkt.failed.indices = sort(union(kkt.failed.indices, fullset.test$failed.indices))
      }
    } else{
      kkt.outer = F
    }
    
    
    
    while(kkt.outer == F){
      i.glmnet = i.glmnet + 1
      if(i.glmnet > nlambda){ #i.e. if none of the subsets work, fit naive ADMM
        stop == TRUE
        fit.admm1 = onestepADMM(x = x, y = y, censor = censor, 
                                family = family, tol = tolerance, 
                                maxit = 20000, lambda = lambda,
                                intercept = F, C = C, d = d,
                                penalty.factor = penalty.factor)
        new_guess = fit.admm1$z.out
        break
      }
      if(family == "logistic"){
        index_i = par_keep[[unique_index[min(i.glmnet, nlambda)]]] - 1
      } else if(family == "cox"){
        index_i = par_keep[[unique_index[min(i.glmnet, nlambda)]]]
      }
      
      index_i = sort(union(index_i, kkt.failed.indices))
      
      #xdata.new = x[, index_i]
      xdata.new = as.matrix(x[, index_i], ncol = length(index_i))
      
      fit.admm1 = onestepADMM(x = xdata.new, y = y, censor = censor, 
                              family = family, tol = tolerance, 
                              maxit = 200, lambda = lambda,
                              intercept = F, C = matrix(C[,index_i], nrow = nrow(C)), d = d,
                              penalty.factor = penalty.factor[index_i])
      
      fullrank = (rankMatrix(C[,index_i, drop = F]) == nrow(C))
      if(fullrank == T){
        new_guess = rep(0, p)
        new_guess[index_i] = inner.solution(admm.fit = fit.admm1, family = family, x = xdata.new, 
                                            y = y, intercept = F, censor = censor, 
                                            lambda = lambda, C = C[,index_i, drop = F], 
                                            penalty.factor = penalty.factor[index_i], d = d)
        
        
        fullset.test = KKTtestfct(par_guess = new_guess, family = family, lambda = lambda, 
                                  x = x, y = y, C = C, censor = censor, trace = T)
        kkt.outer = fullset.test$kkt.logical
        if(kkt.outer == F){
          kkt.failed.indices = sort(union(kkt.failed.indices, fullset.test$failed.indices))
        }
      } else{
        kkt.outer = F
      }
    }
    if(stop == T){break}
    solution = cbind(solution, new_guess)
    time.iter = c(time.iter, (proc.time() - ptm)[3])
    i.glmnet.iter = c(i.glmnet.iter, i.glmnet)
  }
  
  if(trace == F){
    return(list(solution = solution))
  } else if(trace == T){
    return(list(solution = solution, time = time.iter, subsets = i.glmnet.iter, lambda_seq = lambda_seq[1:depth]))
  }
}

#' @rdname solvingfcts
fast_no_nr <- function(x, y, censor = NULL, C = NULL, d = NULL,
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
  i.glmnet = 3
  eta_list = list(iterate_1 = list())
  kkt.failed.indices = NULL
  
  
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
    
    index_i = sort(union(index_i, kkt.failed.indices))
    
    xdata.new = as.matrix(x[, index_i], ncol = length(index_i))
    
    if(family == "logistic"){
      fit.admm1 = onestepADMM(x = xdata.new, y = y,
                              family = "logistic", tol = tolerance, 
                              maxit = 50000, lambda = lambda,
                              intercept = F, C = matrix(C[,index_i], nrow = nrow(C)), d = d,
                              penalty.factor = penalty.factor[index_i])
    } else if(family == "cox"){
      fit.admm1 = onestepADMM(x = xdata.new, y = y, censor = censor, 
                              family = "cox", tol = tolerance, 
                              maxit = 50000, lambda = lambda,
                              intercept = F, C = matrix(C[,index_i], nrow = nrow(C)), d = d,
                              penalty.factor = penalty.factor[index_i])
    }
    
    #if the subsetted constraint matrix is not full rank, use the next unconstrained subset
    fullrank = (rankMatrix(C[,index_i, drop = F]) == nrow(C))
    if(fullrank == T){
      new_guess = rep(0, p)
      new_guess[index_i] = fit.admm1$z.out
      
      fullset.test = KKTtestfct(par_guess = new_guess, family = family, lambda = lambda, 
                                x = x, y = y, C = C, censor = censor, trace = T)
      kkt.outer = fullset.test$kkt.logical
      
      if(kkt.outer == F){
        kkt.failed.indices = sort(union(kkt.failed.indices, fullset.test$failed.indices))
      }
    } else{
      kkt.outer = F
    }
    
    while(kkt.outer == F){
      i.glmnet = i.glmnet + 1
      if(i.glmnet > nlambda){
        stop == TRUE
        fit.admm1 = onestepADMM(x = x, y = y, censor = censor, 
                                family = family, tol = tolerance, 
                                maxit = 50000, lambda = lambda,
                                intercept = F, C = C, d = d,
                                penalty.factor = penalty.factor)
        new_guess = fit.admm1$z.out
        break
      }
      if(family == "logistic"){
        index_i = par_keep[[unique_index[min(i.glmnet, nlambda)]]] - 1
      } else if(family == "cox"){
        index_i = par_keep[[unique_index[min(i.glmnet, nlambda)]]]
      }
      
      index_i = sort(union(index_i, kkt.failed.indices))
      
      #xdata.new = x[, index_i]
      xdata.new = as.matrix(x[, index_i], ncol = length(index_i))
      
      fullrank = (rankMatrix(C[,index_i, drop = F]) == nrow(C))
      if(fullrank == T){
        new_guess = rep(0, p)
        new_guess[index_i] = onestepADMM(x = xdata.new, y = y, censor = censor, 
                                         family = family, tol = tolerance, 
                                         maxit = 50000, lambda = lambda,
                                         intercept = F, C = matrix(C[,index_i], nrow = nrow(C)), d = d,
                                         penalty.factor = penalty.factor[index_i])$z.out
        
        
        fullset.test = KKTtestfct(par_guess = new_guess, family = family, lambda = lambda, 
                                  x = x, y = y, C = C, censor = censor, trace = T)
        kkt.outer = fullset.test$kkt.logical
        if(kkt.outer == F){
          kkt.failed.indices = sort(union(kkt.failed.indices, fullset.test$failed.indices))
        }
      } else{
        kkt.outer = F
      }
    }
    if(stop == T){break}
    solution = cbind(solution, new_guess)
    time.iter = c(time.iter, (proc.time() - ptm)[3])
    i.glmnet.iter = c(i.glmnet.iter, i.glmnet)
  }
  
  if(trace == F){
    return(list(solution = solution))
  } else if(trace == T){
    return(list(solution = solution, time = time.iter, subsets = i.glmnet.iter, lambda_seq = lambda_seq[1:depth]))
  }
}


#' @rdname solvingfcts
naive_fit <- function(x, y, censor = NULL, C = NULL, d = NULL,
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
  
  solution = c()
  time.iter = c() #tracking the time to fit for a lambda
  
  i.max.glmnet = 0
  i.glmnet = 3
  eta_list = list(iterate_1 = list())
  
  for(i in 1:depth){
    stop = F
    print(i)
    ptm <- proc.time()
    lambda = lambda_seq[i] 
    
    new_guess = onestepADMM(x = x, y = y, censor = censor, 
                            family = family, tol = tolerance, 
                            maxit = 50000, lambda = lambda,
                            intercept = F, C = C, d = d,
                            penalty.factor = penalty.factor)
    
    solution = cbind(solution, new_guess$z.out)
    time.iter = c(time.iter, (proc.time() - ptm)[3])
  }
  
  if(trace == F){
    return(list(solution = solution))
  } else if(trace == T){
    return(list(solution = solution, time = time.iter, lambda_seq = lambda_seq[1:depth]))
  }
} 
