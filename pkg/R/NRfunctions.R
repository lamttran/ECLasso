#' Unconstrained/constrained Newton-Raphson, with gradient and Hessian functions
#' 
#' @param par_vec numeric vector of regression coefficients
#' @param xdata design matrix
#' @param lambda lasso penalty
#' @param censor if using a Cox model, the censoring vector
#' @param ydata the response vector
#' @param penalty.factor The penalization vector that multiplies the lasso penalty lambda e.g. elements = 1
#' are lasso penalized and elements = 0 are unpenalized
#' @param tolerance If successive N-R iterations differ by less than this, terminate
#' @param family the type of model, one of Cox or logistic
#' @name NRfunctions
NULL
#> NULL

#' @rdname NRfunctions
grad_sign_fct = function(par_vec, xdata, lambda, censor = NULL, 
                         ydata = NULL, penalty.factor = 1, 
                         family = c("cox", "logistic")){
  
  if(length(par_vec) == 1){
    xdata = as.matrix(xdata, ncol = 1)
  }
  
  family = match.arg(family)
  if (length(penalty.factor) == 1) penalty.factor = rep(penalty.factor, length(par_vec))
  
  if(family == "cox"){
    theta = exp(xdata %*% par_vec)
    x.theta = xdata * c(theta)
    thetafrac = apply(x.theta, 2, revcumsum)/c(cumsum_rev(theta))
    grad_vec = colSums(censor * (-xdata + thetafrac))
  } else if(family == "logistic"){
    pi = 1/(1 + exp(-xdata %*% par_vec))
    grad_vec = t(xdata) %*% (pi - ydata) #right one
  }
  
  grad_vec = grad_vec + (sign(par_vec) * penalty.factor * lambda)
  return(grad_vec)
}

#' @rdname NRfunctions
hessian_fct = function(par_vec, xdata, censor = NULL,
                       family = c("cox", "logistic")){
  family = match.arg(family)
  
  if(length(par_vec) == 1){
    xdata = as.matrix(xdata, ncol = 1)
  }
  
  if(family == "cox"){
    theta = exp(xdata %*% par_vec)
    hazard = as.numeric(theta)
    risk = c(cumsum_rev(hazard))
    P = outer(hazard, risk, '/')
    P[upper.tri(P)] = 0
    #for log-likelihood
    W = -P %*% diag(censor) %*% t(P)
    diag(W) = diag(P %*% diag(censor) %*% t(1 - P))
    hessian = t(xdata) %*% W %*% xdata
  } else if(family == "logistic"){
    pi = 1/(1 + exp(-xdata %*% par_vec))
    hessian = t(xdata) %*% diag(c(pi * (1 - pi))) %*% xdata
  }
  
  return(hessian)
}

#' @rdname NRfunctions
NR_constrained = function(par_vec, xdata, lambda, 
                          censor = NULL, ydata = NULL, C = NULL,
                          penalty.factor = 1, tolerance = 1e-6,
                          family = c("cox", "logistic")){
  
  family = match.arg(family)
  if (is.null(C)) {
    #C is not specified, set it to a row of 1.
    C = matrix(1, nrow = 1, ncol = nvars)
  }
  
  if (length(penalty.factor) == 1) penalty.factor = rep(penalty.factor, length(par_vec))
  
  ##iteration counter
  iter = 1
  grad = grad_sign_fct(par_vec, xdata, lambda = lambda, 
                       censor = censor, ydata = ydata, 
                       penalty.factor = penalty.factor, family = family)
  hess = hessian_fct(par_vec, xdata, censor = censor, family = family)
  
  
  hess.augmented = rbind(hess, C)
  hess.augmented = cbind(hess.augmented, c(C, 0))
  
  solution = par_vec + (solve(hess.augmented) %*% c(-grad, 0))[-(ncol(hess) + 1)] 
  
  solution_matrix = cbind(par_vec, solution)
  
  while(norm(c(solution_matrix[,(iter + 1)] - solution_matrix[,iter]), type = "2") > tolerance){
    iter = iter + 1
    par_vec = solution_matrix[,iter]
    
    grad = grad_sign_fct(par_vec, xdata, lambda = lambda, 
                         censor = censor, ydata = ydata, 
                         penalty.factor = penalty.factor, family = family)
    hess = hessian_fct(par_vec, xdata, censor = censor, family = family)
    
    hess.augmented = rbind(hess, C)
    hess.augmented = cbind(hess.augmented, c(C, 0))
    
    solution = par_vec + (solve(hess.augmented) %*% c(-grad, 0))[-(ncol(hess) + 1)] 
    solution_matrix = cbind(solution_matrix, solution)
    
    if(iter == 10){
      break
    }
  }
  
  return(list(solution = solution, solution_matrix = solution_matrix, iter = iter))
}

#' @rdname NRfunctions
NR_unconstrained = function(par_vec, xdata, lambda, 
                            censor = NULL, ydata = NULL,
                            penalty.factor = 1, tolerance = 1e-6,
                            family = c("cox", "logistic")){
  family = match.arg(family)
  
  if (length(penalty.factor) == 1) penalty.factor = rep(penalty.factor, length(par_vec))
  
  ##iteration counter
  iter = 1
  grad = grad_sign_fct(par_vec, xdata, lambda = lambda, 
                       censor = censor, ydata = ydata, 
                       penalty.factor = penalty.factor, family = family)
  hess = hessian_fct(par_vec, xdata, censor = censor, family = family)
  
  solution = par_vec + (solve(hess) %*% -grad)
  
  solution_matrix = cbind(par_vec, solution)
  
  while(norm(c(solution_matrix[,(iter + 1)] - solution_matrix[,iter]), type = "2") > tolerance){
    iter = iter + 1
    par_vec = solution_matrix[,iter]
    
    grad = grad_sign_fct(par_vec, xdata, lambda = lambda, 
                         censor = censor, ydata = ydata, 
                         penalty.factor = penalty.factor, family = family)
    hess = hessian_fct(par_vec, xdata, censor = censor, family = family)
    
    solution = par_vec + (solve(hess) %*% -grad)
    
    solution_matrix = cbind(solution_matrix, solution)
    
  }
  return(list(solution = solution, solution_matrix = solution_matrix, iter = iter))
  
}