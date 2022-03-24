#' Check if KKT conditions hold
#'
#' This function checks if the proposed solution satisfies the KKT conditions for optimality
#' @param par_guess The guess of the regression coefficients
#' @param family The type of model. For now, only logistic and Cox models are supported
#' @param lambda The value of the lasso penalty
#' @param x The design matrix
#' @param y The response vector
#' @param censor If using a Cox model, the censoring vector
#' @param C The linear constraint, as a 1 by p matrix
#' @param trace If false, returns true/false output for KKT conditions. If true, returns that and the following:
#' whether the KKT conditions hold for each of unconstrained and constrained terms, and the eta vectors.
#' Note for KKT to hold the max of the eta_lower vector must be leq the min of the eta_upper vector
#' @keywords KKTtest
#' @export

KKTtestfct = function(par_guess, family = c("logistic", "poisson", "cox"), lambda,
                      x, y = NULL, censor = NULL, C = matrix(1, nrow = 1, ncol = ncol(X)),
                      trace = F){
  
  
  if(family == "logistic"){
    pi = 1/(1 + exp(-x %*% par_guess))
    gradvec = t(x) %*% (pi - y)
  } else if(family == "cox"){
    theta = exp(x %*% par_guess)
    x.theta = x * c(theta)
    thetafrac = apply(x.theta, 2, revcumsum)/c(cumsum_rev(theta))
    gradvec = colSums(censor * (-x + thetafrac))
  }
  
  #Conditions for constrained and unconstrained variables
  constrained.index = which(C != 0)
  
  if(length(constrained.index) > 0){
    sup_eta_lower = round(max(-gradvec[constrained.index] - lambda), 5)
    inf_eta_upper = round(min(-gradvec[constrained.index] + lambda), 5)
    
    constrained.logical = (sup_eta_lower <= inf_eta_upper)
  } else{
    constrained.logical = T
  }
  
  unconstrained.index = which(C == 0)
  if(length(unconstrained.index) > 0){
    unconstrained.grad = round(gradvec[unconstrained.index], 4)
    unconstrained.logical = (max(abs(unconstrained.grad)) == 0)
  } else{
    unconstrained.logical = T
  }
  
  if(constrained.logical == T & unconstrained.logical == T){
    kkt.logical = T
  } else{
    kkt.logical = F
  }
  
  if(trace == T){
    return(list(kkt.logical = kkt.logical, 
                eta_lower = -gradvec[constrained.index] - lambda, 
                eta_upper = -gradvec[constrained.index] + lambda,
                unconstrained.logical = unconstrained.logical, 
                constrained.logical = constrained.logical,
                gradvec = gradvec))
  } else{
    return(kkt.logical)
  }
}