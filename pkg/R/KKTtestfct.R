#' Check if KKT conditions hold
#'
#' This function checks if the proposed solution satisfies the KKT conditions for optimality
#' @param par_guess The guess of the regression coefficients
#' @param family The type of model. For now, only logistic and Cox models are supported
#' @param lambda The value of the lasso penalty
#' @param x The design matrix
#' @param y The response vector
#' @param censor If using a Cox model, the censoring vector
#' @param C The linear constraint(s), by default a 1 by p matrix
#' @param trace If false, returns true/false output for KKT conditions. If true, returns that and
#' which predictors failed the KKT conditions, which are added to later candidate subsets.
#' @keywords KKTtest
#' @export

KKTtestfct = function(par_guess, family = c("logistic", "cox"), lambda, x, y = NULL, 
                      censor = NULL, C = matrix(1, nrow = 1, ncol = ncol(x)), trace = F){
  if(trace == T){
    failed.indices = NULL
  }
  
  
  if(family == "logistic"){
    pi = 1/(1 + exp(-x %*% par_guess))
    gradvec = t(x) %*% (pi - y)
  } else if(family == "cox"){
    theta = exp(x %*% par_guess)
    x.theta = x * c(theta)
    thetafrac = apply(x.theta, 2, revcumsum)/c(cumsum_rev(theta))
    gradvec = colSums(censor * (-x + thetafrac))
  }
  
  if(class(C)[1] == "numeric"){
    C = matrix(C, nrow = 1)
  }
  
  #Conditions for constrained and unconstrained variables
  constrained.index = which(colMeans(C) != 0)
  
  if(length(constrained.index) > 0){
    if(nrow(C) == 1){
      if(trace == F){
        sup_eta_lower = round(max((-gradvec[constrained.index] - lambda)/c(C[, constrained.index])), 1)
        inf_eta_upper = round(min((-gradvec[constrained.index] + lambda)/c(C[, constrained.index])), 1)
        constrained.logical = (sup_eta_lower <= inf_eta_upper)
      } else{
        eta_lower = round((-gradvec[constrained.index] - lambda)/c(C[, constrained.index]), 1)
        eta_upper = round((-gradvec[constrained.index] + lambda)/c(C[, constrained.index]), 1)
        constrained.logical = (max(eta_lower) <= min(eta_upper))
        failed.indices = which(eta_lower > min(eta_upper) | eta_upper < max(eta_lower))
      }
    } else{
      no_intersect = F
      nonzero.index = which(par_guess != 0)
      unique.constraints = unique(C, MARGIN = 2)
      kkt.indices = c()
      for(i in 1:ncol(unique.constraints)){
        intersection = intersect(which(colSums(C == c(unique.constraints[, i])) == nrow(C)), 
                                 nonzero.index)
        if(length(intersection) == 0 & i < ncol(unique.constraints)){
          # no_intersect = T
          # break
          next
        } else if(length(intersection) == 0 & i == ncol(unique.constraints)){
          break
        }
        kkt.indices = c(kkt.indices, intersection[1])
      }
      
      if(length(kkt.indices) > nrow(C)){
        kkt.indices = kkt.indices[1:nrow(C)]
      } else if(length(kkt.indices) < nrow(C)){
        no_intersect = T
      }
      
      if(no_intersect == T){
        constrained.logical = F
      } else if(no_intersect == F){
        kkt.matrix = t(C)[kkt.indices, ]
        kkt.y = -gradvec[kkt.indices] - lambda * sign(par_guess[kkt.indices])
        eta.vec = solve(kkt.matrix, kkt.y)
        s.vec = round((-gradvec - c(matrix(eta.vec, nrow = 1) %*% C))/lambda, 1)
        constrained.logical = ifelse(min(s.vec) >= -1 & max(s.vec) <= 1, T, F)
        
        if(trace == T){
          failed.indices = which(s.vec < -1 | s.vec > 1)
        }
      }
    }
  } else{
    constrained.logical = T
  }
  
  unconstrained.index = which(colMeans(C) == 0)
  if(length(unconstrained.index) > 0){
    unconstrained.grad = round(gradvec[unconstrained.index], 1)
    unconstrained.logical = (max(abs(unconstrained.grad)) == 0)
  } else{
    unconstrained.logical = T
  }
  
  #If both the constrained and unconstrained KKT conditions are true, output true
  if(constrained.logical == T & unconstrained.logical == T){
    kkt.logical = T
  } else{
    kkt.logical = F
  }
  
  if(trace == F){
    return(kkt.logical)
  } else{
    return(list(kkt.logical = kkt.logical, failed.indices = failed.indices, eta_lower = eta_lower,
                eta_upper = eta_upper))
  }
}
