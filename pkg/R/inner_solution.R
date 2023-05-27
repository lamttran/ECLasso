#' Function for constrained solution of subset (inner solution)
#'
#' This function uses repeated ADMM and Newton-Raphson to find a solution for the predictor subset
#' @param admm.fit The ADMM fit outputted by onestepADMM
#' @param family The type of model. Note for now only logistic and Cox models are supported.
#' @param x The design matrix
#' @param y The response vector
#' @param lambda The lasso penalty
#' @param censor If using a Cox model, the censoring vector
#' @param tolerance Tolerance of ADMM, influences the primal and dual criteria for termination
#' @param intercept Whether an unpenalized an unconstrained intercept should be added to the model
#' @param C The constraint that linearly combines the regression coefficients. Should be a matrix
#' @param penalty.factor The penalization vector that multiplies the lasso penalty lambda e.g. elements = 1
#' are lasso penalized and elements = 0 are unpenalized
#' @param d The scalar value that C times the regrssion coefficients should equal
#' @param trace If false, returns regression coefficients. If true, also returns number of ADMM iterations required
#' @keywords inner_solution
#' @export

inner.solution = function(admm.fit, family = c("logistic", "cox"), x, y = NULL, lambda,
                          censor = NULL, tolerance = 1e-7, intercept = F, C, penalty.factor, d = 0,
                          trace = F){
  family = match.arg(family)
  kkt.logical = F
  
  if(class(C)[1] == "numeric"){
    C = matrix(C, nrow = 1)
  }
  
  #which params are constrained and which are unconstrained
  constrained.index = which(colSums(C) != 0)
  unconstrained.index = which(colSums(C) == 0)
  
  num.iter = 100
  while(kkt.logical == F){
    admm.fit = onestepADMM(x = x, y = y, censor = censor, 
                           family = family, tol = tolerance, 
                           maxit = num.iter, lambda = lambda,
                           intercept = intercept, C = matrix(C, nrow = nrow(C)), d = d,
                           penalty.factor = penalty.factor, 
                           initial.guess = admm.fit$solution[, 2],
                           delta_0 = admm.fit$delta.out, u1_0 = admm.fit$u1.out,
                           u2_0 = admm.fit$u2.out, z_0 = admm.fit$z.out, hessian = admm.fit$hessian)
    
    #to pass to N-R, must first make the solution on constrained predictors feasible
    #but only in the predictors that were already nonzero
    
    par.nr = admm.fit$z.out
    nonzero.index = which(par.nr != 0)
    stop.zero.soln = F
    stop.one.soln = F
    stop.rank.soln = F
    stop.nr.soln = F
    while(length(nonzero.index) == 0 & admm.fit$iteration[2] == 2 * num.iter){
      num.iter = num.iter * 2
      admm.fit = onestepADMM(x = x, y = y, censor = censor, 
                             family = family, tol = tolerance, 
                             maxit = num.iter, lambda = lambda,
                             intercept = intercept, C = matrix(C, nrow = nrow(C)), d = d,
                             penalty.factor = penalty.factor, 
                             initial.guess = admm.fit$solution[,2],
                             delta_0 = admm.fit$delta.out, u1_0 = admm.fit$u1.out,
                             u2_0 = admm.fit$u2.out, z_0 = admm.fit$z.out, hessian = admm.fit$hessian)
      
      #to pass to N-R, must first make the solution on constrained predictors feasible
      #but only in the predictors that were already nonzero
      
      par.nr = admm.fit$z.out
      nonzero.index = which(par.nr != 0)
      
      if(admm.fit$iteration[2] != 2 * num.iter){
        stop.zero.soln = T
        break
      }
    }
    
    if(stop.zero.soln == T){break}
    
    constrained.nonzero = constrained.index[which(constrained.index %in% nonzero.index)]
    
    ##we can't make a guess for 1 non-zero predictor in sum-to-zero constraint
    if(length(d) == 1){
      if(length(constrained.nonzero) == 1 & d == 0){
        while(length(constrained.nonzero) == 1){
          num.iter = num.iter * 2
          admm.fit = onestepADMM(x = x, y = y, censor = censor, 
                                 family = family, tol = tolerance, 
                                 maxit = num.iter, lambda = lambda,
                                 intercept = intercept, C = matrix(C, nrow = nrow(C)), d = d,
                                 penalty.factor = penalty.factor, 
                                 initial.guess = admm.fit$solution[,2],
                                 delta_0 = admm.fit$delta.out, u1_0 = admm.fit$u1.out,
                                 u2_0 = admm.fit$u2.out, z_0 = admm.fit$z.out, hessian = admm.fit$hessian)
          
          if(admm.fit$iteration[2] != 2 * num.iter){
            stop.one.soln = T
            break
          }
          
          #to pass to N-R, must first make the solution on constrained predictors feasible
          #but only in the predictors that were already nonzero
          
          par.nr = admm.fit$z.out
          nonzero.index = which(par.nr != 0)
          constrained.nonzero = constrained.index[which(constrained.index %in% nonzero.index)]
        }
      }
    }
    if(stop.one.soln == T){break}
    
    
    
    if(length(constrained.nonzero) > 0){
      C_nonzero = C[,constrained.nonzero, drop = F]
      rank_condition = (rankMatrix(C_nonzero) == nrow(C_nonzero))
      
      #To find a feasible solution for N-R, we need to project into a potentially affine space
      #This requires a matrix inversion, and we need a constraint matrix with full rank
      if(rank_condition == F){
        while(rank_condition == F){
          num.iter = num.iter * 2
          admm.fit = onestepADMM(x = x, y = y, censor = censor, 
                                 family = family, tol = tolerance, 
                                 maxit = num.iter, lambda = lambda,
                                 intercept = intercept, C = matrix(C, nrow = nrow(C)), d = d,
                                 penalty.factor = penalty.factor, 
                                 initial.guess = admm.fit$solution[,2],
                                 delta_0 = admm.fit$delta.out, u1_0 = admm.fit$u1.out,
                                 u2_0 = admm.fit$u2.out, z_0 = admm.fit$z.out, hessian = admm.fit$hessian)
          if(admm.fit$iteration[2] != 2 * num.iter){
            stop.rank.soln = T
            break
          }
          par.nr = admm.fit$z.out
          nonzero.index = which(par.nr != 0)
          constrained.nonzero = constrained.index[which(constrained.index %in% nonzero.index)]
          C_nonzero = C[, constrained.nonzero, drop = F]
          rank_condition = (rankMatrix(C_nonzero) == nrow(C_nonzero))
        }
      }
      
      
      par_nonzero = par.nr[constrained.nonzero]
      par_nonzero = (diag(nrow = length(constrained.nonzero)) - t(C_nonzero) %*% solve(C_nonzero %*% t(C_nonzero)) %*%
                       C_nonzero) %*% par_nonzero + t(C_nonzero) %*% 
        solve(C_nonzero %*% t(C_nonzero)) %*% d
      par.nr[constrained.nonzero] = par_nonzero
    }
    if(stop.rank.soln == T){break}
    
    
    par.iter = 1
    #if any constrained predictors are nonzero, fit constrained N-R
    #if not, fit unconstrained N-R
    if(any(nonzero.index %in% constrained.index)){
      par.nr.out = try(NR_constrained(par.nr[nonzero.index], xdata = x[,nonzero.index], 
                                      censor = censor, lambda = lambda, 
                                      ydata = y, family = family,
                                      penalty.factor = penalty.factor[nonzero.index],
                                      C = C[, nonzero.index, drop = F], tol = tolerance), silent = T)
      
      if(!("try-error" %in% class(par.nr.out)) == T){
        par.iter = par.nr.out$iter
        par.nr.inner = par.nr.out$solution
      }
    } else{
      par.nr.out = try(NR_unconstrained(par.nr[nonzero.index], xdata = x[,nonzero.index], 
                                        censor = censor, lambda = lambda, 
                                        ydata = y, family = family,
                                        penalty.factor = penalty.factor[nonzero.index], tol = 1e-7))
      par.iter = par.nr.out$iter
      par.nr.inner = par.nr.out$solution
    }
    
    
    #if N-R fails to converge, we need to run more ADMM
    while("try-error" %in% class(par.nr.out) | par.iter == 10){
      num.iter = num.iter * 2
      admm.fit = onestepADMM(x = x, y = y, censor = censor, 
                             family = family, tol = tolerance, 
                             maxit = num.iter, lambda = lambda,
                             intercept = intercept, C = matrix(C, nrow = nrow(C)), d = d, 
                             penalty.factor = penalty.factor, 
                             initial.guess = admm.fit$solution[,2],
                             delta_0 = admm.fit$delta.out, u1_0 = admm.fit$u1.out,
                             u2_0 = admm.fit$u2.out, z_0 = admm.fit$z.out, hessian = admm.fit$hessian)
      
      if(admm.fit$iteration[2] != 2 * num.iter){
        stop.nr.soln = T
        break
      }
      
      par.nr = admm.fit$z.out
      nonzero.index = which(par.nr != 0)
      constrained.nonzero = constrained.index[which(constrained.index %in% nonzero.index)]
      nonzero.parameters = par.nr[constrained.nonzero]
      
      if(length(constrained.nonzero) > 0){
        C_nonzero = C[,constrained.nonzero, drop = F]
        rank_condition = (rankMatrix(C_nonzero) == nrow(C_nonzero))
        
        #To find a feasible solution for N-R, we need to project into a potentially affine space
        #This requires a matrix inversion, and we need a constraint matrix with full rank
        if(rank_condition == F){
          while(rank_condition == F){
            num.iter = num.iter * 2
            admm.fit = onestepADMM(x = x, y = y, censor = censor, 
                                   family = family, tol = tolerance, 
                                   maxit = num.iter, lambda = lambda,
                                   intercept = intercept, C = matrix(C, nrow = nrow(C)), d = d,
                                   penalty.factor = penalty.factor, 
                                   initial.guess = admm.fit$solution[,2],
                                   delta_0 = admm.fit$delta.out, u1_0 = admm.fit$u1.out,
                                   u2_0 = admm.fit$u2.out, z_0 = admm.fit$z.out, hessian = admm.fit$hessian)
            
            par.nr = admm.fit$z.out
            nonzero.index = which(par.nr != 0)
            constrained.nonzero = constrained.index[which(constrained.index %in% nonzero.index)]
            C_nonzero = C[,constrained.nonzero, drop = F]
            rank_condition = (rankMatrix(C_nonzero) == nrow(C_nonzero))
          }
        }
        
        
        
        par_nonzero = par.nr[constrained.nonzero]
        par_nonzero = (diag(nrow = length(constrained.nonzero)) - t(C_nonzero) %*% solve(C_nonzero %*% t(C_nonzero)) %*%
                         C_nonzero) %*% par_nonzero + t(C_nonzero) %*% 
          solve(C_nonzero %*% t(C_nonzero)) %*% d
        par.nr[constrained.nonzero] = par_nonzero
        
      }
      
      par.nr.out = try(NR_constrained(par.nr[nonzero.index], xdata = x[,nonzero.index], 
                                      censor = censor, lambda = lambda, 
                                      ydata = y, family = family,
                                      penalty.factor = penalty.factor[nonzero.index],
                                      C = C[,nonzero.index], tol = tolerance), silent = T)
      
      if(!("try-error" %in% class(par.nr.out)) == T){
        par.iter = par.nr.out$iter
        par.nr.inner = par.nr.out$solution
      }
    }
    if(stop.nr.soln == T){break}
    #our guess is the N-R solution, plus 0s for terms not used in N-R
    new_guess = rep(0, ncol(x))
    new_guess[nonzero.index] = par.nr.inner
    
    #checking KKT conditions for the inner subset
    kkt.logical = KKTtestfct(new_guess, family = family, lambda = lambda, 
                             x = x, y = y, C = C, censor = censor, trace = F)
    
    #incrementing the number of ADMM iterations of the next attempt
    num.iter = num.iter * 2
    if(num.iter > 100000){
      break
    }
  }
  if(stop.zero.soln == T){
    new_guess = rep(0, ncol(x))
  } else if (stop.one.soln == T){
    new_guess = admm.fit$z.out
  } else if (stop.rank.soln == T){
    new_guess = admm.fit$z.out
  } else if (stop.nr.soln == T){
    new_guess = admm.fit$z.out
  }
  
  if(trace == F){
    return(par = new_guess)
  } else{
    return(list(par = new_guess, numiter = admm.fit$iteration))
  }
}
