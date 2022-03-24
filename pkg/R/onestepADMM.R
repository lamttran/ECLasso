#' ADMM function
#'
#' This function allows you to use ADMM to fit the constrained lasso
#' @param x Design matrix, data for the predictors
#' @param y Response vector
#' @param family The type of model to fit; at this point,
#' only logistic and cox models have been updated with the subset approach outlined in the paper
#' @param rho proximal mapping parameter for z-update in ADMM
#' @param censor If a cox model is fit, the censoring vector for the survival data
#' @param C the vector that linearly constrains the regression coefficients, e.g. a 1 by p
#' vector constrains the sum of the coefficients
#' @param d the value that the constrained regression coefficients should equal e.g. the above
#' constraint vector and d = 0 would be a sum-to-zero constraint
#' @param inexact Whether ADMM should be fit with an inexact speedup, required for all non-gaussian models
#' @param initial.guess initial regression coefficient guess
#' @param tol tolerance for ADMM, determining primal and dual conditions for termination
#' @param maxit maximum ADMM iterations, required do the slowness of ADMM for high accuracy convergence
#' @import glmnet
#' @import survival
#' @import MASS
#' @param lambda Lasso penalty parameter
#' @param intercept Whether an unpenalized, unconstrained intercept should be included in the model
#' @param min.abs The minimum value such that terms smaller in absolute magnitude are shrunk to 0
#' @param lower.limits If values are less than lower.limits, they are set to lower.limits
#' @param upper.limits If values are greater than upper.limits, they are set to upper.limits
#' @param penalty.factor The vector that multiplies the lasso penalty, allowing for different penalties
#' for each term e.g. 1 is penalization equal to lambda and 0 is no penalization
#' @param delta_0 output from prior ADMM as new start point for rerun ADMM
#' @param u1_0 output from prior ADMM as new start point for rerun ADMM
#' @param u2_0 output from prior ADMM as new start point for rerun ADMM
#' @param z_0 output from prior ADMM as new start point for rerun ADMM
#' @param hessian output from prior ADMM as new start point for rerun ADMM
#' @param exclude predictors to exclude from the model, coefficient forced to equal 0
#' @param trace if true, also outputs time to fit, solution at each ADMM iteration, Hessian matrix, and
#' signs of the regression coefficients
#' @keywords ADMM
#' @export
#' @examples
#' onestepADMM()

onestepADMM <- function(x, y, family = c("gaussian", "cox", "poisson", "logistic"),
                        rho = 1, censor = NULL, C = NULL, d = NULL, inexact = TRUE, initial.guess = NULL,
                        tol = 1e-7, maxit = 1000, lambda = NULL, intercept = F,
                        min.abs = 0, lower.limits = -Inf, upper.limits = Inf,
                        penalty.factor = 1, delta_0 = NULL, u1_0 = NULL, u2_0 = NULL,
                        z_0 = NULL, hessian = NULL, exclude = rep(FALSE, nvars), trace = T){

  family = match.arg(family)

  abs.tol = ifelse(inexact, tol * 1e-4, tol * 1e-2) #may have to change this
  rel.tol = ifelse(inexact, tol * 1e-2, tol)
  max.iter = ifelse(inexact, maxit * 2, maxit)

  A = as.matrix(x)
  b = y

  nobs = nrow(A)
  nvars = ncol(A)

  if (is.null(d)) {d = 0}
  if (is.null(C)) {
    #C is not specified, set it to a row of 1.
    C = matrix(1, nrow = 1, ncol = nvars)
  }

  if (length(min.abs) == 1) min.abs = rep(min.abs, nvars)
  if (length(lower.limits) == 1) lower.limits = rep(lower.limits, nvars)
  if (length(upper.limits) == 1) upper.limits = rep(upper.limits, nvars)
  if (length(penalty.factor) == 1) penalty.factor = rep(penalty.factor, nvars)

  if (!is.logical(exclude)) {

    exc = rep(F, nvars)
    exc[exclude] = T
    exclude = exc

  }

  if (intercept == TRUE) {

    A = cbind(1, A)
    C = cbind(0, C)
    min.abs = c(0, min.abs)
    lower.limits = c(-Inf, lower.limits)
    upper.limits = c(Inf, upper.limits)
    penalty.factor = c(0, penalty.factor)
    exclude = c(F, exclude)

  }

  p = ncol(A)
  if(is.null(initial.guess)){
    x = as.matrix(rep(0, p), nrow = p, ncol = 1)
  } else{
    x = as.matrix(initial.guess, nrow = p, ncol = 1)
  }

  if (is.null(z_0)){
    z = matrix(0, nrow = p, ncol = 1)
  } else{
    z = z_0
  }


  if (is.null(u1_0)){
    u1 = 0
  } else{
    u1 = u1_0
  }

  if(is.null(u2_0)){
    u2 = matrix(0, nrow = p, ncol = 1)
  } else{
    u2 = u2_0
  }

  if(is.null(hessian)){

    if(family == "gaussian"){
      Atb = t(A) %*% b

      if (inexact == TRUE) {

        J = t(A) %*% A

      } else {

        E = rbind(A, sqrt(rho) * C)
        LU = Choleski_factors(E, rho)
        L = LU[[1]]
        #m*m lower triangular
        U = LU[[2]]

      }
    } else if(family == "cox"){

      ord = order(b)
      b1 = b[ord]
      censor1 = censor[ord]
      A1 = as.matrix(A[ord, ])

      hazard = as.numeric(exp(A1 %*% x))
      risk = c(cumsum_rev(hazard))
      P = outer(hazard, risk, '/')
      P[upper.tri(P)] = 0
      W = -P %*% diag(censor1) %*% t(P)
      diag(W) = diag(P %*% diag(censor1) %*% t(1 - P))
      H = t(A1) %*% W %*% A1

      #J = -H
      J = H

      ##new stuff here, sort by decreasing again like old algorithm
      ord = order(b, decreasing = T)
      b1 = b[ord]
      censor1 = censor[ord]
      A1 = as.matrix(A[ord, ])

      ##
    } else if (family == "poisson") {
      J = t(A) %*% diag(exp(A %*% x)[, 1]) %*% A;
    } else if (family == "logistic") {
      #Func = -diag(c(b)) %*% A
      #J = t(Func) %*% diag((exp(Func %*% x) / (1 + exp(Func %*% x)) ^ 2)[, 1]) %*% Func
      pi = 1/(1 + exp(-A %*% x))
      J = t(A) %*% diag(c(pi * (1 - pi))) %*% A
    }

  } else{
    J = hessian
    if(family == "cox"){
      ord = order(b, decreasing = T)
      A1 = as.matrix(A[ord, ])
      censor1 = censor[ord]
    }
  }

  if(inexact == T){
    eigenvalues = eigen(J)$values
    h = max(eigenvalues)
    P = 1 / (h + rho) * (diag(p) - rho / (h + rho * (p + 1)) * t(C) %*% C)
  }



  if(trace == T){
    z_iterate = matrix(NA, nrow = p, ncol = max.iter)
    solution_iterate = matrix(NA, nrow = p, ncol = max.iter)
    signs = matrix(NA, nrow = p, ncol = max.iter)
    time_iterate = matrix(NA, nrow = 5, ncol = ceiling(max.iter / 100))
  }

  iter_num = c()
  solution = c()

  iter_num = c(iter_num, 1)
  solution = cbind(solution, x)


  ptm = proc.time()

  ##delta part of the term to update beta equals (E'E + rho * I)^-1(x'y + delta)
  for(k in 1:max.iter){
    if(k == 1 & is.null(delta_0) == F){
      delta = delta_0
    } else{
      delta = rho * (z - u2 + t(C) %*% (d - u1))
    }


    if(family == "gaussian"){
      if(inexact == T){
        x = P %*% (delta + h * x - J %*% x + Atb)
      } else {
        q = Atb + delta
        x = backsolve(U, forwardsolve(L, q))
      }
    }

    if(family == "cox"){

      thetaIT = exp(A1 %*% x)
      A2 = as.matrix(A1) * as.vector(thetaIT)
      numeratorIT = cumsum_col(A2)
      denominatorIT = cumsum(thetaIT)
      w = c(Arma_colSums((A1 - numeratorIT/denominatorIT) * censor1))
      x = P %*% (h * x + delta + w)
    }

    if(family == "poisson"){
      x = P %*% (h * x + delta + t(A) %*% b - t(A) %*% exp(A %*% x))
    }

    if(family == "logistic"){
      #x = P %*% (h * x + delta - t(Func) %*% (exp(Func %*% x) / (1 + exp(Func %*% x))))
      pi = 1/(1 + exp(-A %*% x))
      x = P %*% (h * x + delta - t(A) %*% (pi - b))

    }

    x[abs(x) < min.abs] = 0
    x[x < lower.limits] = lower.limits
    x[x > upper.limits] = upper.limits
    x[exclude] = 0

    # z-update
    zold = z
    z = x + u2

    z = shrinkage(z, lambda * penalty.factor / rho)

    #u-update
    Cx = as.numeric(C %*% x)
    u1 = u1 + Cx - d
    u2 = u2 + (x - z)


    if(k %% 100 == 0){
      time_iterate[, (k / 100)] = proc.time() - ptm
    }

    if(trace == T){
      z_iterate[, k] = z
      solution_iterate[, k] = x
      signs[, k] = sign(z)
    }

    #diagnostics, reporting, termination checks
    r_norm = vector.2.norm(rbind(Cx - d, x - z))
    #r_norm
    s_norm = vector.2.norm( - rho * (z - zold))
    #s_norm
    eps_pri = sqrt(p + 1) * abs.tol + rel.tol *
      max(vector.2.norm(rbind(Cx, x)), vector.2.norm( - z), abs(d))
    #eps_pri
    eps_dual = sqrt(p) * abs.tol + rel.tol * vector.2.norm(rho * (u1 + u2))
    #eps_dual

    #termination criterion
    if (r_norm < eps_pri && s_norm < eps_dual) {
      break
    }
  }

  iter_num = c(iter_num, k)

  solution = cbind(solution, t(t(x)))
  ##rounding off terms bc you can get stuff that's of order -8 or even smaller
  if(trace == F){
    return(list(iteration = iter_num, solution = solution,
                delta.out = delta, u1.out = u1, u2.out = u2, z.out = z))
  } else {
    return(list(iteration = iter_num, solution = solution,
                delta.out = delta, u1.out = u1, u2.out = u2, z.out = z,
                z_iterate = z_iterate, signs = signs,
                solution_iterate = solution_iterate,
                time_iterate = time_iterate, hessian = J))
  }
}
