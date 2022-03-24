#If you do not have the spatstat.utils package, uncomment the following line and install it
#install.packages("spatstat.utils")

library(MASS)
library(glmnet)
library(survival)
library(spatstat.utils)
set.seed(2024) #2024 scenario 1, 2025 scenario 2, 2026 scenario 3

n = 100
beta = c(1, -0.8, 0.4, 0, -0.6, rep(0, 45)) #for scenario 1
#beta = c(0.3, 0.6, -0.5, -0.1, 0.5, 0.6, 0, 0.2, 0.7, -0.3, rep(0, 90)) # for scenario 2
#beta = c(0.4, 0.45, 0.15, 0.3, 0.2, rep(0, 95)) #for scenario 3
p = length(beta)
x.mean = rnorm(p)

corr.mat = zapsmall(0.3 ^ abs(outer(1:p, 1:p, FUN = "-")))
x = mvrnorm(n, x.mean, corr.mat)

eta = x %*% beta + 0.5*rnorm(n)
y.logistic = rbinom(n, 1, 1/(1 + exp(-eta)))

times = rexp(n, exp(as.numeric(x %*% beta)))
time.censor = rexp(n, exp(as.numeric(x %*% beta) + x[,1] * beta[1]))
censorv = ifelse(times < time.censor, 1, 0) #censoring indicator
y.cox <- ifelse(times < time.censor, times, time.censor) #observed time
ord = order(y.cox, decreasing = F)
rm(times, time.censor)

x.cox = x[ord,]
times_final = y.cox[ord]
censor_final = censorv[ord]
rm(ord, censorv)

penalty.factor.full = rep(1, p)
constraint.full = rep(1, p) #for scenarios 1 and 2
#constraint.full = c(rep(1, 3), rep(2, 97)) #for scenario 3

d = 0 #for scenario 1
#d = 1 #for scenario 2
#d = 2 #for scenario 3
##fit for logistic
tolerance = 1e-7
rho = 1

fit1 = glmnet(x = x, y = y.logistic, family = "binomial", intercept = F)
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


for(i in 1:20){
  stop = F
  ptm <- proc.time()
  
  lambda = lambda_seq[i] 
  
  
  index_i = par_keep[[unique_index[min(i.glmnet, nlambda)]]] - 1
  while(length(index_i) < 1){
    i.glmnet = i.glmnet + 1
    index_i = par_keep[[unique_index[min(i.glmnet, nlambda)]]]
  }
  
  #xdata.new = x[, index_i]
  xdata.new = as.matrix(x[, index_i], ncol = length(index_i))
  

  fit.admm1 = onestepADMM(x = xdata.new, y = y.logistic,
                          family = "logistic", tol = tolerance, 
                          maxit = 100, lambda = lambda,
                          intercept = F, C = matrix(constraint.full[index_i], nrow = 1), d = d,
                          penalty.factor = constraint.full[index_i])
  
  
  new_guess = rep(0, p)
  
  new_guess[index_i] = inner.solution(admm.fit = fit.admm1, family = "logistic", x = xdata.new, 
                                      y = y.logistic, intercept = F,
                                      lambda = lambda, C = constraint.full[index_i], 
                                      penalty.factor = penalty.factor.full[index_i], d = d)
  
  
  kkt.outer = KKTtestfct(par_guess = new_guess, family = "logistic", lambda = lambda, 
                         x = x, y = y.logistic, C = constraint.full)
  
  while(kkt.outer == F){
    i.glmnet = i.glmnet + 1
    if(i.glmnet > nlambda){
      stop == TRUE
      break
    }
    index_i = par_keep[[unique_index[min(i.glmnet, nlambda)]]] - 1

    #xdata.new = x[, index_i]
    xdata.new = as.matrix(x[, index_i], ncol = length(index_i))
    
    fit.admm1 = onestepADMM(x = xdata.new, y = y.logistic,
                            family = "logistic", tol = tolerance, 
                            maxit = 200, lambda = lambda,
                            intercept = F, C = matrix(constraint.full[index_i], nrow = 1), d = d,
                            penalty.factor = constraint.full[index_i])
    
    new_guess = rep(0, p)
    new_guess[index_i] = inner.solution(admm.fit = fit.admm1, family = "logistic", x = xdata.new, 
                                        y = y.logistic, intercept = F,
                                        lambda = lambda, C = constraint.full[index_i], 
                                        penalty.factor = penalty.factor.full[index_i], d = d)
    
    kkt.outer = KKTtestfct(par_guess = new_guess, family = "logistic", lambda = lambda, 
                           x = x, y = y.logistic, C = constraint.full)
  }
  if(stop == T){break}
  
  solution = cbind(solution, new_guess)
  time.iter = c(time.iter, (proc.time() - ptm)[3])
  i.glmnet.iter = c(i.glmnet.iter, i.glmnet)
}

solnNorm <- (abs(solution[, ncol(solution)]) - min(abs(solution[, ncol(solution)]))) / 
  (max(abs(solution[, ncol(solution)])) - min(abs(solution[, ncol(solution)])))

##for full plot
plotlambdas <- plot(log(lambda_seq)[1:ncol(solution)], solution[1, ], type = "l",
                    xlim = c(log(lambda_seq[ncol(solution)]), log(lambda_seq[1])),
                    ylim = c(-1.2 * max(abs(solution)),
                             1.2 * max(abs(solution))), xlab = "Log lambda",
                    ylab = "Param. Estimate", col = "red", main = "Logistic constrained model, scenario 3")

for (k in 2 : nrow(solution)){
  lines(log(lambda_seq)[1:ncol(solution)], solution[k, ], col = ifelse(beta[k] != 0, 
                                                                       ifelse(beta[k] > 0, "red", "blue"), "black"))
  #col = addalpha(rep("red4", 100), solnNorm[k])
}


#fit for cox
tolerance = 1e-7
rho = 1

fit1 = glmnet(x = x.cox, y = Surv(times_final, censor_final), family = "cox", intercept = F)
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

for(i in 1:20){
  stop = F
  ptm <- proc.time()
  
  lambda = lambda_seq[i] 
  
  
  index_i = par_keep[[unique_index[min(i.glmnet, nlambda)]]]
  while(length(index_i) < 2){
    i.glmnet = i.glmnet + 1
    index_i = par_keep[[unique_index[min(i.glmnet, nlambda)]]]
  }
  
  #xdata.new = x[, index_i]
  xdata.new = as.matrix(x.cox[, index_i], ncol = length(index_i))
  
  
  fit.admm1 = onestepADMM(x = xdata.new, y = times_final, censor = censor_final,
                          family = "cox", tol = tolerance, 
                          maxit = 100, lambda = lambda,
                          intercept = F, C = matrix(constraint.full[index_i], nrow = 1), d = d,
                          penalty.factor = constraint.full[index_i])
  
  
  new_guess = rep(0, p)
  
  new_guess[index_i] = inner.solution(admm.fit = fit.admm1, family = "cox", x = xdata.new, 
                                      y = times_final, censor = censor_final, intercept = F,
                                      lambda = lambda, C = matrix(constraint.full[index_i], nrow = 1), 
                                      penalty.factor = constraint.full[index_i], d = d)
  
  
  kkt.outer = KKTtestfct(par_guess = new_guess, family = "cox", lambda = lambda, 
                         x = x.cox, y = times_final, censor = censor_final,
                         C = matrix(constraint.full, nrow = 1))
  
  while(kkt.outer == F){
    i.glmnet = i.glmnet + 1
    if(i.glmnet > nlambda){
      stop == TRUE
      break
    }
    index_i = par_keep[[unique_index[min(i.glmnet, nlambda)]]]
    
    #xdata.new = x[, index_i]
    xdata.new = as.matrix(x.cox[, index_i], ncol = length(index_i))
    
    fit.admm1 = onestepADMM(x = xdata.new, y = times_final, censor = censor_final,
                            family = "cox", tol = tolerance, 
                            maxit = 200, lambda = lambda,
                            intercept = F, C = matrix(constraint.full[index_i], nrow = 1), d = d,
                            penalty.factor = constraint.full[index_i])
    
    new_guess = rep(0, p)
    new_guess[index_i] = inner.solution(admm.fit = fit.admm1, family = "cox", x = xdata.new, 
                                        y = times_final, censor = censor_final, intercept = F,
                                        lambda = lambda, C = matrix(constraint.full[index_i], nrow = 1), 
                                        penalty.factor = constraint.full[index_i], d = d)
    
    kkt.outer = KKTtestfct(par_guess = new_guess, family = "cox", lambda = lambda, 
                           x = x.cox, y = times_final, censor = censor_final,
                           C = matrix(constraint.full, nrow = 1))
  }
  
  if(stop == T){break}
  
  solution = cbind(solution, new_guess)
  time.iter = c(time.iter, (proc.time() - ptm)[3])
  i.glmnet.iter = c(i.glmnet.iter, i.glmnet)
}

solnNorm <- (abs(solution[, ncol(solution)]) - min(abs(solution[, ncol(solution)]))) / 
  (max(abs(solution[, ncol(solution)])) - min(abs(solution[, ncol(solution)])))

##for full plot
plotlambdas <- plot(log(lambda_seq)[1:ncol(solution)], solution[1, ], type = "l",
                    xlim = c(log(lambda_seq[ncol(solution)]), log(lambda_seq[1])),
                    ylim = c(-1.2 * max(abs(solution)),
                             1.2 * max(abs(solution))), xlab = "Log lambda",
                    ylab = "Param. Estimate", col = "red", main = "Cox constrained model, scenario 3")

for (k in 2 : nrow(solution)){
  lines(log(lambda_seq)[1:ncol(solution)], solution[k, ], col = ifelse(beta[k] != 0, 
                                                                       ifelse(beta[k] > 0, "red", "blue"), "black"))
}
