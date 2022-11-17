library(ECLasso)
library(intasymmRcppArma)
library(MASS)


n_simulations = 20 
depth = 20
time.matrix = matrix(nrow = n_simulations * 3, ncol = depth)

#Un-comment lines according to the simulation scenario to replicate

n = 100 #Scenarios 1-3
#n = 75 #Scenario 4

#True parameter vectors for simulation scenarios 1-4

#Scenario 1
beta = c(1, -0.8, 0.4, 0, -0.6, rep(0, 45))

#Scenario 2
#beta = c(0.5, -0.4, 0.2, 0, -0.3, rep(0, 45)) 

#Scenario 3
#beta = c(0.1, 0.35, -0.2, 0.5, -0.4, -0.3, 0.4, 0.7, rep(0, 92))

#Scenario 4
#beta = c(1, -0.8, 0.4, 0, -0.6, rep(0, 95))

p = length(beta)

##Constraint matrices for simulation scenarios 1-4

#Scenario 1
C = matrix(1, ncol = 50)
d = matrix(0)

#Scenario 2
#C = matrix(rbind(c(rep(1, p)), c(1, 1, 2, rep(0, 47))), nrow = 2)
#d = matrix(c(0, 0.5), ncol = 1)

#Scenario 3
#C = matrix(rbind(c(rep(1, p)), c(0, 0, 1, 0, 1, rep(0, 95))), nrow = 2)
#d = matrix(c(1, -1), ncol = 1)

#Scenario 4
#C = matrix(1, ncol = 100)
#d = matrix(0)

for(k in 1:n_simulations){
  print(k)
  set.seed(k + 2023)
  time.vector = matrix(nrow = 3, ncol = depth)
  
  x.mean = rnorm(p)
  corr.mat = zapsmall(0.3 ^ abs(outer(1:p, 1:p, FUN = "-")))
  x = mvrnorm(n, x.mean, corr.mat)
  
  
  
  #logistic code
  eta = x %*% beta + 0.5*rnorm(n)
  y = rbinom(n, 1, 1/(1 + exp(-eta)))
  rm(x.mean, eta, corr.mat)

  testhold = ECLasso(x = x, y = y, C = C, d = d, family = "logistic", 
                     depth = depth, method = "combined_subset", trace = T)
  time.vector[1, ] = as.vector(cumsum(testhold$time))

  testhold2 = ECLasso(x = x, y = y, C = C, d = d, family = "logistic", 
                      depth = depth, method = "admm_subset", trace = T)
  time.vector[2, ] = as.vector(cumsum(testhold2$time))

  testhold3 = ECLasso(x = x, y = y, C = C, d = d, family = "logistic", 
                      depth = depth, method = "naive", trace = T)
  time.vector[3, ] = as.vector(cumsum(testhold3$time))
  
  #cox code
  # times = rexp(n, exp(as.numeric(x %*% beta)))
  # time.censor = rexp(n, exp(as.numeric(x %*% beta) + x[, 1] * beta[1]))
  # censorv = ifelse(times < time.censor, 1, 0) #censoring indicator
  # y.cox <- ifelse(times < time.censor, times, time.censor) #observed time
  # ord = order(y.cox, decreasing = F)
  # rm(times, time.censor)
  # 
  # x.cox = x[ord,]
  # times_final = y.cox[ord]
  # censor_final = censorv[ord]
  # rm(ord, censorv)
  # 
  # testhold = ECLasso(x = x.cox, y = times_final, censor = censor_final, C = C, d = d,
  #                             family = "cox", depth = depth, method = "combined_subset", trace = T)
  # time.vector[1, ] = as.vector(cumsum(testhold$time))
  # 
  # testhold2 = ECLasso(x = x.cox, y = times_final, censor = censor_final, C = C, d = d,
  #                        family = "cox", depth = depth, method = "admm_subset", trace = T)
  # time.vector[2, ] = as.vector(cumsum(testhold2$time))
  # 
  # testhold3 = ECLasso(x = x.cox, y = times_final, censor = censor_final, C = C, d = d,
  #                       family = "cox", depth = depth, method = "naive", trace = T)
  # time.vector[3, ] = as.vector(cumsum(testhold3$time))
  
  
  equal12 = all.equal(round(testhold$solution[, 10], 4), round(testhold2$solution[, 10], 4))
  equal23 = all.equal(round(testhold3$solution[, 10], 4), round(testhold2$solution[, 10], 4))
  if(equal12 != T | equal23 != T){
    break
  }
  
  time.matrix[(3 * k - 2):(3 * k), ] = time.vector
}

time.df = data.frame(time.matrix)
time.df = cbind(rep(c("Combined Subset", "ADMM Subset", "ADMM Naive"), n_simulations), time.df)
colnames(time.df) = c("Method", paste0(1:depth))