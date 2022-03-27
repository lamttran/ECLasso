library(MASS)
library(glmnet)
library(survival)
library(spatstat.utils)
set.seed(1111) #2024 scenario 1, 2025 scenario 2, 2026 scenario 3

n = 100
beta = c(1, -0.8, 0.4, 0, -0.6, rep(0, 45)) #for scenario 1
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
constraint.full = rep(1, p)

#a sum-to-zero constrained lasso example for a logistic model w/ binary data
fast_constrained(x = x, y = y.logistic, family = "logistic")

#a sum-to-zero constrained lasso example for a cox w/ survival data 
fast_constrained(x = x.cox, y = times_final, censor = censor_final, family = "cox")
