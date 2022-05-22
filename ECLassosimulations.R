##Note, red lines in plot represent solution path for predictors w/ true positive parameter
##Blue lines represent path for predictors w/ true negative parameter
##Gray lines represent path for noise predictors


library(MASS)
library(glmnet)
library(survival)
library(spatstat.utils)
library(Matrix)
set.seed(118)  

n = 100

beta = c(0.5, -0.4, 0.2, 0, -0.3, rep(0, 45)) 
p = length(beta)
x.mean = rnorm(p)

corr.mat = zapsmall(0.3 ^ abs(outer(1:p, 1:p, FUN = "-")))
x = mvrnorm(n, x.mean, corr.mat)

eta = x %*% beta + 0.5*rnorm(n)
y.logistic = rbinom(n, 1, 1/(1 + exp(-eta)))
rm(x.mean, eta, corr.mat)

C = matrix(rbind(c(rep(1, p)), c(1, 1, 2, rep(0, 47))), nrow = 2)
d = matrix(c(0, 0.5), ncol = 1)

y = y.logistic

##Cox terms
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

#Constrained logistic model
logistic_fit = ECLasso(x = x, y = y.logistic, C = C, d = d, family = "logistic", depth = 20, trace = T)
zapsmall(C%*%logistic_fit$solution)

solution = logistic_fit$solution
lambda_seq = logistic_fit$lambda_seq

solnNorm <- (abs(solution[, ncol(solution)]) - min(abs(solution[, ncol(solution)]))) / 
  (max(abs(solution[, ncol(solution)])) - min(abs(solution[, ncol(solution)])))

plotlambdas <- plot(log(lambda_seq)[1:ncol(solution)], solution[1, ], type = "l",
                    xlim = c(log(lambda_seq[ncol(solution)]), log(lambda_seq[1])),
                    ylim = c(-1.2 * abs(min(solution)),
                             1.2 * abs(max(solution))), xlab = "Log lambda",
                    ylab = "Param. Estimate", col = "red", cex.lab=1.5, 
                    cex.axis=1.5, cex.main=2, cex.sub=1.5, lwd = 1.5)

for (k in 2 : nrow(solution)){
  lines(log(lambda_seq)[1:ncol(solution)], solution[k, ], 
        col = ifelse(beta[k] != 0, ifelse(beta[k] > 0, "red", "blue"), "#0000002D"), lwd = 1.5)
  #col = addalpha(rep("black", 100), solnNorm[k])
}

abline(h = 0, lwd = 2)

#Constrained Cox model
cox_fit = ECLasso(x = x.cox, y = times_final, censor = censor_final, C = C, d = d,
                             family = "cox", depth = 20, trace = T)

solution = cox_fit$solution
lambda_seq = cox_fit$lambda_seq

solnNorm <- (abs(solution[, ncol(solution)]) - min(abs(solution[, ncol(solution)]))) / 
  (max(abs(solution[, ncol(solution)])) - min(abs(solution[, ncol(solution)])))

plotlambdas <- plot(log(lambda_seq)[1:ncol(solution)], solution[1, ], type = "l",
                    xlim = c(log(lambda_seq[ncol(solution)]), log(lambda_seq[1])),
                    ylim = c(-1.2 * abs(min(solution)),
                             1.2 * abs(max(solution))), xlab = "Log lambda",
                    ylab = "Param. Estimate", col = "red",  cex.lab=1.5, cex.axis=1.5,
                    cex.main=2, cex.sub=1.5, lwd = 1.5)

for (k in 2 : nrow(solution)){
  lines(log(lambda_seq)[1:ncol(solution)], solution[k, ], 
        col = ifelse(beta[k] != 0, ifelse(beta[k] > 0, "red", "blue"), "#0000002D"), lwd = 1.5)
}

abline(h = 0, lwd = 2)
