library(ECLasso)
library(intasymmRcppArma)
library(MASS)

#Simulate some simple logistic and survival data
set.seed(13)
n = 100 #number of samples
depth = 20 #number of lambda values to consider from unconstrained sequence
beta = c(1, -0.8, 0.4, 0, -0.6, rep(0, 45)) #True sparse parameter vector
p = length(beta)

C = matrix(1, ncol = 50)
d = matrix(0) #A sum-to-zero constraint on regression coefficients

x.mean = rnorm(p)
corr.mat = zapsmall(0.3 ^ abs(outer(1:p, 1:p, FUN = "-")))
x = mvrnorm(n, x.mean, corr.mat) #Simulated design matrix with correlation

#Simulation of binary response for logistic regression
eta = x %*% beta + 0.5*rnorm(n)
y = rbinom(n, 1, 1/(1 + exp(-eta)))
rm(x.mean, eta, corr.mat)

#Simulation of survival data for Cox regression
times = rexp(n, exp(as.numeric(x %*% beta)))
time.censor = rexp(n, exp(as.numeric(x %*% beta) + x[, 1] * beta[1]))
censorv = ifelse(times < time.censor, 1, 0) #censoring indicator
y.cox <- ifelse(times < time.censor, times, time.censor) #observed time
ord = order(y.cox, decreasing = F)
rm(times, time.censor)

x.cox = x[ord,]
times_final = y.cox[ord]
censor_final = censorv[ord]
rm(ord, censorv)

#Calling ECLasso to solve the constrained lasso
sim_logistic = ECLasso(x = x, y = y, C = C, d = d, family = "logistic", 
                       depth = depth, trace = T, method = "combined_subset")

sim_cox = ECLasso(x = x.cox, y = times_final, censor = censor_final, C = C, d = d,
                  family = "cox", depth = depth, trace = T, method = "combined_subset")


#Plots of the logistic and Cox solution paths
#Red lines denote true positive predictors, blue lines true negative predictors
plotlambdas <- plot(log(sim_logistic$lambda_seq)[1:ncol(sim_logistic$solution)], sim_logistic$solution[1, ], type = "l",
                    xlim = c(log(sim_logistic$lambda_seq[ncol(sim_logistic$solution)]), log(sim_logistic$lambda_seq[1])),
                    ylim = c(-1.2 * abs(min(sim_logistic$solution)),
                             1.2 * abs(max(sim_logistic$solution))), 
                    main = "Logistic solution path", xlab = "Log lambda",
                    ylab = "Param. Estimate", col = "red", cex.lab = 1.5, 
                    cex.axis = 1.5, cex.main = 2, cex.sub = 1.5, lwd = 1.5)

for (k in 2 : nrow(sim_logistic$solution)){
  lines(log(sim_logistic$lambda_seq)[1:ncol(sim_logistic$solution)], sim_logistic$solution[k, ], 
        col = ifelse(beta[k] != 0, ifelse(beta[k] > 0, "red", "blue"), "#0000002D"), lwd = 1.5)
}

abline(h = 0, lwd = 2)


plotlambdas <- plot(log(sim_cox$lambda_seq)[1:ncol(sim_cox$solution)], sim_cox$solution[1, ], type = "l",
                    xlim = c(log(sim_cox$lambda_seq[ncol(sim_cox$solution)]), log(sim_cox$lambda_seq[1])),
                    ylim = c(-1.2 * abs(min(sim_cox$solution)),
                             1.2 * abs(max(sim_cox$solution))), 
                    main = "Cox solution path", xlab = "Log lambda",
                    ylab = "Param. Estimate", col = "red", cex.lab = 1.5, 
                    cex.axis = 1.5, cex.main = 2, cex.sub = 1.5, lwd = 1.5)

for (k in 2 : nrow(sim_cox$solution)){
  lines(log(sim_cox$lambda_seq)[1:ncol(sim_cox$solution)], sim_cox$solution[k, ], 
        col = ifelse(beta[k] != 0, ifelse(beta[k] > 0, "red", "blue"), "#0000002D"), lwd = 1.5)
}

abline(h = 0, lwd = 2)