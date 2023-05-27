library(glmnet)

##logistic outcome of interest: aapperio4_bin
##Load in data below
load("~/ORIGINS_processed_data.Rdata")

##Remove missing data
OTU.new = OTU.mat
OTU.new = OTU.new[-which(is.na(resp.mat$age) | is.na(resp.mat$sex) | is.na(resp.mat$aapperio4_bin) |
                           is.na(resp.mat$raceethn) | is.na(resp.mat$bmi)),]
resp.mat = resp.mat[-which(is.na(resp.mat$age) | is.na(resp.mat$sex) | is.na(resp.mat$aapperio4_bin) |
                             is.na(resp.mat$raceethn) | is.na(resp.mat$bmi)),]

##Re-categorize race variable

resp.mat$hispanic = (resp.mat$raceethn == 0)
resp.mat$black = (resp.mat$raceethn == 2)
resp.mat$other.race = (resp.mat$raceethn == 3)

##Remove taxa with more than 80% 0 counts
##Add a small increment as in microbiome literature
col.zeroes = colSums(OTU.new == 0)
zero.cutoff = floor(nrow(resp.mat) * 0.8)
cols.remove = which(col.zeroes > zero.cutoff)
OTU.new = OTU.new[,-cols.remove]

rm(col.zeroes, zero.cutoff, cols.remove)

OTU.new = OTU.new + 0.5

##compositionalize the data, log transform, and center scale the taxa
for(i in 1:nrow(OTU.new)){
  OTU.new[i,] = OTU.new[i,]/rowSums(OTU.new)[i]
}

OTU.new = log(OTU.new)
OTU.new = cbind(OTU.new, resp.mat$age, resp.mat$sex, 
                resp.mat$hispanic, resp.mat$black, 
                resp.mat$other.race, resp.mat$bmi)

for(i in 1:268){
  OTU.new[, i] = OTU.new[, i] - mean(OTU.new[, i])
}


##Scale age and BMI, change scale of interpretation to per 10 unit change
OTU.new[,269] = (OTU.new[,269] - mean(OTU.new[,269]))/10
OTU.new[,274] = (OTU.new[,274] - mean(OTU.new[,274]))/10

##Tolerance and step-size hyper-parameter
tolerance = 1e-4
rho = 1 

##The paper plots used depth = 25, but depth = 15 results in a much faster fit
taxa_fit = ECLasso(x = OTU.new, y = resp.mat$aapperio4_bin, 
                   family = "logistic", penalty.factor = c(rep(1, 268), rep(0, 6)), 
                   C = matrix(c(rep(1, 268), rep(0, 6)), nrow = 1),
                   method = "combined_subset", intercept = T, tolerance = tolerance,
                   depth = 25, trace = T)


solution = taxa_fit$solution
lambda_seq = taxa_fit$lambda_seq

solnNorm <- (abs(solution[, ncol(solution)]) - min(abs(solution[, ncol(solution)]))) / 
  (max(abs(solution[, ncol(solution)])) - min(abs(solution[, ncol(solution)])))

##Generate plots found in paper
##for full plot
plotlambdas <- plot(log(lambda_seq)[1:ncol(solution)], solution[1, ], type = "l",
                    xlim = c(log(lambda_seq[ncol(solution)]), log(lambda_seq[1])),
                    ylim = c(-1.2 * max(abs(solution[, ncol(solution)])),
                             1.2 * max(abs(solution[, ncol(solution)]))), xlab = "Log lambda",
                    ylab = "Param. Estimate", col = "black", main = "Full Solution Path",
                    cex.lab=1.5, cex.axis=1.5, cex.main=2, cex.sub=1.5, lwd = 1.5)

for (k in 2 : nrow(solution)){
  lines(log(lambda_seq)[1:ncol(solution)], solution[k, ], 
        col = ifelse(solution[k, 1] != 0, "black", "gray71"), lwd = 1.5)
}

abline(h = 0, lwd = 2, col = "black")

rect(xleft = 4.2, ybottom = -0.1, xright = 4.95, ytop = 0.1, border = "black", lwd = 2.5)

##for taxa zoomed in plot
plotlambdas <- plot(log(lambda_seq)[1:ncol(solution)], solution[1, ], type = "l",
                    xlim = c(4.2, 4.95),
                    ylim = c(-0.1, 0.1), xlab = "Log lambda",
                    ylab = "Param. Estimate", col = "black", main = "Taxa Solution Path",
                    cex.lab=1.5, cex.axis=1.5, cex.main=2, cex.sub=1.5, lwd = 1.5)




for (k in 2 : nrow(solution)){
  lines(log(lambda_seq)[1:ncol(solution)], solution[k, ], 
        col = ifelse(solution[k, 1] != 0, "gray71", "black"), lwd = 1.5)
}

abline(h = 0, lwd = 2, col = "black")

