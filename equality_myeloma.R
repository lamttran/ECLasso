library(ADMMRcppArma)
library(ECLasso)
library(glmnet)

##Choose the number of high-expression genes to examine
##In our paper, we considered 1000, 2000, and 5000
num_topgenes = 1000


##Random seed to break ties in event times
set.seed(2021)

##Load in the RData file
load("~/MM_data.RData")
gene_expression_final <- t(genemat)

##removing info from 3 patients whose event times were exactly 0
gene_expression_final <- gene_expression_final[, -c(69, 126, 469)] 
gene_expression_final <- t(gene_expression_final)
censor_final = gseClinic$EFS.CENSOR..1.event..JUN2008[-c(69, 126, 469)]
times_final = gseClinic$EFS.TIME.JUN2008[-c(69, 126, 469)]

##Adding some noise to expression levels to break ties
times_final = 1 + times_final + rnorm(length(times_final), 0, 0.01)

##Setting up event time and censoring vectors
ord = order(times_final)
gene_expression_final = gene_expression_final[ord, ]
censor_final = censor_final[ord]
times_final = times_final[ord]

##removing genes with low measurement quality
genes_remove = c("REEP5", "SLC7A7", "SH3BP5", "LBR")
gene_expression_final =  gene_expression_final[,-which(colnames(gene_expression_final) %in% genes_remove)]

##preselection by expression levels and center
genes_to_keep = sort(order(colMeans(gene_expression_final), decreasing = T)[1:num_topgenes])
gene_expression_final <- gene_expression_final[, genes_to_keep]
n = nrow(gene_expression_final)

for(i in 1:ncol(gene_expression_final)){
  mean_center = mean(gene_expression_final[, i])
  gene_expression_final[, i] = gene_expression_final[, i] - mean_center
  rm(mean_center)
}

##The constrained fit, a simple sum-to-zero constraint
tolerance = 1e-4
myeloma_fit = ECLasso(x = gene_expression_final, y = times_final, intercept = F, 
                      censor = censor_final, family = "cox", tolerance = tolerance, 
                      method = "combined_subset", trace = T)

##Plotting
solution = myeloma_fit$solution
lambda_seq = myeloma_fit$lambda_seq

solnNorm <- (abs(solution[, ncol(solution)]) - min(abs(solution[, ncol(solution)]))) / 
  (max(abs(solution[, ncol(solution)])) - min(abs(solution[, ncol(solution)])))

plotlambdas <- plot(log(lambda_seq)[1:ncol(solution)], solution[1, ], type = "l",
                    xlim = c(log(lambda_seq[ncol(solution)]), log(lambda_seq[1])),
                    ylim = c(-1.2 * max(abs(solution[, ncol(solution)])),
                             1.2 * max(abs(solution[, ncol(solution)]))), xlab = "Log lambda",
                    ylab = "Param. Estimate", col = "black", main = "Myeloma Constrained Fit",
                    cex.lab=1.5, cex.axis=1.5, cex.main=2, cex.sub=1.5)

for (k in 2 : nrow(solution)){
  lines(log(lambda_seq)[1:ncol(solution)], solution[k, ], col = "black", lwd = 1.5)
}

abline(h = 0, lwd = 2)