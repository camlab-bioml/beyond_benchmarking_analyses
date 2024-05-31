suppressPackageStartupMessages({
  library(tidyverse) # CRAN
  library(here) # CRAN
  library(DT) # CRAN
  library(here) # CRAN
  library(Matrix) # CRAN
  library(purrr)
  library(Metrics)
  library(caret)
  library(matrixStats)
  library(glmnet)
})

args <- commandArgs(trailingOnly = TRUE)
models <- readRDS(args[[1]])
test <- readRDS(args[[2]])

testMat <- test$data

testNames <- testMat$name
testPipelines <- testMat$pipelines
testMat$name <- NULL
testMat$pipelines <- NULL

form <- as.formula(~ .*.)

testMatModel <- model.matrix(form, testMat)[,-1]

filts <- colnames(testMatModel)[grepl("^filt", colnames(testMatModel))][1:2]
norms <- colnames(testMatModel)[grepl("^norm", colnames(testMatModel))][1:2]

pipelineFeats <- c(filts, norms, "dims", "res")
datasetFeats <- setdiff(colnames(testMat), c("filt", "norm", "dims", "res"))
intTerms <- c(t(outer(pipelineFeats, datasetFeats, paste)))
intTerms <- gsub(" ", ":", intTerms)

intModelFeats <- c(pipelineFeats, datasetFeats, intTerms)

data <- testMatModel[,intModelFeats]


sil <- models$sil$model
print(colnames(data))

silPreds <- predict(sil, data, s=sil$bestTune$lambda)



db <- models$db$model
dbPreds <- predict(db, data)

ch <- models$ch$model
chPreds <- predict(ch, data)

gsea <- models$gsea$model
gseaPreds <- predict(gsea, data)

print(head(silPreds))

saveRDS(list(sil=silPreds, db=dbPreds, ch=chPreds, gsea=gseaPreds), "/home/campbell/cfang/bb-rebuttal/results/models/linearRegression/linear_regression_large_samples_preds.RDS")
