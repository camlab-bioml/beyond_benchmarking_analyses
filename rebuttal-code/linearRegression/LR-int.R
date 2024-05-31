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

args <- commandArgs(trailingOnly=TRUE)

train <- readRDS(args[[1]])
test <- readRDS(args[[2]])

trainMat <- train$data
testMat <- test$data
trainScores <- train$labels
testScores <- test$labels

set.seed(1)
folds <- groupKFold(trainMat$name, k=10)

trainNames <- trainMat$name
trainPipelines <- trainMat$pipelines
trainMat$name <- NULL
trainMat$pipelines <- NULL

testNames <- testMat$name
testPipelines <- testMat$pipelines
testMat$name <- NULL
testMat$pipelines <- NULL

testMat <- cbind(testMat, endBias5prime=rep(0, nrow(testMat)))

form <- as.formula(~ .*.)

trainMatModel <- model.matrix(form, trainMat)[,-1]
testMatModel <- model.matrix(form, testMat)[,-1]

filts <- colnames(trainMatModel)[grepl("^filt", colnames(trainMatModel))][1:2]
norms <- colnames(trainMatModel)[grepl("^norm", colnames(trainMatModel))][1:2]

endBiases <- colnames(trainMatModel)[grepl("^endBias", colnames(trainMatModel))][1:2]

pipelineFeats <- c(filts, norms, "dims", "res")
datasetFeats <- setdiff(colnames(trainMat), c("filt", "norm", "dims", "res", "FromBlood", "endBias", "isFACS"))
datasetFeats <- c(datasetFeats, "FromBloodyes", "isFACSyes", endBiases)
intTerms <- c(t(outer(pipelineFeats, datasetFeats, paste)))
intTerms <- gsub(" ", ":", intTerms)

intModelFeats <- c(pipelineFeats, datasetFeats, intTerms)
# print(head(trainMatModel))
# print(intModelFeats)




trainMatModelInt <- trainMatModel[,intModelFeats]
testMatModelInt <- testMatModel[,intModelFeats]

control <- trainControl(method="repeatedcv", number=10, repeats=3, index=folds)

print(head(trainMatModel))
print(head(testMatModel))

print(dim(trainMatModel))
print(dim(testMatModel))

colnames(testMatModel)[which(colnames(testMatModel)=="isFACSyes:endBias5prime")] <- "endBias5prime:isFACSyes"

print(setdiff(colnames(trainMatModel), colnames(testMatModel)))
print(setdiff(colnames(testMatModel), colnames(trainMatModel)))

silModel <- train(x=trainMatModelInt, y=as.numeric(trainScores$silCorImp), method="glmnet", trControl=control)
print(silModel)
silPreds <- predict(silModel$finalModel, newx=testMatModelInt, s=silModel$bestTune$lambda)
silCoefs <- coef(silModel$finalModel, silModel$bestTune$lambda)
silTrainPreds <- predict(silModel$finalModel, newx=trainMatModelInt, s=silModel$bestTune$lambda)

dbModel <- train(x=trainMatModelInt, y=as.numeric(trainScores$dbCorImp), method="glmnet", trControl=control)
print(dbModel)
dbPreds <- predict(dbModel$finalModel, newx=testMatModelInt, s=dbModel$bestTune$lambda)
dbCoefs <- coef(dbModel$finalModel, dbModel$bestTune$lambda)
dbTrainPreds <- predict(dbModel$finalModel, newx=trainMatModelInt, s=dbModel$bestTune$lambda)

chModel <- train(x=trainMatModelInt, y=as.numeric(trainScores$chCorImp), method="glmnet", trControl=control)
print(chModel)
chPreds <- predict(chModel$finalModel, newx=testMatModelInt, s=chModel$bestTune$lambda)
chCoefs <- coef(chModel$finalModel, chModel$bestTune$lambda)
chTrainPreds <- predict(chModel$finalModel, newx=trainMatModelInt, s=chModel$bestTune$lambda)

gseaModel <- train(x=trainMatModelInt, y=as.numeric(trainScores$gseaScaledImp), method="glmnet", trControl=control)
print(gseaModel)
gseaPreds <- predict(gseaModel$finalModel, newx=testMatModelInt, s=gseaModel$bestTune$lambda)
gseaCoefs <- coef(gseaModel$finalModel, gseaModel$bestTune$lambda)
gseaTrainPreds <- predict(gseaModel$finalModel, newx=trainMatModelInt, s=gseaModel$bestTune$lambda)


testData <- cbind(pipelines = testPipelines, name=testNames, testMatModelInt)
trainData <- cbind(pipelines=trainPipelines, name=trainNames, trainMatModelInt)


silTrainData <- cbind(trainData, silCorImp = as.numeric(trainScores$silCorImp))
dbTrainData <- cbind(trainData, dbCorImp = as.numeric(trainScores$dbCorImp))
chTrainData <- cbind(trainData, chCorImp = as.numeric(trainScores$chCorImp))
gseaTrainData <- cbind(trainData, gseaScaledImp = as.numeric(trainScores$gseaScaledImp))

silTestData <- cbind(testData, silCorImp = as.numeric(testScores$silCorImp))
dbTestData <- cbind(testData, dbCorImp = as.numeric(testScores$dbCorImp))
chTestData <- cbind(testData, chCorImp = as.numeric(testScores$chCorImp))
gseaTestData <- cbind(testData, gseaScaledImp = as.numeric(testScores$gseaScaledImp))

sil <- list(model=silModel, testPreds=silPreds, coefs=silCoefs, testData=silTestData, trainData=silTrainData, trainPreds=silTrainPreds)
db <- list(model=dbModel, testPreds=dbPreds, coefs=dbCoefs, testData=dbTestData, trainData=dbTrainData, trainPreds=dbTrainPreds)
ch <- list(model=chModel, testPreds=chPreds, coefs=chCoefs, testData=chTestData, trainData=chTrainData, trainPreds=chTrainPreds)
gsea <- list(model=gseaModel, testPreds=gseaPreds, coefs=gseaCoefs, testData=gseaTestData, trainData=gseaTrainData, trainPreds=gseaTrainPreds)

models = list(sil=sil, db=db, ch=ch, gsea=gsea)

saveRDS(models, "/home/campbell/cfang/bb-rebuttal/results/models/linearRegression/allTermsInteractions/LR_int_allTerms.RDS")
