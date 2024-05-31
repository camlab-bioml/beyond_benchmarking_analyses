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

trainMat <- model.matrix(~., trainMat)[,-1]
testMat <- model.matrix(~., testMat)[,-1]
testMat <- cbind(testMat, endBias5prime=rep(0, nrow(testMat)))

control <- trainControl(method="repeatedcv", number=10, repeats=3, index=folds)

print(head(testMat))
print(head(trainMat))

print(dim(testMat))
print(dim(trainMat))

silModel <- train(x=trainMat, y=as.numeric(trainScores$silCorImp), method="glmnet", trControl=control)
print(silModel)
silPreds <- predict(silModel, newdata=testMat)
silCoefs <- coef(silModel$finalModel, silModel$bestTune$lambda)
silTrainPreds <- predict(silModel, newdata=trainMat)

dbModel <- train(x=trainMat, y=as.numeric(trainScores$dbCorImp), method="glmnet", trControl=control)
print(dbModel)
dbPreds <- predict(dbModel, newdata=testMat)
dbCoefs <- coef(dbModel$finalModel, dbModel$bestTune$lambda)
dbTrainPreds <- predict(dbModel, newdata=trainMat)

chModel <- train(x=trainMat, y=as.numeric(trainScores$chCorImp), method="glmnet", trControl=control)
print(chModel)
chPreds <- predict(chModel, newdata=testMat)
chCoefs <- coef(chModel$finalModel, chModel$bestTune$lambda)
chTrainPreds <- predict(chModel, newdata=trainMat)

gseaModel <- train(x=trainMat, y=as.numeric(trainScores$gseaScaledImp), method="glmnet", trControl=control)
print(gseaModel)
gseaPreds <- predict(gseaModel, newdata=testMat)
gseaCoefs <- coef(gseaModel$finalModel, gseaModel$bestTune$lambda)
gseaTrainPreds <- predict(gseaModel, newdata=trainMat)


testData <- cbind(pipelines = testPipelines, name=testNames, testMat)
trainData <- cbind(pipelines=trainPipelines, name=trainNames, trainMat)


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

saveRDS(models, "/home/campbell/cfang/bb-rebuttal/results/models/linearRegression/linearRegressionPreds.RDS")
