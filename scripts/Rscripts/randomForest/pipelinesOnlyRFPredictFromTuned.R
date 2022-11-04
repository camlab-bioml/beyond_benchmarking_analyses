suppressPackageStartupMessages({
  library(tidyverse) # CRAN
  library(here) # CRAN
  library(DT) # CRAN
  library(here) # CRAN
  library(Matrix) # CRAN
  library(purrr)
  library(ggplot2)
  library(Metrics)
  library(gridExtra)
  library(caret)
  library(matrixStats)
  library(randomForest)
  library(pcaMethods)
  library(stringr)
  library(stringi)
})


args <- commandArgs(trailingOnly=TRUE)

designMatTrain <- readRDS(args[[1]])
designMatTest <- readRDS(args[[2]])
# Read in trained models
silModel <- readRDS(args[[3]])$finalModel
dbModel <- readRDS(args[[4]])$finalModel
chModel <- readRDS(args[[5]])$finalModel
gseaModel <- readRDS(args[[6]])$finalModel
#modelType <- args[[7]]

print(args)
pipelineFeats <- c("filt", "norm", "res", "dims")

predRF <- function(trainData, testData, trainScores, testScores, tunedRF){
  # Save name and pipeline values but don't include in model
  trainName <- trainData$name
  trainPipelines <- trainData$pipelines
  trainData$name <- NULL
  trainData$pipelines <- NULL
  
  testName <- testData$name
  testPipelines <- testData$pipelines
  testData$name <- NULL
  testData$pipelines <- NULL
  
  # Predict on train and test sets
  print(names(testData))
  testData <- testData[,pipelineFeats]
  print(head(testData))
  preds <- predict(tunedRF, newdata=testData)
  testRsq <- R2(preds, as.numeric(testScores))
  
  trainData <- trainData[,pipelineFeats]
  print(head(trainData))
  trainPreds <- predict(tunedRF, newdata=trainData)
  trainRsq <- R2(trainPreds, as.numeric(trainScores))
  
  # Add pipelines and names back to data
  trainData <- cbind(pipelines=trainPipelines, name=trainName, trainData, Actual=trainScores)
  testData <- cbind(pipelines=testPipelines, name=testName, testData, Actual=testScores)
  print(tunedRF)
  res <- list(Rsqs=c(testRsq, trainRsq), model=tunedRF, testPreds=preds, trainPreds=trainPreds, testData=testData, trainData=trainData, mtry=tunedRF$tuneValue$mtry, ntree=tunedRF$tuneValue$ntree)
  return(res)
}




trainMat <- designMatTrain$data
testMat <- designMatTest$data

print("sil")
silRF <- predRF(trainMat, testMat, designMatTrain$labels$silCorImp, designMatTest$labels$silCorImp, silModel)
print("db")
dbRF <- predRF(trainMat, testMat, designMatTrain$labels$dbCorImp, designMatTest$labels$dbCorImp, dbModel)
print("ch")
chRF <- predRF(trainMat, testMat, designMatTrain$labels$chCorImp, designMatTest$labels$chCorImp, chModel)
print("gsea")
gseaRF <- predRF(trainMat, testMat, designMatTrain$labels$gseaScaledImp, designMatTest$labels$gseaScaledImp, gseaModel)


models <- list(sil=silRF, db=dbRF, ch=chRF, gsea=gseaRF)
saveRDS(models, sprintf("/home/campbell/cfang/automl_scrna/results/models/randomForest/pipelinesOnlyRFPreds.RDS"))

