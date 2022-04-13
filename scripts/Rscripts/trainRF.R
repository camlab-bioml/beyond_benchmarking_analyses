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

numbers_only <- function(x) {
  sum(suppressWarnings(!is.na(as.numeric(as.character(x))))) == length(x)
}

createRFMat <- function(designMat){
  res <- designMat$res
  dims <- designMat$dims
  df <- designMat[,-which(colnames(designMat) %in% c("sil", "db","ch", "nclusts","silCor", "dbCor", "chCor","silCorImp", "dbCorImp", "chCorImp","gseaScaledImp"))]
  df <- df %>% 
    mutate_if(numbers_only, as.numeric) %>%
    mutate_if(numbers_only, scale)
  df$filt <- as.factor(df$filt)
  df$norm <- as.factor(df$norm)
  df$res <- as.numeric(res)
  df$dims <- as.numeric(dims)
  return(df)
}

args <- commandArgs(trailingOnly=TRUE)

designMat <- readRDS(args[[1]])
ari <- readRDS(args[[2]])
trainTestSplit <- readRDS(args[[3]])
exprScores <- readRDS(args[[4]])
print(head(exprScores))

avg_expr_scores <- rbind(exprScores$trainExpr, exprScores$testExpr)
avg_expr_scores <- avg_expr_scores[,1:20]
avg_expr_scores <- cbind(name=rownames(avg_expr_scores), avg_expr_scores)
rownames(avg_expr_scores) <- NULL
avg_expr_scores <- as_tibble(avg_expr_scores)
print(head(avg_expr_scores))

designMat <- merge(designMat, avg_expr_scores, by="name")

# train test split
designMatTrain <- designMat[which(designMat$name %in% trainTestSplit$train),]
designMatTest <- designMat[which(designMat$name %in% trainTestSplit$test),]

trainMat <- createRFMat(designMatTrain)
testMat <- createRFMat(designMatTest)

trainRF <- function(trainData, testData, trainScores, testScores, mtry, ntree=150){
  trainName <- trainData$name
  trainPipelines <- trainData$pipelines
  trainData$name <- NULL
  trainData$pipelines <- NULL
  
  testName <- testData$name
  testPipelines <- testData$pipelines
  testData$name <- NULL
  testData$pipelines <- NULL
  
  rf <- randomForest(data.matrix(trainData), as.numeric(trainScores), mtry=mtry, ntree=ntree)
  preds <- predict(rf, newdata=data.matrix(testData))
  
  testRsq <- R2(preds, as.numeric(testScores))
  trainPreds <- predict(rf, newdata=data.matrix(trainData))
  trainRsq <- R2(trainPreds, as.numeric(trainScores))
  
  trainData <- cbind(pipelines=trainPipelines, name=trainName, trainData, Actual=trainScores)
  testData <- cbind(pipelines=testPipelines, name=testName, testData, Actual=testScores)
  res <- list(Rsqs=c(testRsq, trainRsq), model=rf, testPreds=preds, trainPreds=trainPreds, testData=testData, trainData=trainData, mtry=mtry, ntree=ntree)
  return(res)
}

print("sil")
silRF <- trainRF(trainMat, testMat, designMatTrain$silCorImp, designMatTest$silCorImp, mtry=40)
print("db")
dbRF <- trainRF(trainMat, testMat, designMatTrain$dbCorImp, designMatTest$dbCorImp, mtry=40)
print("ch")
chRF <- trainRF(trainMat, testMat, designMatTrain$chCorImp, designMatTest$chCorImp, mtry=40)
print("gsea")
gseaRF <- trainRF(trainMat, testMat, designMatTrain$gseaScaledImp, designMatTest$gseaScaledImp, mtry=40)


models <- list(sil=silRF, db=dbRF, ch=chRF, gsea=gseaRF)
saveRDS(models, "/home/campbell/cfang/automl_scrna/results/models/randomForest/RFresults.RDS")

