suppressPackageStartupMessages({
  library(tidyverse) # CRAN
  library(here) # CRAN
  library(DT) # CRAN
  library(here) # CRAN
  library(Matrix) # CRAN
  library(purrr)
  library(ggplot2)
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
  df$dim <- as.numeric(dims)
  return(df)
}

args <- commandArgs(trailingOnly=TRUE)

designMat <- readRDS(args[[1]])
ari <- readRDS(args[[2]])
trainTestSplit <- readRDS(args[[3]])
exprScores <- readRDS(args[[4]])

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

# Random forest
tuneRFCV <- function(trainMat, testMat, trainScores, testScores){
  customRF <- list(type = "Regression",
                   library = "randomForest",
                   loop = NULL)
  
  trainMat$name <- NULL
  trainMat$pipelines <- NULL
  testMat$name <- NULL
  testMat$pipelines <- NULL
  
  customRF$parameters <- data.frame(parameter = c("mtry", "ntree"),
                                    class = rep("numeric", 2),
                                    label = c("mtry", "ntree"))
  
  customRF$grid <- function(x, y, len = NULL, search = "grid") {}
  
  customRF$fit <- function(x, y, wts, param, lev, last, weights, classProbs) {
    randomForest(x, y,
                 mtry = param$mtry,
                 ntree=param$ntree)
  }
  
  #Predict label
  customRF$predict <- function(modelFit, newdata, preProc = NULL, submodels = NULL)
    predict(modelFit, newdata)
  
  #Predict prob
  customRF$prob <- function(modelFit, newdata, preProc = NULL, submodels = NULL)
    predict(modelFit, newdata, type = "response")
  
  customRF$sort <- function(x) x[order(x[,1]),]
  customRF$levels <- function(x) x$classes
  
  control <- trainControl(method="repeatedcv", number=10, repeats=3)
  tunegrid <- expand.grid(.mtry=c(20:40),.ntree=c(150))
  # train the model
  model <- train(x=data.matrix(trainMat), y=as.numeric(trainScores), method=customRF, trControl=control, tuneGrid=tunegrid)
  # summarize the model   
  print(model)
  return(model)
}

gseaModel <- tuneRFCV(trainMat, testMat, designMatTrain$gseaScaledImp, designMatTest$gseaScaledImp)
saveRDS(gseaModel, "/home/campbell/cfang/automl_scrna/results/tuneRF/RFTunedgsea.RDS")

