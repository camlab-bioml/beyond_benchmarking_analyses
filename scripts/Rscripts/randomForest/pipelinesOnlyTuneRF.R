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
args <- commandArgs(trailingOnly=TRUE)

designMatTrain <- readRDS(args[[1]])
metric <- args[[2]]

# designMatTrain <- readRDS("corrected_RF_trainMatFactor_unscaled.RDS")
# metric <- "sil"
pipelineFeats <- c("filt", "norm", "res", "dims")

# Random forest
tuneRFCV <- function(trainMat, trainScores){
  customRF <- list(type = "Regression",
                   library = "randomForest",
                   loop = NULL)
  
  #trainMat$name <- NULL
  trainMat$pipelines <- NULL
  
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
  
  # Create dataset-aware folds
  set.seed(1)
  folds <- groupKFold(trainMat$name, k=10)
  trainMat$name <- NULL
  trainMat <- trainMat[,pipelineFeats]
  print(head(trainMat))
  control <- trainControl(method="repeatedcv", number=10, repeats=3, index=folds)
  tunegrid <- expand.grid(.mtry=c(1:ncol(trainMat)),.ntree=c(100, 300, 500))
  # train the model
  model <- train(x=trainMat, y=as.numeric(trainScores), method=customRF, trControl=control, tuneGrid=tunegrid)
  # summarize the model   
  print(model)
  return(model)
}

trainScores <- designMatTrain$labels %>%
  select(contains(metric))%>%
  select(contains("Imp"))

head(trainScores[,1])
trainMat <- designMatTrain$data
tunedModel <- tuneRFCV(trainMat, trainScores[,1])

modelPath <- sprintf("/home/campbell/cfang/automl_scrna/results/tuneRF/pipelinesOnly/RF%s_pipelinesOnly.RDS", metric)
saveRDS(tunedModel, modelPath)