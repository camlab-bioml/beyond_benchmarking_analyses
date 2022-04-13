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

createRFMat <- function(designMat, paramsFactor=FALSE){
  res <- designMat$res
  dims <- designMat$dims
  df <- designMat[,-which(colnames(designMat) %in% c("sil", "db","ch", "nclusts","silCor", "dbCor", "chCor","silCorImp", "dbCorImp", "chCorImp","gseaScaledImp"))]
  df <- df %>% 
    mutate_if(numbers_only, as.numeric) %>%
    mutate_if(numbers_only, scale)
  
  df$filt <- as.factor(df$filt)
  df$norm <- as.factor(df$norm)
  
  if (paramsFactor){
    df$res <- as.factor(res)
    df$dims <- as.factor(dims)
  }else{
    df$res <- as.numeric(res)
    df$dims <- as.numeric(dims)
  }
  labels <- designMat[,which(colnames(designMat) %in% c("silCorImp", "dbCorImp", "chCorImp","gseaScaledImp"))]
  res <- list(data=df, labels=labels)
  return(res)
}

args <- commandArgs(trailingOnly=TRUE)

designMat <- readRDS(args[[1]])
trainTestSplit <- readRDS(args[[2]])

# train test split
designMatTrain <- designMat[which(designMat$name %in% trainTestSplit$train),]
designMatTest <- designMat[which(designMat$name %in% trainTestSplit$test),]

trainMatFactor <- createRFMat(designMatTrain, paramsFactor=TRUE)
trainMatNum <- createRFMat(designMatTrain, paramsFactor=FALSE)

testMatFactor <- createRFMat(designMatTest, paramsFactor=TRUE)
testMatNum <- createRFMat(designMatTest, paramsFactor=FALSE)

saveRDS(trainMatFactor, "/home/campbell/cfang/automl_scrna/data/corrected_RF_trainMatFactor_unscaled.RDS")
saveRDS(testMatFactor, "/home/campbell/cfang/automl_scrna/data/corrected_RF_testMatFactor_unscaled.RDS")

saveRDS(trainMatNum, "/home/campbell/cfang/automl_scrna/data/corrected_RF_trainMatNumeric_unscaled.RDS")
saveRDS(testMatNum, "/home/campbell/cfang/automl_scrna/data/corrected_RF_testMatNumeric_unscaled.RDS")


