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

scale_test <- function(testMat, trainMeans, trainVars){
  for (i in 1:ncol(testMat)){
    mean <- trainMeans[which(names(trainMeans) == colnames(testMat)[[i]])]
    var <- trainVars[which(names(trainVars) == colnames(testMat)[[i]])]
    sd <- sqrt(var)
    col <- (as.numeric(testMat[,i]) - mean)/sd
    testMat[,i] <- col
  }
  return(testMat)
}

createRFMat <- function(designMat, testSet=FALSE, trainMeans=NULL, trainVars=NULL){
  df <- designMat[,-which(colnames(designMat) %in% c("sil", "db","ch", "nclusts","silCor", "dbCor", "chCor","silCorImp", "dbCorImp", "chCorImp","gseaScaledImp"))]
  if (testSet) {
    df <- df %>%
      select(-c("name", "pipelines", "filt", "norm"))

    df <- scale_test(df, trainMeans, trainVars)

    df$name <- designMat$name
    df$pipelines <- designMat$pipelines
    df$filt <- designMat$filt
    df$norm <- designMat$norm
  }else{
    trainMeans <- df %>%
      select(-c("name", "pipelines", "filt", "norm"))%>%
      mutate_if(numbers_only, as.numeric) %>%
      colMeans()

    print(colnames(df))
    trainVars <- df %>%
      select(-c("name", "pipelines", "filt", "norm"))%>%
      mutate_if(numbers_only, as.numeric) %>%
      as.matrix()%>%
      colVars()

    df <- df %>%
      mutate_if(numbers_only, as.numeric) %>%
      mutate_if(numbers_only, scale) # use mean and variance from train on the test set


    names(trainVars) <- names(trainMeans)
    labels <- designMat[,which(colnames(designMat) %in% c("silCorImp", "dbCorImp", "chCorImp","gseaScaledImp"))]
    res <- list(data=df, labels=labels, trainMeans=trainMeans, trainVars=trainVars) # change to numeric
    return(res)
    }

  labels <- designMat[,which(colnames(designMat) %in% c("silCorImp", "dbCorImp", "chCorImp","gseaScaledImp"))]
  res <- list(data=df, labels=labels) # change to numeric
  return(res)
}

args <- commandArgs(trailingOnly=TRUE)

designMat <- readRDS(args[[2]])
designMatLarge <- readRDS(args[[1]])
trainTestSplit <- readRDS(args[[3]])


# train test split
designMatTrain <- designMat[which(designMat$name %in% trainTestSplit$train),]

top50 <- designMatLarge$`percent.top_50`
designMatLarge$pct_counts_in_top_50_features <- top50
designMatLarge$`percent.top_50` <- NULL

top20 <- designMatLarge$pct_counts_in_top_20_features
designMatLarge$pct_counts_top_20_features <- top20
designMatLarge$pct_counts_in_top_20_features <- NULL

trainMatNum <- createRFMat(designMatTrain, testSet = FALSE)
testMatNum <- createRFMat(designMatLarge, testSet=TRUE, trainMeans = trainMatNum$trainMeans, trainVars = trainMatNum$trainVars)

print(colnames(trainMatNum$data))
print(colnames(testMatNum$data))

print(setdiff(colnames(trainMatNum$data), colnames(testMatNum$data)))
testMatNum$data <- testMatNum$data[,colnames(trainMatNum$data)]

print(testMatNum)



# in large:
# pct_counts_in_top_20_features
# percent.top_50

#saveRDS(trainMatNum, "/home/campbell/cfang/automl_scrna/data/corrected_RF_trainMatNumeric_unscaled.RDS")
saveRDS(testMatNum, "/home/campbell/cfang/bb-rebuttal/data/large_samples_corrected_RF_testMatNumeric_unscaled.RDS")


