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
# Read in data
args <- commandArgs(trailingOnly = TRUE)
designMatPath <- args[[1]]
ariPath <- args[[2]]

designMatPath <- "/home/campbell/cfang/automl_scrna/data/uncorrected_designMatrix_scaled.RDS"
ariPath <- "/home/campbell/cfang/automl_scrna/results/supervised_metrics/ari_unscaled.csv"

designMat <- readRDS(designMatPath)
ari <- read.csv(ariPath)

# think about removing collinear predictors
# remove duplicate predictors
designMat <- designMat[,-which(colnames(designMat) %in% c("sum.1", "detected.1", "total.1", "pct_counts_in_top_20_feature", "percent.top_50", "total_counts","total_features"))]

ari$X <- NULL
colnames(ari) <- unlist(lapply(colnames(ari), gsub, pattern="*\\.supervised\\.csv", replacement=""))
colnames(ari) <- unlist(lapply(colnames(ari), gsub, pattern="*\\.", replacement="-"))
ariExps <- colnames(ari)[which(grepl("^E-*",colnames(ari)))]

# subset to only resolution 0.5 pipelines
designMat <- designMat[which(str_detect(designMat$pipelines, "resolution.0.5")),]
ari <- ari[which(str_detect(ari$pipelines, "resolution.0.5")),]
ari_unscaled <- ari

#flip sign of db
designMat$db <- designMat$db * -1
designMat$pipelines <- gsub(pattern="\"", replacement="", designMat$pipelines)

# initialize train test sizes
PCT_TEST <- 0.3
PCT_TRAIN <- 1-PCT_TEST
N_EXPS <- length(unique(designMat$name))
N_TEST <- PCT_TEST*N_EXPS

# put all experiments with ARI into test set
ariTest <- designMat[which(designMat$name %in% ariExps),]
# sample remaining experiments to put in test set
remainingExpsNames <- unique(designMat$name[-which(designMat$name %in% ariExps)])
set.seed(0)
testExpNames <- sample(remainingExpsNames, size=N_TEST-length(ariExps), replace=FALSE)
testExps <- designMat[which(designMat$name %in% testExpNames),]
designMatTest <- rbind(testExps, ariTest)
# remove test set experiments from train set
designMatTrain <- designMat[-which(designMat$name %in% c(testExpNames, ariExps)),]

# pca the train set
trainMetrics <- designMatTrain[,which(colnames(designMatTrain) %in% c("sil", "db", "ch"))]
trainPCA <- ppca(as.matrix(trainMetrics))
trainScores <- trainPCA@scores[,1]
# transform the test set using pca model fitted with train set
testMetrics <- designMatTest[,which(colnames(designMatTest) %in% c("sil", "db", "ch"))]
testMetrics <- predict(trainPCA, newdata = as.matrix(testMetrics))
testScores <- testMetrics$scores[,1]
trainScoresMat <- as_tibble(cbind(designMatTrain, trainScores))
trainScoresMat <- select(trainScoresMat, -c("sil", "db", "ch", "V5"))

# control <- trainControl(method="cv", number=10)
# # train the model
# model <- train(trainScores ~., data=trainScoresMat, method="rf", trControl=control, tuneLength=5)
# # summarize the model
# print(model)


