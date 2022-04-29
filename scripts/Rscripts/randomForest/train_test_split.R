suppressPackageStartupMessages({
  library(tidyverse) # CRAN
  library(here) # CRAN
  library(DT) # CRAN
  library(here) # CRAN
  library(Matrix) # CRAN
  library(stats)
})
args <- commandArgs(trailingOnly=TRUE)
designMat <- readRDS(args[[1]])
ari <- read.csv(args[[2]])

ari$X <- NULL
colnames(ari) <- unlist(lapply(colnames(ari), gsub, pattern="*\\.supervised\\.csv", replacement=""))
colnames(ari) <- unlist(lapply(colnames(ari), gsub, pattern="*\\.", replacement="-"))
ariExps <- colnames(ari)[which(grepl("^E-*",colnames(ari)))]

# initialize train test sizes
PCT_TEST <- 0.3
PCT_TRAIN <- 1-PCT_TEST
N_EXPS <- length(unique(designMat$name))
N_TEST <- floor(PCT_TEST*N_EXPS)

# put all experiments with ARI into test set
ariTest <- designMat[which(designMat$name %in% ariExps),]
# sample remaining experiments to put in test set
remainingExpsNames <- unique(designMat$name[-which(designMat$name %in% ariExps)])
set.seed(1)
testExpNames <- sample(remainingExpsNames, size=N_TEST-length(ariExps), replace=FALSE)
testExps <- designMat[which(designMat$name %in% testExpNames),]
designMatTest <- rbind(testExps, ariTest)

# remove test set experiments from train set
#designMatTrain <- designMat[-which(designMat$name %in% c(testExpNames, ariExps)),]
designMatTrain <- designMat[which(designMat$name %in% setdiff(designMat$name, c(testExpNames, ariExps))),]
print(unique(designMatTest$name))
print(unique(designMatTrain$name))

split <- list(test=unique(designMatTest$name), train=unique(designMatTrain$name))
saveRDS(split, "/home/campbell/cfang/automl_scrna/data/train_test_split.RDS")