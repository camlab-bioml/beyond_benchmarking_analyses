suppressPackageStartupMessages({
  library(tidyverse) # CRAN
  library(here) # CRAN
  library(DT) # CRAN
  library(here) # CRAN
  library(Matrix) # CRAN
  library(purrr)
  library(ggplot2)
  #library(Metrics)
  library(gridExtra)
  library(caret)
  library(matrixStats)
  library(randomForest)
  library(pcaMethods)
  library(stringr)
  library(stringi)
  #library(glmnet)
  #library(ggforce)
})
args <- commandArgs(trailingOnly = TRUE)

imputeMed <- function(x){
  x <- as.numeric(x)
  med <- median(x, na.rm=TRUE)
  x[is.na(as.numeric(x))] <- med
  return(as.numeric(x))
}


designMatPath <- args[[1]]
ariPath <- args[[2]]
numClustersPath <- args[[3]]
# designMatPath <- "uncorrected_designMatrix_scaled.RDS"
# ariPath <- "ari_unscaled.csv"
# numClustersPath <- "num_clusters.csv"

designMat <- readRDS(designMatPath)
ari <- read.csv(ariPath)
nclusts <- read.csv(numClustersPath)

avg_expr <- readRDS(args[[4]])
# avg_expr <- readRDS("avg_expr_pca_scores.RDS")
avg_expr_scores <- avg_expr[,1:20]
avg_expr_scores <- cbind(name=rownames(avg_expr_scores), avg_expr_scores)
rownames(avg_expr_scores) <- NULL
avg_expr_scores <- as_tibble(avg_expr_scores)

# remove duplicate predictors
designMat <- designMat[,-which(colnames(designMat) %in% c("sum.1", "detected.1", "total.1", "pct_counts_in_top_20_feature", "percent.top_50", "total_counts","total_features"))]
designMat <- merge(designMat, avg_expr_scores, by="name")
designMat <- designMat %>%
  rename(filt = V2.x)%>%
  rename(dims = V3.x)%>%
  rename(norm = V4.x)%>%
  rename(res = V5.x)

ari$X <- NULL
colnames(ari) <- unlist(lapply(colnames(ari), gsub, pattern="*\\.supervised\\.csv", replacement=""))
colnames(ari) <- unlist(lapply(colnames(ari), gsub, pattern="*\\.", replacement="-"))
ariExps <- colnames(ari)[which(grepl("^E-*",colnames(ari)))]

#flip sign of db
designMat$db <- designMat$db * -1
designMat$pipelines <- gsub(pattern="\"", replacement="", designMat$pipelines)

colnames(nclusts) <- unlist(lapply(colnames(nclusts), gsub, pattern="*\\.", replacement="-"))
nclusts$X <- gsub(pattern="\"", replacement="", nclusts$X)
nclusts$X <- gsub(pattern="\\;", replacement="\\.", nclusts$X)
nclusts$X <- gsub(pattern="\\=", replacement="\\.", nclusts$X)
colnames(nclusts)[which(colnames(nclusts)=="X")] <- "pipelines"

# Merge nclusts with designMat
nclusts <- as_tibble(nclusts)
nclusts <- pivot_longer(nclusts, cols=which(grepl("^E-*", colnames(nclusts))))
designClusts <- merge(designMat, nclusts, by=c("pipelines", "name"))
colnames(designClusts)[which(colnames(designClusts)=="value")] <- "nclusts"
designMat <- designClusts

# Read in GSEA scores
gsea <- read.csv(args[[5]])
# gsea <- read.csv("all_gsea_means.csv")
gsea$X <- NULL
gsea <- as_tibble(gsea) 
gseaNAs <- gsea %>% 
  group_by(name) %>%
  group_modify(~as.data.frame(sum(is.na(.$means))))%>%
  rename(numNA=2)

gseaImp <- gsea %>%
  group_by(name) %>% 
  group_modify(~ as.data.frame(cbind(pipelines=.$pipelines, means=imputeMed(.$means))))
gseaImp <- gseaImp %>%
  rename(gsea = means)
gseaImp$pipelines <- gsub(pattern="\"", replacement="", gseaImp$pipelines)
gseaImp$pipelines <- gsub(pattern="\\;", replacement="\\.", gseaImp$pipelines)
gseaImp$pipelines <- gsub(pattern="\\=", replacement="\\.", gseaImp$pipelines)

designMat <- merge(designMat, gseaImp, by=c("pipelines", "name")) 

numbers_only <- function(x) {
  sum(suppressWarnings(!is.na(as.numeric(as.character(x))))) == length(x)
}

designMat <- designMat %>% mutate_if(numbers_only, as.numeric)


silScoresCorrected <- designMat %>%
  group_by(name) %>%
  group_modify(~ as.data.frame(cbind(silCor=resid(loess(as.numeric(.$sil) ~ as.numeric(.$nclusts) ), na.action=na.exclude), pipelines=.$pipelines)))

dbScoresCorrected <- designMat %>%
  group_by(name) %>%
  group_modify(~ as.data.frame(cbind(dbCor=resid(loess(as.numeric(.$db) ~ as.numeric(.$nclusts) ), na.action=na.exclude),pipelines=.$pipelines)))

chScoresCorrected <- designMat %>%
  group_by(name) %>%
  group_modify(~ as.data.frame(cbind(chCor=resid(loess(as.numeric(.$ch) ~ as.numeric(.$nclusts) ), na.action=na.exclude),pipelines=.$pipelines)))


designMat <- merge(designMat, silScoresCorrected, by=c("name","pipelines"))
designMat <- merge(designMat, dbScoresCorrected, by=c("name","pipelines"))
designMat <- merge(designMat, chScoresCorrected, by=c("name","pipelines"))

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


silTrain <- designMatTrain %>%
  group_by(name) %>%
  group_modify(~ as.data.frame(cbind(pipelines=.$pipelines, silCorImp=imputeMed(.$silCor))))

silTest <- designMatTest %>%
  group_by(name) %>%
  group_modify(~ as.data.frame(cbind(pipelines=.$pipelines, silCorImp=imputeMed(.$silCor))))

dbTrain <- designMatTrain %>%
  group_by(name) %>%
  group_modify(~ as.data.frame(cbind(pipelines=.$pipelines, dbCorImp=imputeMed(.$dbCor))))

dbTest <- designMatTest %>%
  group_by(name) %>%
  group_modify(~ as.data.frame(cbind(pipelines=.$pipelines, dbCorImp=imputeMed(.$dbCor))))

chTrain <- designMatTrain %>%
  group_by(name) %>%
  group_modify(~ as.data.frame(cbind(pipelines=.$pipelines, chCorImp=imputeMed(.$chCor))))

chTest <- designMatTest %>%
  group_by(name) %>%
  group_modify(~ as.data.frame(cbind(pipelines=.$pipelines, chCorImp=imputeMed(.$chCor))))

mergeScores <- function(designMat, scores){
  merge(designMat, scores, by=c("pipelines", "name"))
}

createLRMat <- function(designMat){
  df <- designMat[,-which(colnames(designMat) %in% c("sil", "db","ch", "nclusts","silCor", "dbCor", "chCor","silCorImp", "dbCorImp", "chCorImp","gsea"))]
  df <- df %>% 
    mutate_if(numbers_only, as.numeric) %>%
    mutate_if(numbers_only, scale)
  return(df)
}


sildesignMatTrain <- mergeScores(designMatTrain, silTrain)
sildesignMatTest <- mergeScores(designMatTest, silTest)
silTrainData <- createLRMat(sildesignMatTrain)
silTestData <- createLRMat(sildesignMatTest)

dbdesignMatTrain <- mergeScores(designMatTrain, dbTrain)
dbdesignMatTest <- mergeScores(designMatTest, dbTest)
dbTrainData <- createLRMat(dbdesignMatTrain)
dbTestData <- createLRMat(dbdesignMatTest)

chdesignMatTrain <- mergeScores(designMatTrain, chTrain)
chdesignMatTest <- mergeScores(designMatTest, chTest)
chTrainData <- createLRMat(chdesignMatTrain)
chTestData <- createLRMat(chdesignMatTest)


# datasetR2s <- function(x){
#   datasetR2s <- x %>%
#     group_by(name) %>%
#     group_modify(~ as.data.frame(R2(as.numeric(.$preds), as.numeric(.$Actual))))
#   colnames(datasetR2s)[[2]] <- "R2"
#   print(datasetR2s[order(datasetR2s$R2),])
#   return(datasetR2s)
# }

# Random forest
tuneRFCV <- function(trainMat, testMat, trainScores, testScores){
  customRF <- list(type = "Regression",
                   library = "randomForest",
                   loop = NULL)
  
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
  tunegrid <- expand.grid(.mtry=c(20:30),.ntree=c(150))
  # train the model
  model <- train(x=data.matrix(trainMat), y=as.numeric(trainScores), method=customRF, trControl=control, tuneGrid=tunegrid)
  # summarize the model   
  print(model)
  return(model)
}

#silModel <- tuneRFCV(silTrainData, silTestData, sildesignMatTrain$silCorImp, sildesignMatTest$silCorImp)
#dbModel <- tuneRFCV(dbTrainData, dbTestData, dbdesignMatTrain$dbCorImp, dbdesignMatTest$dbCorImp)
chModel <- tuneRFCV(chTrainData, chTestData, chdesignMatTrain$chCorImp, chdesignMatTest$chCorImp)

#models <- list(sil=silModel, db=dbModel, ch=chModel)
saveRDS(chModel, "/home/campbell/cfang/automl_scrna/results/tuneRF/RFTunedCH.RDS")