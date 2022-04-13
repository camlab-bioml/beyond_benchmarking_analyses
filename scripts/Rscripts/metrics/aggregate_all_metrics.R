suppressPackageStartupMessages({
  library(tidyverse) # CRAN
  library(here) # CRAN
  library(DT) # CRAN
  library(here) # CRAN
  library(Matrix) # CRAN
  library(scran) # BioConductor
  library(purrr)
})
expsMatrices <- commandArgs(trailingOnly = TRUE)
sil <- c()
db <- c()
ch <- c()
for (i in 1:length(expsMatrices)){
  metrics <- read.csv(expsMatrices[[i]])
  colnames(metrics)[[1]] <- "pipelines"
  expName <- strsplit(expsMatrices[[i]], "-metrics.csv")[[1]][[1]]
  expName <- strsplit(expName, "metrics_matrices/")[[1]][[2]]
  print(expName)
  print(dim(metrics))
  names(metrics)[[1]] <- "index"
  print(head(metrics))
  expNames <- rep(expName, nrow(metrics))
  sil[[i]] <- data.frame("pipelines"=metrics$pipelines, "sil"=metrics$sil, "exp"=expNames)
  db[[i]] <- data.frame("pipelines"=metrics$pipelines, "db"=metrics$db, "exp"=expNames)
  ch[[i]] <- data.frame("pipelines"=metrics$pipelines, "ch"=metrics$ch, "exp"=expNames)
}
# do.call(rbind(sil)), then make it a tibble and pivot
silMat <- as_tibble(do.call(rbind, sil))
dbMat <- as_tibble(do.call(rbind, db))
chMat <- as_tibble(do.call(rbind, ch))

silMat <- pivot_wider(silMat, names_from = "exp", values_from="sil")
dbMat <- pivot_wider(dbMat, names_from = "exp", values_from="db")
chMat <- pivot_wider(chMat, names_from = "exp", values_from="ch")

silName <- "/home/campbell/cfang/automl_scrna/results/metrics/sil_unscaled.csv"
dbName <- "/home/campbell/cfang/automl_scrna/results/metrics/db_unscaled.csv"
chName <- "/home/campbell/cfang/automl_scrna/results/metrics/ch_unscaled.csv"

print(head(silMat))
print(head(dbMat))
print(head(chMat))

write.csv(silMat, silName)
write.csv(dbMat, dbName)
write.csv(chMat, chName)

sils <- rep("sil", nrow(silMat))
dbs <- rep("db", nrow(dbMat))
chs <- rep("ch", nrow(chMat))

silMat <- cbind(silMat, score=sils)
dbMat <- cbind(dbMat, score=dbs)
chMat <- cbind(chMat, score=chs)

# Remove 0.01 resolution pipelines
silMat <- silMat[which(!str_detect(silMat$pipelines, "resolution.0.01")),]
dbMat <- dbMat[which(!str_detect(dbMat$pipelines, "resolution.0.01")),]
chMat <- chMat[which(!str_detect(chMat$pipelines, "resolution.0.01")),]

#Scale each dataset
silMat[,sapply(silMat, is.numeric)] <- scale(silMat[,sapply(silMat, is.numeric)])
dbMat[,sapply(dbMat, is.numeric)] <- scale(dbMat[,sapply(dbMat, is.numeric)])
chMat[,sapply(chMat, is.numeric)] <- scale(chMat[,sapply(chMat, is.numeric)])

# Now pivot metrics to 288*86 x 3
allMetrics <- as_tibble(rbind(silMat, dbMat, chMat))
allMetrics <- pivot_longer(allMetrics, cols=starts_with("E-"))
allMetrics <- pivot_wider(allMetrics, names_from="score")

print(head(allMetrics))
saveRDS(allMetrics, "/home/campbell/cfang/automl_scrna/results/metrics/all_metrics.RDS")

# # Create design matrix
# pipelines <- unique(allMetrics$pipelines)
# filt <- lapply(pipelines, function(pipeline){
#   stringency <- strsplit(pipeline, "filt.")[[1]][[3]]
#   stringency <- strsplit(stringency, ".norm")[[1]][[1]]
#   return(stringency)
# })
# filt <- unlist(filt)
# dims  <- lapply(pipelines, function(pipeline){
#   dim <- strsplit(pipeline, "dims.")[[1]][[2]]
#   dim <- strsplit(dim, ".k")[[1]][[1]]
#   return(dim)
# })
# dims <- unlist(dims)
# norms  <- lapply(pipelines, function(pipeline){
#   norm <- strsplit(pipeline, "norm.")[[1]][[3]]
#   norm <- strsplit(norm, ".sel")[[1]][[1]]
#   return(norm)
# })
# norms <- unlist(norms)
# res <- lapply(pipelines, function(pipeline){
#   clust.res <- strsplit(pipeline, "resolution.")[[1]][[2]]
#   clust.res <- strsplit(clust.res, ".min")[[1]][[1]]
#   return(clust.res)
# })
# res <- unlist(res)
# design <- as.data.frame(cbind(pipelines=pipelines, filt=filt, dims=dims, norms=norms, res=res))
# 
# expData <- read.csv("regression_design.csv")
# expData <- expData[,!(colnames(expData) %in% c("phenoid"))]
# expData <- as_tibble(expData)
# 
# 
# designMatrix <- lapply(pipelines, function(pipeline){
#   pipeline <- matrix(rep(design[which(design$pipelines == pipeline),], each = nrow(expData)), nrow=nrow(expData))
#   pipeline <- cbind(as.data.frame(pipeline), expData)
#   return(pipeline)
# })
# 
# designMatrix <- do.call("rbind", designMatrix)
# designMatrix <- sapply(designMatrix, unlist)
# colnames(designMatrix)[[1]] <- "pipelines"
# colnames(designMatrix)[which(colnames(designMatrix) == "X")] <- "name"
# designMatrix <- as.data.frame(designMatrix)
# designMatrix <- merge(designMatrix, allMetrics, by=c("pipelines", "name"))
# 
# saveRDS(designMatrix, "/home/campbell/cfang/automl_scrna/data/uncorrected_designMatrix_scaled.RDS")

