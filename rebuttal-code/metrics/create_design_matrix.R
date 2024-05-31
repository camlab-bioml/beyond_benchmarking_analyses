suppressPackageStartupMessages({
  library(tidyverse) # CRAN
  library(here) # CRAN
  library(DT) # CRAN
  library(here) # CRAN
  library(Matrix) # CRAN
  library(scran) # BioConductor
  library(purrr)
})
args <- commandArgs(trailingOnly = TRUE)
allMetrics <- readRDS(args[[1]])
expData <- readRDS(args[[2]])

print(head(allMetrics))

pipelines <- unique(allMetrics$pipelines)
filt <- lapply(pipelines, function(pipeline){
  stringency <- strsplit(pipeline, "filt.")[[1]][[3]]
  stringency <- strsplit(stringency, ".norm")[[1]][[1]]
  return(stringency)
})
filt <- unlist(filt)
dims  <- lapply(pipelines, function(pipeline){
  dim <- strsplit(pipeline, "dims.")[[1]][[2]]
  dim <- strsplit(dim, ".k")[[1]][[1]]
  return(dim)
})
dims <- unlist(dims)
norms  <- lapply(pipelines, function(pipeline){
  norm <- strsplit(pipeline, "norm.")[[1]][[3]]
  norm <- strsplit(norm, ".sel")[[1]][[1]]
  return(norm)
})
norms <- unlist(norms)
res <- lapply(pipelines, function(pipeline){
  clust.res <- strsplit(pipeline, "resolution.")[[1]][[2]]
  clust.res <- strsplit(clust.res, ".min")[[1]][[1]]
  return(clust.res)
})
res <- unlist(res)
design <- as.data.frame(cbind(pipelines=pipelines, filt=filt, dims=dims, norms=norms, res=res))

expData <- as_tibble(expData)
print(expData)

designMatrix <- lapply(pipelines, function(pipeline){
  pipeline <- matrix(rep(design[which(design$pipelines == pipeline),], each = nrow(expData)), nrow=nrow(expData))
  pipeline <- cbind(as.data.frame(pipeline), expData)
  return(pipeline)
})

designMatrix <- do.call("rbind", designMatrix)
designMatrix <- sapply(designMatrix, unlist)
colnames(designMatrix)[[1]] <- "pipelines"
colnames(designMatrix)[which(colnames(designMatrix) == "X")] <- "name"
designMatrix <- as.data.frame(designMatrix)
designMatrix <- merge(designMatrix, allMetrics, by=c("pipelines", "name"))

saveRDS(designMatrix, "/home/campbell/cfang/bb-rebuttal/data/large-samples_uncorrected_designMatrix_scaled.RDS")
