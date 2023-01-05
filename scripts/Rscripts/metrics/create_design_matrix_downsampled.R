suppressPackageStartupMessages({
  library(tidyverse) # CRAN
  library(here) # CRAN
  library(DT) # CRAN
  library(here) # CRAN
  library(Matrix) # CRAN
  library(scran) # BioConductor
  library(purrr)
  library(dplyr)
})
removeDuplicateColumns <- function(x) {
  x <- x[!duplicated(as.list(x))]
  if ("pct_counts_in_top_20_features" %in% colnames(x)){
    print(colnames(x))
    x <- x %>%
      dplyr::rename(pct_counts_top_20_features = pct_counts_in_top_20_features)
  }
  return(x)
}

args <- commandArgs(trailingOnly = TRUE)
output <- args[[1]]
args <- args[-1]
allMetrics <- readRDS(args[[1]])
allMetricsDown <- readRDS(args[[2]])
expData <- as.data.frame(readRDS(args[[3]]))
downsampleExpData <- lapply(args[4:6], readRDS)

### Merge the coldata for each type of experiment
downsampleExpData <- lapply(downsampleExpData, removeDuplicateColumns)
print('downsample ok')
expData <- removeDuplicateColumns(expData)
print('expData ok ')

downsampleExpData <- do.call(rbind, downsampleExpData)
expData$downsampleType <- "None"

allExpData <- rbind(downsampleExpData, expData)
print(head(allExpData))

### Now merge the metrics for each type of experiment

allMetrics <- rbind(allMetrics, allMetricsDown)

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

allExpData <- as_tibble(allExpData)

designMatrix <- lapply(pipelines, function(pipeline){
  pipeline <- matrix(rep(design[which(design$pipelines == pipeline),], each = nrow(expData)), nrow=nrow(expData))
  pipeline <- cbind(as.data.frame(pipeline), allExpData)
  return(pipeline)
})

designMatrix <- do.call("rbind", designMatrix)
designMatrix <- sapply(designMatrix, unlist)
colnames(designMatrix)[[1]] <- "pipelines"
colnames(designMatrix)[which(colnames(designMatrix) == "X")] <- "name"
designMatrix <- as.data.frame(designMatrix)
designMatrix <- merge(designMatrix, allMetrics, by=c("pipelines", "name"))

print(designMatrix)
saveRDS(designMatrix, output)
