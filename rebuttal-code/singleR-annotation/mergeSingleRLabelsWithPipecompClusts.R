suppressPackageStartupMessages({
  #library(scater) # BioConductor
  library(SingleCellExperiment) # BioConductor
  #library(DropletUtils) # BioConductor
  library(tidyverse) # CRAN
  library(here) # CRAN
  library(DT) # CRAN
  #library(pheatmap) # CRAN
  library(here)
  library(Matrix)
})

args <- commandArgs(trailingOnly = TRUE)
singleRAnnot_path <- args[[1]]
clustersPaths <- args[-1]

labels <- readRDS(singleRAnnot_path)
cellNames <- rownames(labels)

labels <- labels[,c("labels")]
labels <- as.data.frame(labels)
labels$X <- cellNames

print(head(labels))

print(args[[1]])
expName <- strsplit(args[[1]], split="/")[[1]][[3]]
expName <- strsplit(expName, "-singleR-results.RDS")[[1]][[1]]
print(expName)

for (path in clustersPaths){
  clusters <- read.csv(path)

  # if (any(labels[,2] == "")){
  #   labels <- labels[-grep("^$", labels[,2]),]
  # }

  labels[,1] <- as.numeric(as.factor(labels[,1]))
  print(head(labels))
  print("clusters")
  print(head(clusters))

  clust_labels <- merge(clusters, labels, by="X")
  names(clust_labels)[3] <- "true_labels"
  print("clust_labels")
  print(head(clust_labels))

  run <- strsplit(path, split=".csv")[[1]][[1]]
  run <- strsplit(run, split="/")[[1]][[5]]
  print(run)

  clust_labels_name <- paste0(run, "-singleR-labelled.csv")
  clustDir <- paste0("/home/campbell/cfang/bb-rebuttal/results/pipecomp_outputs/singleR_labelled_clusters/", expName)
  clust_labels_name <- paste(clustDir, clust_labels_name, sep="/")
  print(clust_labels_name)

  write.csv(clust_labels, clust_labels_name)
}
