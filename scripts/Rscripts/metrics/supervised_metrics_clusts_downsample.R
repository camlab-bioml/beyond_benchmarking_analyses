suppressPackageStartupMessages({
  library(scater) # BioConductor
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
downsampleType <- args[[length(args)]]
print(downsampleType)
expDesignPath <- args[[1]]
clustersPaths <- args[-1]
clustersPaths <- clustersPaths[-length(clustersPaths)]

labels <- read.table(expDesignPath, sep="\t")
names <- labels[1,]
print(names)
cellTypes <- names[grepl("cell type", names)]
print("-----cell type------")
print(cellTypes)
authorLabels <- cellTypes[grepl("author", cellTypes)]
print(length(authorLabels))
if (length(authorLabels)!=1){
  print("empty authorLabels")
  authorLabels <- cellTypes[grepl("Factor", cellTypes)]
  authorLabels <- authorLabels[!grepl("Ontology Term", authorLabels)]
}
print("------authors------")
print(authorLabels)
print("final author label")
print(authorLabels[[1]])
authorLabels <- labels[,which(labels[1,]==authorLabels[[1]])]
labels <- cbind(labels[,1], authorLabels)
print(labels)
colnames(labels) <- c("X", "V2")

print(labels)
labels <- labels[-1,]

expName <- strsplit(expDesignPath, split="labels-")[[1]][[2]]
expName <- strsplit(expName, split=".tsv")[[1]][[1]]

for (path in clustersPaths){
  clusters <- read.csv(path)
  
  if (any(labels[,2] == "")){
    labels <- labels[-grep("^$", labels[,2]),]
  }
  
  labels[,2] <- as.numeric(as.factor(labels[,2]))
  #print(head(labels))
  #print("clusters")
  #print(head(clusters))
  
  clust_labels <- merge(clusters, labels, by="X")
  names(clust_labels)[3] <- "true_labels"
  #print("clust_labels")
  #print(head(clust_labels))
  
  run <- strsplit(path, split=".csv")[[1]][[1]]
  run <- strsplit(run, split="/")[[1]][[6]]
  print(run)
  
  clust_labels_name <- paste0(run, "-labelled.csv")
  clustDir <- paste("/home/campbell/cfang/automl_scrna/results/pipecomp_outputs/labelled_clusters", expName, downsampleType, sep="/")
  clust_labels_name <- paste(clustDir, clust_labels_name, sep="/")
  #print(clust_labels_name)
  
  if(nrow(clust_labels)==0){
    print(clusters[1,2])
    clust_labels <- as.data.frame(clust_labels)
    print(colnames(clust_labels))
    colnames(clust_labels)[[2]] <- clusters[1,2]
    #print(clust_labels)
  }
  print(colnames(clust_labels))
  write.csv(clust_labels, clust_labels_name)
}
