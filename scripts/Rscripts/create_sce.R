suppressPackageStartupMessages({
  library(scater) # BioConductor
  library(SingleCellExperiment) # BioConductor
  library(DropletUtils) # BioConductor
  library(tidyverse) # CRAN
  library(here) # CRAN
  library(DT) # CRAN
  library(pheatmap) # CRAN
  library(here)
  library(Matrix)
  library(DESeq2)
  library(pipeComp)
})
knitr::opts_chunk$set(echo = TRUE,
                      cache = TRUE)
source(system.file("extdata", "scrna_alternatives.R", package="pipeComp"))
args = commandArgs(trailingOnly = TRUE)
cmmd <- paste("ls", args[1], sep=" ")


# reading in the data
file_address <- args[1]
list_exp <- paste("ls", file_address, sep = " ")
exp_files <- system(list_exp, intern = TRUE)
data_dir <- paste(file_address, exp_files[2], sep = "/")
data_col_dir <- paste(file_address, grep("col", exp_files, value=TRUE, fixed=TRUE), sep="/")
data_row_dir <- paste(file_address, grep("row", exp_files, value=TRUE, fixed=TRUE), sep="/")
cluster_data_dir <- paste(file_address, grep("clusters", exp_files, value=TRUE, fixed=TRUE), sep="/")
data_design <- paste(file_address, grep("experiment-design", exp_files, value=TRUE, fixed=TRUE), sep="/")
print(data_dir)
print(data_design)
print(data_col_dir)

data <- readMM(data_dir)
data_col <- read_tsv(data_col_dir, col_names = FALSE)
data_row <- read_tsv(data_row_dir, col_names = FALSE)
cluster_data <- read_tsv(cluster_data_dir, col_names = TRUE)
expdesign <- read_tsv(data_design)
expdesign = semi_join(expdesign[,1], data_col[,1], by = c("Assay" = "X1"))

# creating the SCE
data <- as.matrix(data)
data <- Matrix(data, sparse=TRUE)
sce <- SingleCellExperiment(list(counts = data), colData = as.vector(expdesign))
colnames(sce) <- as.vector(data_col$X1)
rownames(sce) <- as.vector(data_row$X1)

# formatting the clustering data
sel_k <- which(cluster_data$sel.K)
clusters <- cluster_data[sel_k, -c(1:2)]
clusters <- as.vector(as.matrix(clusters))
cell_names <- names(cluster_data)[-c(1,2)]
names(clusters) <- cell_names
sce$phenoid <- as.vector(clusters[colnames(sce)])

# computing QC metrics and metadata
sce <- addPerCellQC(sce)
sce <- addPerFeatureQC(sce)
sce <- add_meta(sce)

# renaming some colData to match what's used in pipeComp
pct_counts_in_top_50_features <- sce$pct_counts_top_50_features
sce$pct_counts_in_top_50_features <- pct_counts_in_top_50_features
sce$pct_counts_top_50_features <- NULL
sce$pct_counts_Mt <- sce$pct_Mt
sce$pct_Mt <- NULL
sce$pct_counts_in_top_20_features <- sce$percent.top_20
sce$percent.top_20 <- NULL

# saving the SCE as an RDS
file_address <- strsplit(file_address, "/")[[1]][[8]]
file_address <- paste(file_address, file_address, sep="/")
sce_filename <- paste(file_address, "-SCE.RDS", sep="")
sce_filename <- paste("/home/campbell/cfang/automl_scrna/data/experiments/", sce_filename, sep="")
print(sce_filename)
saveRDS(sce, sce_filename)




