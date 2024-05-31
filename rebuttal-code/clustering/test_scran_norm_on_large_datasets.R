suppressPackageStartupMessages({
  library(SingleCellExperiment) # BioConductor
  library(DropletUtils) # BioConductor
  library(tidyverse) # CRAN
  library(here) # CRAN
  library(DT) # CRAN
  library(pheatmap) # CRAN
  library(here) # CRAN
  library(Matrix) # CRAN
  library(DESeq2) # BioConductor
  library(pipeComp) # BioConductor
  library(data.table)
  library(variancePartition) # BioConductor
  library(Cairo)
  library(Seurat) # BioConductor
  library(yaml) # CRAN
})

args = commandArgs(trailingOnly = TRUE)

source(system.file("extdata", "scrna_alternatives.R", package="pipeComp"))
pip_def <- scrna_pipeline(pipeClass = "seurat")

evalDummy <- function(...){
  return(1)
}
pip_def@evaluation$doublet <- evalDummy
pip_def@evaluation$filtering <- evalDummy
pip_def@evaluation$normalization <- evalDummy
pip_def@evaluation$selection <- evalDummy
pip_def@evaluation$dimreduction <- evalDummy
pip_def@evaluation$clustering <- evalDummy

pip_def@aggregation$doublet <- evalDummy
pip_def@aggregation$filtering <- evalDummy
pip_def@aggregation$normalization <- evalDummy
pip_def@aggregation$selection <- evalDummy
pip_def@aggregation$dimreduction <- evalDummy
pip_def@aggregation$clustering <- evalDummy
# Read in alternative parameters from yaml file
# #alternatives <- read_yaml(args[[1]])
# comb <- buildCombMatrix(alternatives)
#
# #print(comb)
# comb <- comb[ (comb$norm =="norm.scran") ,]
# print(comb)
alternatives <- list(
  doubletmethod=c("none"),
  filt=c("filt.lenient"),
  norm=c("norm.scran"),
  sel=c("sel.vst"),
  selnb=2000,
  dr=c("seurat.pca"),
  clustmethod=c("clust.seurat"),
  dims=c(10),
  resolution=c(0.01, 0.1)
)

print(names(alternatives))
#print(colnames(comb))
print(args)
# Take a file name from the command line and run pipecomp on it
filename <- strsplit(args[[2]], "/")[[1]][[4]]
filename <- strsplit(filename, "-SCE.RDS")[[1]][[1]]
print(filename)
filename <- paste(filename, filename, sep="/")
print(filename)

Sys.time()
Sys.Date()
res <- runPipeline(args[[2]], alternatives, pip_def,
                   output.prefix = paste0("/home/campbell/cfang/bb-rebuttal/results/pipecomp_outputs/test-scran",filename),
                   nthreads=1,
                   skipErrors=TRUE,
                   debug=TRUE,
                   saveEndResults=TRUE)

Sys.time()
Sys.Date()
