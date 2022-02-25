suppressPackageStartupMessages({
  library(scater) # BioConductor
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
alternatives <- read_yaml(args[[1]])

print(args)
# Take a file name from the command line and run pipecomp on it
filename <- strsplit(args[[2]], "/")[[1]][[9]]
filename <- strsplit(filename, "-SCE.RDS")[[1]][[1]]
print(filename)
filename <- paste(filename, filename, sep="/")
print(filename)
res <- runPipeline(args[[2]], alternatives, pip_def, output.prefix = paste0("/home/campbell/cfang/automl_scrna/results/pipecomp_outputs/",filename), nthreads=1, skipErrors=TRUE, debug=TRUE, saveEndResults=TRUE)
#saveRDS(res, file=paste0(filename, "-backup.RDS"))
