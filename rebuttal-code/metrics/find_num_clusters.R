suppressPackageStartupMessages({
  library(here) # CRAN
  library(DT) # CRAN
  library(here) # CRAN
  library(Matrix) # CRAN
  library(purrr)
  library(ggplot2)
})


exps = commandArgs(trailingOnly = TRUE)
nclusts <- matrix(,nrow=324, ncol=length(exps))
expnames <- c()
for (i in 1:length(exps)){
  exp.name <- strsplit(exps[[i]], split="/")[[1]][[3]]
  print(exp.name)
  expnames <- c(expnames, exp.name)
  clusts <- readRDS(exps[[i]])
  clusts <- lapply(clusts, unique)
  clusts <- lapply(clusts, length)
  pipelines <- names(clusts)
  names(clusts) <- NULL
  clusts <- unlist(clusts)
  print(length(clusts))
  nclusts[,i] <- t(clusts)
}
print(nclusts)
colnames(nclusts) <- expnames
rownames(nclusts) <- pipelines
print(nclusts)
write.csv(nclusts, file="/home/campbell/cfang/bb-rebuttal/results/pipecomp_outputs/num_clusters.csv")
