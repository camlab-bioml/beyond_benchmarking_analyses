suppressPackageStartupMessages({
  library(tidyverse) # CRAN
  library(here) # CRAN
  library(DT) # CRAN
  library(here) # CRAN
  library(Matrix) # CRAN
  library(plyr)
  library(dplyr)
})

args = commandArgs(trailingOnly = TRUE)

exp <- args[[1]]

cmmd <- paste("ls /home/campbell/cfang/automl_scrna/results/pipecomp_outputs/clusters/", exp, sep="")

clusts <- system(cmmd, intern=TRUE)
i <- 0
for (clust in clusts){
  i <- i+1
  clust_path <- paste("/home/campbell/cfang/automl_scrna/results/pipecomp_outputs/clusters/", exp, sep="")
  clust_path <- paste(clust_path, clust, sep="/")
  print(clust_path)
  clust <- read.csv(clust_path)
  pip <- colnames(clust)[[2]]
  print(pip)
  res <- strsplit(pip, "resolution.")[[1]][[2]]
  res <- strsplit(res, ".min")[[1]][[1]]
  print(res)
  if (res != "0.01"){
    write.csv(clust, sprintf("/home/campbell/cfang/automl_scrna/results/pipecomp_outputs/paperData/clusters/%s/%s-%s.csv", exp, exp, i))
  }
}

