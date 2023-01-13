suppressPackageStartupMessages({
  library(tidyverse) # CRAN
  library(here) # CRAN
  library(DT) # CRAN
  library(here) # CRAN
  library(Matrix) # CRAN
  library(purrr)
  library(stringr)
  library(stringi)
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)
gsea <- read.csv(args[[1]])

numbers_only <- function(x) {
  print(sum(suppressWarnings(!is.na(as.numeric(as.character(x))))) == length(x))
}


#gsea <- read.csv("E-HCAD-9-gsea.csv")
gsea$X1 <- NULL
gsea$X <- NULL
gseatidy <- as_tibble(gsea)
  
pipelines <- gseatidy$pipelines
name <- gseatidy$name
gseatidy$pipelines <- NULL
gseatidy$name <- NULL

gsea_tidy <- gseatidy %>%
  mutate(across(starts_with("X", ignore.case=FALSE),as.numeric))%>%
  mutate(means=rowMeans(., na.rm=TRUE))

gsea_tidy <- cbind(pipelines = pipelines, name=name, gsea_tidy)
exp.name <- strsplit(args[[1]], split="/")[[1]][[9]]
exp.name <- strsplit(exp.name, split="-gsea.csv")[[1]][[1]]
file.name <- paste0("/home/campbell/cfang/automl_scrna/results/gsea_results/gsea_matrices/", exp.name)
file.name <- paste0(file.name, "-gsea-averaged.csv")
print(file.name)
write.csv(gsea_tidy, file.name)
  