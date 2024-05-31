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


gsea$X1 <- NULL
gsea$X <- NULL
gseatidy <- as_tibble(gsea)

pipelines <- gseatidy$pipelines
name <- gseatidy$name
gseatidy$pipelines <- NULL
gseatidy$name <- NULL

print(head(gseatidy))
gseatidy$x <- NULL
gsea_tidy <- gseatidy %>%
  mutate(across(starts_with("X", ignore.case=FALSE),as.numeric))%>%
  mutate(means=rowMeans(., na.rm=TRUE))

gsea_tidy <- cbind(pipelines = pipelines, name=name, gsea_tidy)
exp.name <- strsplit(args[[1]], split="/")[[1]][[4]]
exp.name <- strsplit(exp.name, split="-gsea.csv")[[1]][[1]]
file.name <- paste0("/home/campbell/cfang/bb-rebuttal/results/gsea_results/gsea_matrices/", exp.name)
file.name <- paste0(file.name, "-gsea-averaged.csv")
print(file.name)
write.csv(gsea_tidy, file.name)
