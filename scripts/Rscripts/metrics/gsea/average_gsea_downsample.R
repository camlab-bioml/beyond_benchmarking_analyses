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
output <- args[[1]]
gsea <- read.csv(args[[2]])

numbers_only <- function(x) {
  print(sum(suppressWarnings(!is.na(as.numeric(as.character(x))))) == length(x))
}

print(gsea)

gsea$X1 <- NULL
gsea$X <- NULL

gsea$`X.1` <- NULL
gseatidy <- as_tibble(gsea)


pipelines <- gseatidy$pipelines
name <- gseatidy$name
gseatidy$pipelines <- NULL
gseatidy$name <- NULL

gseatidy$x <- NULL

print(gseatidy)
gsea_tidy <- gseatidy %>%
  mutate(across(starts_with("X", ignore.case=FALSE),as.numeric))%>%
  mutate(means=rowMeans(., na.rm=TRUE))
print(gsea_tidy)
  

gsea_tidy <- cbind(pipelines = pipelines, name=name, gsea_tidy)
print(gsea_tidy)

write.csv(gsea_tidy, output)
