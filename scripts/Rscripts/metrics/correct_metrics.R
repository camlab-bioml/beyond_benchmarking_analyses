suppressPackageStartupMessages({
  library(tidyverse) # CRAN
  library(here) # CRAN
  library(DT) # CRAN
  library(here) # CRAN
  library(Matrix) # CRAN
  library(purrr)
  library(ggplot2)
  library(Metrics)
  library(gridExtra)
  library(caret)
  library(matrixStats)
  library(randomForest)
  library(pcaMethods)
  library(stringr)
  library(stringi)
})
args <- commandArgs(trailingOnly=TRUE)

designMat <- readRDS(args[[1]])
ari <- readRDS(args[[2]])
nclusts <- readRDS(args[[3]])
gsea <- read.csv(args[[4]])

# Flip sign of db
designMat$db <- designMat$db * -1

numbers_only <- function(x) {
  sum(suppressWarnings(!is.na(as.numeric(as.character(x))))) == length(x)
}

print(head(designMat))
designMat <- designMat %>% mutate_if(numbers_only, as.numeric)

print(head(designMat))
imputeMed <- function(x){
  x <- as.numeric(x)
  med <- median(x, na.rm=TRUE)
  x[is.na(as.numeric(x))] <- med
  return(as.numeric(x))
}

gsea$X <- NULL
gsea <- as_tibble(gsea)

gseaImp <- gsea %>%
  group_by(name) %>% 
  group_modify(~ as.data.frame(cbind(pipelines=.$pipelines, means=imputeMed(.$means))))%>%
  dplyr::rename(gsea = means)

gseaImp$pipelines <- gsub(pattern="\"", replacement="", gseaImp$pipelines)
gseaImp$pipelines <- gsub(pattern="\\;", replacement="\\.", gseaImp$pipelines)
gseaImp$pipelines <- gsub(pattern="\\=", replacement="\\.", gseaImp$pipelines)

gseaImp <- gseaImp %>%
  group_by(name) %>%
  group_modify(~ as.data.frame(cbind(pipelines=.$pipelines, scale(as.numeric(.$gsea)))))
colnames(gseaImp)[[3]] <- "gseaScaledImp"

print(head(gseaImp))
print(head(designMat))

designMat <- merge(designMat, gseaImp, by=c("pipelines", "name")) 

print(head(designMat))
silScoresCorrected <- designMat %>%
  group_by(name) %>%
  group_modify(~ as.data.frame(cbind(silCor=resid(loess(as.numeric(.$sil) ~ as.numeric(.$nclusts), na.action=na.exclude), na.action=na.exclude), pipelines=.$pipelines)))

dbScoresCorrected <- designMat %>%
  group_by(name) %>%
  group_modify(~ as.data.frame(cbind(dbCor=resid(loess(as.numeric(.$db) ~ as.numeric(.$nclusts), na.action=na.exclude), na.action=na.exclude),pipelines=.$pipelines)))

chScoresCorrected <- designMat %>%
  group_by(name) %>%
  group_modify(~ as.data.frame(cbind(chCor=resid(loess(as.numeric(.$ch) ~ as.numeric(.$nclusts), na.action=na.exclude), na.action=na.exclude),pipelines=.$pipelines)))

designMat <- merge(designMat, silScoresCorrected, by=c("name","pipelines"))
designMat <- merge(designMat, dbScoresCorrected, by=c("name","pipelines"))
designMat <- merge(designMat, chScoresCorrected, by=c("name","pipelines"))

silDesignMat <- designMat %>%
  group_by(name) %>%
  group_modify(~ as.data.frame(cbind(pipelines=.$pipelines, silCorImp=imputeMed(.$silCor))))

dbDesignMat <- designMat %>% 
  group_by(name) %>%
  group_modify(~ as.data.frame(cbind(pipelines=.$pipelines, dbCorImp=imputeMed(.$dbCor))))

chDesignMat <- designMat %>%
  group_by(name) %>%
  group_modify(~ as.data.frame(cbind(pipelines=.$pipelines, chCorImp=imputeMed(.$chCor))))

designMat <- merge(designMat, silDesignMat, by=c("name", "pipelines"))
designMat <- merge(designMat, dbDesignMat, by=c("name", "pipelines"))
designMat <- merge(designMat, chDesignMat, by=c("name", "pipelines"))
print(head(designMat))

saveRDS(designMat, "/home/campbell/cfang/automl_scrna/data/corrected_designMat_unscaled.RDS")
