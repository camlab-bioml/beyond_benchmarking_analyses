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
nclusts <- readRDS(args[[2]])
gsea <- read.csv(args[[3]])
gseaDownsample <- read.csv(args[[4]])

# combine the gsea results from both regular and downsampled experiments
#head(gsea)
#head(gseaDownsample)
gsea <- rbind(gsea, gseaDownsample)

# Flip sign of db
designMat$db <- designMat$db * -1

numbers_only <- function(x) {
  sum(suppressWarnings(!is.na(as.numeric(as.character(x))))) == length(x)
}

designMat <- designMat %>% mutate_if(numbers_only, as.numeric)

imputeMed <- function(x){
  x <- as.numeric(x)
  med <- median(x, na.rm=TRUE)
  x[is.na(as.numeric(x))] <- med
  return(as.numeric(x))
}

gsea$X <- NULL
gsea <- as_tibble(gsea)

print(any(is.na(gsea$means)))


gseaImp <- gsea %>%
  group_by(name) %>% 
  group_modify(~ as.data.frame(cbind(pipelines=.$pipelines, means=imputeMed(.$means))))%>%
  dplyr::rename(gsea = means)

expsWithNoClusters <- unique(gsea$name[which(is.na(gsea$means))])

gseaImp$pipelines <- gsub(pattern="\"", replacement="", gseaImp$pipelines)
gseaImp$pipelines <- gsub(pattern="\\;", replacement="\\.", gseaImp$pipelines)
gseaImp$pipelines <- gsub(pattern="\\=", replacement="\\.", gseaImp$pipelines)

gseaImp <- gseaImp %>%
  group_by(name) %>%
  group_modify(~ as.data.frame(cbind(pipelines=.$pipelines, scale(as.numeric(.$gsea)))))
colnames(gseaImp)[[3]] <- "gseaScaledImp"


designMat <- merge(designMat, gseaImp, by=c("pipelines", "name")) 

silScoresCorrected <- designMat %>%
  group_by(name) %>%
  group_modify(~ as.data.frame(cbind(silCor=resid(loess(as.numeric(.$sil) ~ as.numeric(.$nclusts), na.action=na.exclude), na.action=na.exclude), pipelines=.$pipelines, sil=.$sil)))

dbScoresCorrected <- designMat %>%
  group_by(name) %>%
  group_modify(~ as.data.frame(cbind(dbCor=resid(loess(as.numeric(.$db) ~ as.numeric(.$nclusts), na.action=na.exclude), na.action=na.exclude),pipelines=.$pipelines, db = .$db)))

chScoresCorrected <- designMat %>%
  group_by(name) %>%
  group_modify(~ as.data.frame(cbind(chCor=resid(loess(as.numeric(.$ch) ~ as.numeric(.$nclusts), na.action=na.exclude), na.action=na.exclude),pipelines=.$pipelines, ch=.$ch)))

silDesignMat <- silScoresCorrected %>%
  group_by(name) %>%
  group_modify(~ as.data.frame(cbind(pipelines=.$pipelines, sil = .$sil, silCor = .$silCor, silCorImp=imputeMed(.$silCor))))

dbDesignMat <- dbScoresCorrected %>% 
  group_by(name) %>%
  group_modify(~ as.data.frame(cbind(pipelines=.$pipelines, db = .$db, dbCor = .$dbCor, dbCorImp=imputeMed(.$dbCor))))

chDesignMat <- chScoresCorrected %>%
  group_by(name) %>%
  group_modify(~ as.data.frame(cbind(pipelines=.$pipelines, ch = .$ch, chCor = .$chCor, chCorImp=imputeMed(.$chCor))))


# Find the experiments that have NAs after correction and imputation
naexps <-unique(silDesignMat$name[which(is.na(silDesignMat$silCorImp))])
print(naexps)
print(length(naexps))

# correct these scores using y-mean(y) and then impute the median like above 
nascoresSil <- silDesignMat %>%
  filter(name %in% naexps)%>%
  select(name, ends_with("Imp"),  sil, pipelines)%>%
  group_by(name) %>%
  group_modify(~ as.data.frame(cbind(pipelines=.$pipelines, sil = .$sil, silCor = as.numeric(.$sil) - mean(as.numeric(.$sil), na.rm=TRUE))))%>%
  group_modify(~ as.data.frame(cbind(pipelines=.$pipelines, sil = .$sil, silCor = .$silCor, silCorImp=imputeMed(as.numeric(.$silCor)))))

nascoresDb <- dbDesignMat %>%
  filter(name %in% naexps)%>%
  select(name, ends_with("Imp"),  db, pipelines)%>%
  group_by(name) %>%
  group_modify(~ as.data.frame(cbind(pipelines=.$pipelines, db = .$db, dbCor = as.numeric(.$db) - mean(as.numeric(.$db), na.rm=TRUE))))%>%
  group_modify(~ as.data.frame(cbind(pipelines=.$pipelines, db = .$db, dbCor = .$dbCor, dbCorImp=imputeMed(as.numeric(.$dbCor)))))

nascoresCh <- chDesignMat%>%
  filter(name %in% naexps)%>%
  select(name, ends_with("Imp"),  ch, pipelines)%>%
  group_by(name) %>%
  group_modify(~ as.data.frame(cbind(pipelines=.$pipelines, ch = .$ch, chCor = as.numeric(.$ch) - mean(as.numeric(.$ch), na.rm=TRUE))))%>%
  group_modify(~ as.data.frame(cbind(pipelines=.$pipelines, ch = .$ch, chCor = .$chCor, chCorImp=imputeMed(as.numeric(.$chCor)))))


# rbind the datasets that have NAs with the original datasets and drop the 
# rows that have NAs in the imputed column because these are duplicates 
silScores <- rbind(silDesignMat, nascoresSil) %>%
  drop_na(silCorImp)%>%
  select(name, pipelines, silCorImp, silCor)

dbScores <- rbind(dbDesignMat, nascoresDb) %>%
  drop_na(dbCorImp)%>%
  select(name, pipelines, dbCorImp, dbCor)

chScores <- rbind(chDesignMat, nascoresCh) %>%
  drop_na(chCorImp)%>%
  select(name, pipelines, chCorImp, chCor)

# merge all corrected scores with the design matrix
designMat <- merge(designMat, silScores, by=c("name", "pipelines"))
designMat <- merge(designMat, dbScores, by=c("name", "pipelines"))
designMat <- merge(designMat, chScores, by=c("name", "pipelines"))

designMat <- designMat %>%
  filter(!(name %in% expsWithNoClusters))

saveRDS(designMat, "/home/campbell/cfang/automl_scrna/data/corrected_designMat_unscaled_downsampled.RDS")
