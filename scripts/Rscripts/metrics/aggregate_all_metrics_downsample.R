suppressPackageStartupMessages({
  library(tidyverse) # CRAN
  library(here) # CRAN
  library(DT) # CRAN
  library(here) # CRAN
  library(Matrix) # CRAN
  library(scran) # BioConductor
  library(purrr)
})
expsMatrices <- commandArgs(trailingOnly = TRUE)
sil <- c()
db <- c()
ch <- c()

standardizePipeNames <- function(x){
  x <- gsub("\\+", ".", x)
  x <- gsub("\\=", ".", x)
  return(x)
}

for (i in 1:length(expsMatrices)){
  metrics <- read.csv(expsMatrices[[i]])
  colnames(metrics)[[1]] <- "pipelines"
  expName <- strsplit(expsMatrices[[i]], "-metrics.csv")[[1]][[1]]
  expName <- strsplit(expName, "metrics_matrices/")[[1]][[2]]
  #print(expName)
  #print(dim(metrics))
  names(metrics)[[1]] <- "index"
  #print(head(metrics))
  #print(dim(metrics))
  expNames <- rep(expName, nrow(metrics))
  if(any(lapply(metrics$pipelines, nchar) <= 20)){
    print(expName)
    print(which(lapply(metrics$pipelines, nchar) <= 20))
    print(metrics$pipelines[which(lapply(metrics$pipelines, nchar) <= 20)])
  }
  metrics$pipelines <- unlist(lapply(metrics$pipelines, standardizePipeNames))
  
  sil[[i]] <- data.frame("pipelines"=metrics$pipelines, "sil"=metrics$sil, "exp"=expNames)
  db[[i]] <- data.frame("pipelines"=metrics$pipelines, "db"=metrics$db, "exp"=expNames)
  ch[[i]] <- data.frame("pipelines"=metrics$pipelines, "ch"=metrics$ch, "exp"=expNames)
}
# do.call(rbind(sil)), then make it a tibble and pivot
silMat <- as_tibble(do.call(rbind, sil)) #%>%
  #mutate(row=row_number())

#print(silMat$pipelines)

dbMat <- as_tibble(do.call(rbind, db))#%>%
  #mutate(row=row_number())

chMat <- as_tibble(do.call(rbind, ch))#%>%
  #mutate(row=row_number())

silMat <- pivot_wider(silMat, names_from = "exp", values_from="sil")
dbMat <- pivot_wider(dbMat, names_from = "exp", values_from="db")
chMat <- pivot_wider(chMat, names_from = "exp", values_from="ch")

silName <- "/home/campbell/cfang/automl_scrna/results/metrics/sil_unscaled_downsample.csv"
dbName <- "/home/campbell/cfang/automl_scrna/results/metrics/db_unscaled_downsample.csv"
chName <- "/home/campbell/cfang/automl_scrna/results/metrics/ch_unscaled_downsample.csv"

print(head(silMat))
print(head(dbMat))
print(head(chMat))

write.csv(silMat, silName)
write.csv(dbMat, dbName)
write.csv(chMat, chName)

sils <- rep("sil", nrow(silMat))
dbs <- rep("db", nrow(dbMat))
chs <- rep("ch", nrow(chMat))

silMat <- cbind(silMat, score=sils)
dbMat <- cbind(dbMat, score=dbs)
chMat <- cbind(chMat, score=chs)

# Remove 0.01 resolution pipelines
silMat <- silMat[which(!str_detect(silMat$pipelines, "resolution.0.01")),]
dbMat <- dbMat[which(!str_detect(dbMat$pipelines, "resolution.0.01")),]
chMat <- chMat[which(!str_detect(chMat$pipelines, "resolution.0.01")),]

#Scale each dataset
silMat[,sapply(silMat, is.numeric)] <- scale(silMat[,sapply(silMat, is.numeric)])
dbMat[,sapply(dbMat, is.numeric)] <- scale(dbMat[,sapply(dbMat, is.numeric)])
chMat[,sapply(chMat, is.numeric)] <- scale(chMat[,sapply(chMat, is.numeric)])

# Now pivot metrics to 288*86 x 3
allMetrics <- as_tibble(rbind(silMat, dbMat, chMat))
allMetrics <- pivot_longer(allMetrics, cols=starts_with("E-"))
allMetrics <- pivot_wider(allMetrics, names_from="score")

print(head(allMetrics))
saveRDS(allMetrics, "/home/campbell/cfang/automl_scrna/results/metrics/all_metrics_downsample.RDS")


# [1] "E-CURD-10-SCE-50"
# [1] "E-ENAD-17-SCE-50"
# [1] "E-ENAD-17-SCE-75"
# [1] "E-HCAD-25-SCE-downsampleCounts-50"
# [1] "E-HCAD-35-SCE-downsampleCounts-50"
# [1] "E-MTAB-4850-SCE-50"
# [1] "E-MTAB-4850-SCE-75"
# [1] "E-MTAB-6379-SCE-50"
# [1] "E-MTAB-7303-SCE-50"