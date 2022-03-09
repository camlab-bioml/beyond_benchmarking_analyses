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
  library(matrixStats)
  library(stringr)
  library(stringi)
  library(ggforce)
  library(cowplot)
})

models_path <- "RFmodels.RDS"
designMatPath <- "uncorrected_designMatrix_scaled.RDS"
ariPath <- "ari_unscaled.csv"
numClustersPath <- "num_clusters.csv"

models <- readRDS(models_path)
designMat <- readRDS(designMatPath)
ari <- read.csv(ariPath)
nclusts <- read.csv(numClustersPath)

silRF <- models$silRF
dbRF <- models$dbRF
chRF <- models$chRF
gseaRF <- models$gseaRF

# Clean up designMat
# remove duplicate predictors
designMat <- designMat[,-which(colnames(designMat) %in% c("sum.1", "detected.1", "total.1", "pct_counts_in_top_20_feature", "percent.top_50", "total_counts","total_features"))]
designMat <- merge(designMat, avg_expr_scores, by="name")
designMat <- designMat %>%
  dplyr::rename(filt = V2.x)%>%
  dplyr::rename(dims = V3.x)%>%
  dplyr::rename(norm = V4.x)%>%
  dplyr::rename(res = V5.x)

#flip sign of db
designMat$db <- designMat$db * -1
designMat$pipelines <- gsub(pattern="\"", replacement="", designMat$pipelines)

colnames(nclusts) <- unlist(lapply(colnames(nclusts), gsub, pattern="*\\.", replacement="-"))
nclusts$X <- gsub(pattern="\"", replacement="", nclusts$X)
nclusts$X <- gsub(pattern="\\;", replacement="\\.", nclusts$X)
nclusts$X <- gsub(pattern="\\=", replacement="\\.", nclusts$X)
colnames(nclusts)[which(colnames(nclusts)=="X")] <- "pipelines"

# Merge nclusts with designMat
nclusts <- as_tibble(nclusts)
nclusts <- pivot_longer(nclusts, cols=which(grepl("^E-*", colnames(nclusts))))
designClusts <- merge(designMat, nclusts, by=c("pipelines", "name"))
colnames(designClusts)[which(colnames(designClusts)=="value")] <- "nclusts"
designMat <- designClusts

# Clean up ARI mat
ari$X <- NULL
colnames(ari) <- unlist(lapply(colnames(ari), gsub, pattern="*\\.supervised\\.csv", replacement=""))
colnames(ari) <- unlist(lapply(colnames(ari), gsub, pattern="*\\.", replacement="-"))

# Read in design mats for each metric
sildesignMatTest <- readRDS("silDesignTest.RDS")
sildesignMatTrain <- readRDS("silDesignTrain.RDS")
dbdesignMatTest <- readRDS("dbDesignTest.RDS")
dbdesignMatTrain <- readRDS("dbDesignTrain.RDS")
chdesignMatTest <- readRDS("chDesignTest.RDS")
chdesignMatTrain <- readRDS("chDesignTrain.RDS")
gseadesignMatTest <- readRDS("gseaDesignTest.RDS")
gseadesignMatTrain <- readRDS("gseaDesignTrain.RDS")

datasetR2s <- function(x){
  datasetR2s <- x %>%
    group_by(name) %>%
    group_modify(~ as.data.frame(R2(as.numeric(.$preds), as.numeric(.$Actual))))
  colnames(datasetR2s)[[2]] <- "R2"
  return(datasetR2s)
}

datasetCors <- function(x){
  datasetCors <- x %>%
    group_by(name) %>%
    group_modify(~ as.data.frame(cor(as.numeric(.$preds), as.numeric(.$Actual))))
  colnames(datasetCors)[[2]] <- "Cor"
  return(datasetCors)
}

createARIPlottingDfs <- function(preds, dataMat, ariMat){
  df <- cbind(dataMat, preds=preds)
  df <- df %>%
    merge(ariMat, by=c("pipelines", "name"))%>%
    dplyr::rename(Actual = ends_with("Imp")) %>%
    dplyr::rename(ARI = value)
  return(df)
}

createPlottingDfs <- function(preds, dataMat){
  df <- cbind(dataMat, preds=preds)
  df <- df %>%
    dplyr::rename(Actual = ends_with("Imp"))
  return(df)
}

createARIPlotMat <- function(df, ariMat, metric="", model=""){
  df <- df %>%
    merge(ariMat, by=c("pipelines", "name"))%>%
    dplyr::rename(ARI = value)
  ariCors <- df %>%
    group_by(name) %>%
    group_modify(~as.data.frame(cor(.$preds, .$ARI)))
  ariCors <- as.data.frame(ariCors)
  return(ariCors)
}

plotPredsARI <- function(df, ariMat, metric="", model=""){
  df <- df %>%
    merge(ariMat, by=c("pipelines", "name"))%>%
    dplyr::rename(ARI = value)
  print(df)
  ariCors <- df %>%
    group_by(name) %>%
    group_modify(~as.data.frame(cor(.$preds, .$ARI)))
  ariCors <- as.data.frame(ariCors)
  print(ariCors)
  ggplot(df, aes(x=ARI, y=preds))+
    geom_point(aes(colour=as.factor(res)))+
    facet_wrap(~name, scales="free")+
    xlab("Scaled ARI")+
    ylab(sprintf("Predicted %s", metric))+
    ggtitle(sprintf("Predicted %s from %s vs. ARI", metric, model))+
    geom_text(
      data    = ariCors,
      mapping = aes(x = Inf, y = Inf, label = paste("Cor:", round(V1, digits=3))),
      hjust   = 1.05,
      vjust   = 1.5
    )
}

plotPredActuals <- function(df, metric, model, R2labels,npage=1){
  for (i in 1:npage){
    p <- ggplot(data=df, aes(x=as.numeric(Actual), y=as.numeric(preds)))+
      geom_point(aes(colour=nclusts))+
      facet_wrap_paginate(~name, scales="free", page=i, nrow=4, ncol=4)+
      xlab("Actual")+
      ylab(sprintf("Predicted %s", metric))+
      ggtitle(sprintf("Predicted vs. Actual %s from %s", metric, model))+
      geom_text(
        data    = R2labels,
        mapping = aes(x = Inf, y = Inf, label = paste("Cor:", round(Cor, digits=3))),
        hjust   = 1.05,
        vjust   = 1.5,
        color = "red"
      )+ 
      theme(strip.background = element_rect(fill="white"),
            strip.text = element_text(face="bold"))
    print(p)
  }
  
}

silTestPlot <- createPlottingDfs(silRF$testPreds, sildesignMatTest)
dbTestPlot <- createPlottingDfs(dbRF$testPreds,dbdesignMatTest)
chTestPlot <- createPlottingDfs(chRF$testPreds, chdesignMatTest)
gseaTestPlot <- createPlottingDfs(gseaRF$testPreds, gseadesignMatTest) %>%
  dplyr::rename(Actual= gseaImpScale)

silTrainPlot <- createPlottingDfs(silRF$trainPreds, sildesignMatTrain)
dbTrainPlot <- createPlottingDfs(dbRF$trainPreds, dbdesignMatTrain)
chTrainPlot <- createPlottingDfs(chRF$trainPreds, chdesignMatTrain)
gseaTrainPlot <- createPlottingDfs(gseaRF$trainPreds, gseadesignMatTrain) %>%
  dplyr::rename(Actual= gseaImpScale)

silTestDatasetCors <- datasetCors(silTestPlot)
dbTestDatasetCors <- datasetCors(dbTestPlot)
chTestDatasetCors <- datasetCors(chTestPlot)
gseaTestDatasetCors <- datasetCors(gseaTestPlot)

silTrainDatasetCors <- datasetCors(silTrainPlot)
dbTrainDatasetCors <- datasetCors(dbTrainPlot)
chTrainDatasetCors <- datasetCors(chTrainPlot)
gseaTrainDatasetCors <- datasetCors(gseaTrainPlot)


makeColdataDF <- function(R2df){
  dMat <- designMat %>% 
    as_tibble() %>% 
    mutate_if(numbers_only, as.numeric) %>%
    mutate_if(is.numeric, scale) %>%
    dplyr::slice(which(.$name %in% R2df$name)) %>%
    select(c("name", "detected", "total", "ncells", "ngenes", contains("_"))) %>%
    distinct() %>%
    pivot_longer(!name, names_to="dataType") %>%
    print()%>%
    merge(R2df, by="name") #%>%
  #mutate(PredPower = if_else(Cor >= 0.3, "R2 >= 0.1", "R2 < 0.1"))
  return(dMat)
}

## Predicted vs Actual plots ###
# pdf("tunedRFsPredsActualsCorsMar2.pdf")
# plotPredActuals(silTestPlot, metric="Sil (Test set)", model="RF (mtry=35, ntree=150)", silTestDatasetCors, npage=2)
# plotPredActuals(silTrainPlot, metric="Sil (Train set)", model="RF (mtry=35, ntree=150)", silTrainDatasetCors, npage=4)
# plotPredActuals(dbTestPlot, metric="DB (Test set)", model="RF (mtry=33, ntree=150)", dbTestDatasetCors, npage=2)
# plotPredActuals(dbTrainPlot, metric="DB (Train set)", model="RF (mtry=33, ntree=150)", dbTrainDatasetCors, npage=4)
# plotPredActuals(chTestPlot, metric="CH (Test set)", model="RF (mtry=32, ntree=150)", chTestDatasetCors,npage=2)
# plotPredActuals(chTrainPlot, metric="CH (Train set)", model="RF (mtry=32, ntree=150)", chTrainDatasetCors,npage=4)
# plotPredActuals(gseaTestPlot, metric="GSEA (Test set)", model="RF (mtry=32, ntree=150)", gseaTestDatasetCors,npage=2)
# plotPredActuals(gseaTrainPlot, metric="GSEA (Train set)", model="RF (mtry=32, ntree=150)", gseaTrainDatasetCors,npage=4)
# dev.off()


### Predicted vs ARI plots ###
# pdf("tunedRFsPredsARIMar8.pdf")
# plotPredsARI(silTestPlot, ari_tbl, "Sil", "RF (mtry=35, ntree=150)")
# plotPredsARI(dbTestPlot, ari_tbl, "DB", "RF (mtry=33, ntree=150)")
# plotPredsARI(chTestPlot, ari_tbl, "CH", "RF (mtry=32, ntree=150)")
# plotPredsARI(gseaTestPlot, ari_tbl, "GSEA", "RF (mtry=32, ntree=150)")
# dev.off()
### ColData ARI Correlation heatmap ###
ari_tbl <- ari %>%
  as_tibble() %>%
  mutate_if(is.numeric, scale) %>%
  pivot_longer(starts_with("E"))
colnames(ari_tbl)[[3]] <- "value"

silARICor <- createARIPlotMat(silTestPlot, ari_tbl)
dbARICor <- createARIPlotMat(dbTestPlot, ari_tbl)
chARICor <- createARIPlotMat(chTestPlot, ari_tbl)
gseaARICor <- createARIPlotMat(gseaTestPlot, ari_tbl)

allARIdf <- dplyr::bind_rows(list(sil=silARICor, db=dbARICor, ch=chARICor, gsea=gseaARICor), .id = 'metric')

silARI <- makeColdataDF(silARICor)
dbARI <- makeColdataDF(dbARICor)
chARI <- makeColdataDF(chARICor)
gseaARI <- makeColdataDF(gseaARICor)

silColdataARI <- silARI %>%
  group_by(dataType) %>%
  group_modify(~ as.data.frame(cor(as.numeric(.$value), as.numeric(.$V1))))
colnames(silColdataARI)[[2]] <- "silCor"

dbColdataARI <- dbARI %>%
  group_by(dataType) %>%
  group_modify(~ as.data.frame(cor(as.numeric(.$value), as.numeric(.$V1))))
colnames(dbColdataARI)[[2]] <- "dbCor"

chColdataARI <- chARI %>%
  group_by(dataType) %>%
  group_modify(~ as.data.frame(cor(as.numeric(.$value), as.numeric(.$V1))))
colnames(chColdataARI)[[2]] <- "chCor"

gseaColdataARI <- gseaARI %>%
  group_by(dataType) %>%
  group_modify(~ as.data.frame(cor(as.numeric(.$value), as.numeric(.$V1))))
colnames(gseaColdataARI)[[2]] <- "gseaCor"



corHmARI <- merge(silColdataARI, dbColdataARI, by="dataType")
corHmARI <- merge(corHmARI, chColdataARI, by="dataType")
corHmARI<- merge(corHmARI, gseaColdataARI, by="dataType")
rownames(corHmARI) <- corHmARI$dataType
corHmARI$dataType <- NULL
col_ha <- HeatmapAnnotation(avgCor = anno_barplot(rowMeans(corHmARI)))
ComplexHeatmap::Heatmap(t(as.matrix(corHmARI)), top_annotation = col_ha,column_title="Correlation between colData and dataset ARI correlations")


### ColData ARI Cor boxplot ###
wilcoxLabelsARI <- c(wilcox.test(chARICor$V1)$p.value,
                     wilcox.test(dbARICor$V1)$p.value,
                     wilcox.test(gseaARICor$V1)$p.value,
                     wilcox.test(silARICor$V1)$p.value )
wilcoxLabelsARI <- as.data.frame(wilcoxLabels)
wilcoxLabelsARI <- cbind(wilcoxLabelsARI, metric=c("ch","db","gsea","sil"))
wilcoxLabelsARI[,2] <-NULL


ggplot(allARIdf, aes(x=metric, y=V1))+
  geom_boxplot()+
  xlab("Metirc")+
  ylab("Correlation with ARI")+
  ggtitle("Correlation of RF Predictions with ARI")+
  geom_text(
    data    = wilcoxLabelsARI,
    mapping = aes(x = metric, y = Inf, label=round(wilcoxLabels,5)),
    hjust   = 1.05,
    vjust   = 1.5,
    color = "red"
  )+
  cowplot::theme_cowplot()

### ColData Correlation Heatmap ###
silCol <- makeColdataDF(silTestDatasetCors)
dbCol <- makeColdataDF(dbTestDatasetCors)
chCol <- makeColdataDF(chTestDatasetCors)
gseaCol <- makeColdataDF(gseaTestDatasetCors)

alldf <- dplyr::bind_rows(list(sil=silCol, db=dbCol, ch=chCol, gsea=gseaCol), .id = 'metric')
alldf <- alldf #%>%
  #mutate(xLabel = paste(metric, PredPower, sep=" "))

silColdataCor <- silCol %>%
  group_by(dataType) %>%
  group_modify(~ as.data.frame(cor(as.numeric(.$value), as.numeric(.$Cor))))
colnames(silColdataCor)[[2]] <- "silCor"

dbColdataCor <- dbCol %>%
  group_by(dataType) %>%
  group_modify(~ as.data.frame(cor(as.numeric(.$value), as.numeric(.$Cor))))
colnames(dbColdataCor)[[2]] <- "dbCor"

chColdataCor <- chCol %>%
  group_by(dataType) %>%
  group_modify(~ as.data.frame(cor(as.numeric(.$value), as.numeric(.$Cor))))
colnames(chColdataCor)[[2]] <- "chCor"

gseaColdataCor <- gseaCol %>%
  group_by(dataType) %>%
  group_modify(~ as.data.frame(cor(as.numeric(.$value), as.numeric(.$Cor))))
colnames(gseaColdataCor)[[2]] <- "gseaCor"

corHmDf <- merge(silColdataCor, dbColdataCor, by="dataType")
corHmDf <- merge(corHmDf, chColdataCor, by="dataType")
corHmDf <- merge(corHmDf, gseaColdataCor, by="dataType")
rownames(corHmDf) <- corHmDf$dataType
corHmDf$dataType <- NULL
col_ha <- HeatmapAnnotation(avgCor = anno_barplot(rowMeans(corHmDf)))
ComplexHeatmap::Heatmap(t(as.matrix(corHmDf)), top_annotation = col_ha,column_title="Correlation between colData and dataset correlations")

### ColData Correlation Boxplot ###
corbpdf <- merge(silTestDatasetCors, dbTestDatasetCors, by="name", suffix=c("Sil","Db"))
chgsea <- merge(chTestDatasetCors, gseaTestDatasetCors, by="name", suffix=c("Ch", "Gsea"))
corbpdf <- merge(corbpdf, chgsea, by="name")

corbpdf <- corbpdf %>%
  as_tibble(corbpdf) %>%
  pivot_longer(cols=starts_with("Cor"), names_to="metric")
  
wilcoxLabels <- c(wilcox.test(chTestDatasetCors$Cor)$p.value,
                    wilcox.test(dbTestDatasetCors$Cor)$p.value,
                    wilcox.test(gseaTestDatasetCors$Cor)$p.value,
                    wilcox.test(silTestDatasetCors$Cor)$p.value )
wilcoxLabels <- as.data.frame(wilcoxLabels)
wilcoxLabels <- cbind(wilcoxLabels, metric=c("Ch","Db","Gsea","Sil"))
corbpdf$metric <- sub("^Cor", "", corbpdf$metric)
ggplot(corbpdf, aes(x=metric, y=value))+
  geom_boxplot()+
  xlab("Metric")+
  ylab("Correlation between predictions and labels")+
  ggtitle("Boxplot of correlations from RF on test datasets")+
  geom_text(
    data    = wilcoxLabels,
    mapping = aes(x = metric, y = Inf, label=round(wilcoxLabels,5)),
    hjust   = 1.05,
    vjust   = 1.5,
    color = "red"
  )+
  cowplot::theme_cowplot()








# pdf("RFPredictiveColdataBP.pdf")
# for (i in 1:4){
#   p <- ggplot(silCol, aes(x=PredPower, y=value))+
#     geom_boxplot()+
#     facet_wrap(~dataType, scale = "free")+
#     facet_wrap_paginate(~dataType, scale = "free",ncol=2, nrow=3, page=i)+
#     ggtitle("Coldata for Silhouette RF (mtry=15, ntree=150")+
#     xlab("R2 value")
#   print(p)
# }
# for (i in 1:4){
#   p <- ggplot(dbCol, aes(x=PredPower, y=value))+
#     geom_boxplot()+
#     facet_wrap(~dataType, scale = "free")+
#     facet_wrap_paginate(~dataType, scale = "free",ncol=2, nrow=3, page=i)+
#     ggtitle("Coldata for DB RF (mtry=15, ntree=150")+
#     xlab("R2 value")
#   print(p)
# }
# for (i in 1:4){
#   p <- ggplot(chCol, aes(x=PredPower, y=value))+
#     geom_boxplot()+
#     facet_wrap(~dataType, scale = "free")+
#     facet_wrap_paginate(~dataType, scale = "free",ncol=2, nrow=3, page=i)+
#     ggtitle("Coldata for CH RF (mtry=15, ntree=150")+
#     xlab("R2 value")
#   print(p)
# }
# dev.off()
# # silCol <- makeColdataDF(silTestWell)
# # dbCol <- makeColdataDF(dbTestWell)
# # chCol <- makeColdataDF(chTestWell)
# 
# designtbl <- designMat %>% 
#   as_tibble() %>% 
#   select(c("name", "detected", "total", "ncells", "ngenes", contains("_"))) %>%
#   distinct() %>%
#   pivot_longer(!name, names_to="dataType")

# pdf("RFPredictiveColdata.pdf")
# for (i in 1:4){
#   p<- ggplot(designtbl, aes(x=name, y=value))+
#     geom_point()+
#     facet_wrap_paginate(~dataType, scale = "free",ncol=2, nrow=3, page=i)+
#     geom_point(data=silCol, aes(x=name, y=value), color="red")+
#     #geom_text(merge(silTestWell, silCol, by="name"), mapping=aes(x = name, y= value, label=round(R2,digits=2), colour="R2 >= 0.1"), color="blue")+
#     theme(axis.text.x = element_text(angle = 45, size=2))+
#     ggtitle("Coldata for Silhouette RF (mtry=15, ntree=150)")
#   print(p)
# }
# for (i in 1:4){
#   p <- ggplot(designtbl, aes(x=name, y=value))+
#     geom_point()+
#     facet_wrap_paginate(~dataType, scale = "free",ncol=2, nrow=3, page=i)+
#     geom_point(data=dbCol, aes(x=name, y=value), color="red")+
#     #geom_text(merge(dbTestWell, dbCol, by="name"), mapping=aes(x = name, y= value, label=round(R2,digits=2), colour="R2 >= 0.1"), color = "blue")+
#     theme(axis.text.x = element_text(angle = 45, size=2))+
#     ggtitle("Coldata for DB RF (mtry=15, ntree=150)")
#   print(p)
# }
# 
# for (i in 1:4){
#   p <- ggplot(designtbl, aes(x=name, y=value))+
#     geom_point()+
#     facet_wrap_paginate(~dataType, scale = "free",ncol=2, nrow=3, page=i)+
#     geom_point(data=chCol, aes(x=name, y=value), color="red")+
#     #geom_text(merge(chTestWell, chCol, by="name"), mapping=aes(x = name, y= value, label=round(R2,digits=2), colour="R2 >= 0.1"), color = "blue")+
#     theme(axis.text.x = element_text(angle = 45, size=2))+
#     ggtitle("Coldata for CH RF (mtry=15, ntree=150)")
#   print(p)
# }
# dev.off()


# sil <- readRDS("RFTunedSil.RDS")
# db <- readRDS("RFTunedDB.RDS")
# ch <- readRDS("RFTunedCH.RDS")
# 
# pdf("RFtuningRMSE.pdf")
# plot(sil, main="Silhouette")
# plot(db, main="DB")
# plot(ch, main="CH")
# dev.off()