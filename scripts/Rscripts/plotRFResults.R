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
  library(ComplexHeatmap)
})
numbers_only <- function(x) {
  sum(suppressWarnings(!is.na(as.numeric(as.character(x))))) == length(x)
}


models <- readRDS(args[[1]])
designMat <- readRDS(args[[2]])
ari <- readRDS(args[[3]])
RFType <- args[[4]]

nclusts <- designMat %>%
  select(pipelines, name, nclusts)


silRF <- models$sil
dbRF <- models$db
chRF <- models$ch
gseaRF <- models$gsea


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
  df <- cbind(dataMat, preds=preds) %>%
    merge(nclusts, by=c("pipelines", "name"))
  return(df)
}

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

silTestPlot <- createPlottingDfs(silRF$testPreds, silRF$testData)
dbTestPlot <- createPlottingDfs(dbRF$testPreds, dbRF$testData)
chTestPlot <- createPlottingDfs(chRF$testPreds, chRF$testData)
gseaTestPlot <- createPlottingDfs(gseaRF$testPreds, gseaRF$testData)

silTrainPlot <- createPlottingDfs(silRF$trainPreds, silRF$trainData)
dbTrainPlot <- createPlottingDfs(dbRF$trainPreds, dbRF$trainData)
chTrainPlot <- createPlottingDfs(chRF$trainPreds, chRF$trainData)
gseaTrainPlot <- createPlottingDfs(gseaRF$trainPreds, gseaRF$trainData)

silTestDatasetCors <- datasetCors(silTestPlot)
dbTestDatasetCors <- datasetCors(dbTestPlot)
chTestDatasetCors <- datasetCors(chTestPlot)
gseaTestDatasetCors <- datasetCors(gseaTestPlot)

silTrainDatasetCors <- datasetCors(silTrainPlot)
dbTrainDatasetCors <- datasetCors(dbTrainPlot)
chTrainDatasetCors <- datasetCors(chTrainPlot)
gseaTrainDatasetCors <- datasetCors(gseaTrainPlot)

testCors <- dplyr::bind_rows(list(sil=silTestDatasetCors, db=dbTestDatasetCors, ch=chTestDatasetCors, gsea=gseaTestDatasetCors), .id="metric")%>%
  mutate(metric = toupper(metric))

filename <- sprintf("/home/campbell/cfang/scrna_automl/results/figures/RFPlots/RFPlots_{{RFType}}Params.RDS")
pdf(filename)
ggplot(testCors, aes(x=tidytext::reorder_within(name, -Cor, metric), y=Cor))+
  geom_bar(stat="identity")+
  ggtitle(sprintf("Correlation between RF Predictions and Labels on Test"))+
  xlab("Dataset")+
  ylab("Correlation")+
  cowplot::theme_cowplot()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=6))+
  tidytext::scale_x_reordered()+
  facet_wrap(~ metric, scale="free")+
  theme(strip.background = element_rect(fill="white"),
        strip.text = element_text(face="bold"))
#dev.off()

silParams = sprintf("RF (mtry=%s, ntree=%s)", silRF$mtry, silRF$ntree)
dbParams = sprintf("RF (mtry=%s, ntree=%s)", dbRF$mtry, dbRF$ntree)
chParams = sprintf("RF (mtry=%s, ntree=%s)", chRF$mtry, chRF$ntree)
gseaParams = sprintf("RF (mtry=%s, ntree=%s)", gseaRF$mtry, gseaRF$ntree)

# Predicted vs Actual plots ###
plotPredActuals(silTestPlot, metric="Sil (Test set)", model=silParams, silTestDatasetCors, npage=2)
plotPredActuals(silTrainPlot, metric="Sil (Train set)", model=silParams, silTrainDatasetCors, npage=4)
plotPredActuals(dbTestPlot, metric="DB (Test set)", model=dbParams, dbTestDatasetCors, npage=2)
plotPredActuals(dbTrainPlot, metric="DB (Train set)", model=dbParams, dbTrainDatasetCors, npage=4)
plotPredActuals(chTestPlot, metric="CH (Test set)", model=chParams, chTestDatasetCors,npage=2)
plotPredActuals(chTrainPlot, metric="CH (Train set)", model=chParams, chTrainDatasetCors,npage=4)
plotPredActuals(gseaTestPlot, metric="GSEA (Test set)", model=gseaParams, gseaTestDatasetCors,npage=2)
plotPredActuals(gseaTrainPlot, metric="GSEA (Train set)", model=gseaParams, gseaTrainDatasetCors,npage=4)

## Predicted vs ARI plots ###
plotPredsARI(silTestPlot, ari_tbl, "Sil", silParams)
plotPredsARI(dbTestPlot, ari_tbl, "DB", dbParams)
plotPredsARI(chTestPlot, ari_tbl, "CH", chParams)
plotPredsARI(gseaTestPlot, ari_tbl, "GSEA", gseaParams)

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

ComplexHeatmap::Heatmap(t(as.matrix(corHmARI)), top_annotation = col_ha,column_title="Correlation between colData and dataset ARI correlations", cluster_rows = FALSE)


### ColData ARI Cor boxplot ###
wilcoxLabelsARI <- c(wilcox.test(chARICor$V1)$p.value,
                     wilcox.test(dbARICor$V1)$p.value,
                     wilcox.test(gseaARICor$V1)$p.value,
                     wilcox.test(silARICor$V1)$p.value )
wilcoxLabelsARI <- as.data.frame(wilcoxLabelsARI)
wilcoxLabelsARI <- cbind(wilcoxLabelsARI, metric=c("ch","db","gsea","sil"))


ggplot(allARIdf, aes(x=metric, y=V1))+
  geom_boxplot()+
  xlab("Metric")+
  ylab("Correlation with ARI")+
  ggtitle("Correlation of RF Predictions with ARI")+
  geom_text(
    data    = wilcoxLabelsARI,
    mapping = aes(x = metric, y = Inf, label=round(wilcoxLabelsARI,5)),
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

ComplexHeatmap::Heatmap(t(as.matrix(corHmDf)), top_annotation = col_ha,column_title="Correlation between colData and dataset correlations",cluster_rows = FALSE)


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
  ylab("Correlation with labels")+
  ggtitle("Correlation between RF predictions and labels")+
  geom_text(
    data    = wilcoxLabels,
    mapping = aes(x = metric, y = Inf, label=round(wilcoxLabels,5)),
    hjust   = 1.05,
    vjust   = 1.5,
    color = "red"
  )+
  cowplot::theme_cowplot()
dev.off()

