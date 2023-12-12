suppressPackageStartupMessages({
  library(tidyverse) # CRAN
  library(here) # CRAN
  library(DT) # CRAN
  library(here) # CRAN
  library(Matrix) # CRAN
  library(purrr)
  library(ggplot2)
  #library(Metrics)
  library(gridExtra)
  library(matrixStats)
  library(stringr)
  library(stringi)
  library(ggforce)
  library(cowplot)
  library(ComplexHeatmap)
  library(tidytext)
})
numbers_only <- function(x) {
  sum(suppressWarnings(!is.na(as.numeric(as.character(x))))) == length(x)
}


models <- readRDS("./Documents/2022_autoMLscRNA_Cindy/data/pipelinesOnlyRFPreds_Params.RDS")
designMat <- readRDS("./Documents/2022_autoMLscRNA_Cindy/data/corrected_designMat_unscaled.RDS")
ari <- readRDS("./Documents/2022_autoMLscRNA_Cindy/data/ari_unscaled_cleaned.RDS")
#RFType <- "Numeric"

# args <- commandArgs(trailingOnly=TRUE)
# print(args)
# 
# models <- readRDS(args[[1]])
# designMat <- readRDS(args[[2]])
# ari <- readRDS(args[[3]])
# filename <- args[[4]]
#RFType <- args[[4]]

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
    select(c("name", "detected", "sum", "ncells", "ngenes", contains("_"))) %>%
    distinct() %>%
    pivot_longer(!name, names_to="dataType") %>%
    merge(R2df, by="name")
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

plotPredsARI <- function(df, ariMat, metric="", model="", npage=2){
  df <- df %>%
    merge(ariMat, by=c("pipelines", "name"))%>%
    dplyr::rename(ARI = value)
  print(sprintf("Predicted %s from %s", metric, model))
  ariCors <- df %>%
    group_by(name) %>%
    group_modify(~as.data.frame(cor(.$preds, .$ARI)))
  ariCors <- as.data.frame(ariCors)
  for (i in 1:npage) {
    p <- ggplot(df, aes(x=ARI, y=preds))+
      geom_point(aes(colour=as.factor(res)))+
      labs(colour = "Clustering\nresolution")+
      facet_wrap(~name, scales="free")+
      xlab("Scaled ARI")+
      facet_wrap_paginate(~name, scales="free", page=i, nrow=3, ncol=2)+
      ylab(sprintf("Predicted %s from %s", metric, model))+
      geom_text(
        data    = ariCors,
        mapping = aes(x = Inf, y = Inf, label = paste("Cor:", signif(V1, digits=3))),
        hjust   = 1.05,
        vjust   = 1.5
      )+ 
      cowplot::theme_cowplot()+
      theme(strip.background = element_rect(fill="white"),
            strip.text = element_text(face="bold"))
    print(p)
  }
  
}

plotPredActuals <- function(df, metric, model, R2labels,npage=1){
  print(sprintf("Predicted %s from %s", metric, model))
  for (i in 1:npage){
    p <- ggplot(data=df, aes(x=as.numeric(Actual), y=as.numeric(preds)))+
      geom_point(aes(colour=as.factor(res)))+
      labs(colour = "Clustering\nresolution")+
      facet_wrap_paginate(~name, scales="free", page=i, nrow=3, ncol=3)+
      xlab(sprintf("Ground truth %s", metric))+
      ylab(sprintf("Predicted %s from %s", metric, model))+
      geom_text(
        data    = R2labels,
        mapping = aes(x = Inf, y = Inf, label = paste("Cor:", signif(Cor, digits=3))),
        hjust   = 1.05,
        vjust   = 1.5,
        color = "black"
      )+ 
      cowplot::theme_cowplot()+
      theme(strip.background = element_rect(fill="white"),
            strip.text = element_text(face="bold"))
    print(p)
  }
  
}

ari_tbl <- ari %>%
  as_tibble() %>%
  mutate_if(is.numeric, scale) %>%
  pivot_longer(starts_with("E"))
colnames(ari_tbl)[[3]] <- "value"

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

testCors$metricf <- factor(testCors$metric, levels=c("CH", "DB", "SIL", "GSEA"))
#pdf("PaperFigures/RFOriginalResultsPlots_pipelinesOnly.pdf")
ggplot(testCors, aes(x=tidytext::reorder_within(name, -Cor, metric), y=Cor))+
  geom_bar(stat="identity")+
  xlab("Dataset")+
  ylab("Correlation between RF predictions and ground truth on test set")+
  cowplot::theme_cowplot()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=6))+
  tidytext::scale_x_reordered()+
  facet_wrap(~ metricf, scale="free_x")+
  theme(strip.background = element_rect(fill="white"),
        strip.text = element_text(face="bold", size=20))
#dev.off()

silParams = sprintf("RF (mtry=%s, ntree=%s)", silRF$mtry, silRF$ntree)
dbParams = sprintf("RF (mtry=%s, ntree=%s)", dbRF$mtry, dbRF$ntree)
chParams = sprintf("RF (mtry=%s, ntree=%s)", chRF$mtry, chRF$ntree)
gseaParams = sprintf("RF (mtry=%s, ntree=%s)", gseaRF$mtry, gseaRF$ntree)

filename <-"PaperFigures/RFPlots_NumericParams_plots_pipelinesOnly.pdf"
#pdf(filename)
# Predicted vs Actual plots ###
plotPredActuals(silTestPlot, metric="Sil (Test set)", model=silParams, silTestDatasetCors, npage=3)
plotPredActuals(silTrainPlot, metric="Sil (Train set)", model=silParams, silTrainDatasetCors, npage=7)
plotPredActuals(dbTestPlot, metric="DB (Test set)", model=dbParams, dbTestDatasetCors, npage=3)
plotPredActuals(dbTrainPlot, metric="DB (Train set)", model=dbParams, dbTrainDatasetCors, npage=7)
plotPredActuals(chTestPlot, metric="CH (Test set)", model=chParams, chTestDatasetCors,npage=3)
plotPredActuals(chTrainPlot, metric="CH (Train set)", model=chParams, chTrainDatasetCors,npage=7)
plotPredActuals(gseaTestPlot, metric="GSEA (Test set)", model=gseaParams, gseaTestDatasetCors,npage=3)
plotPredActuals(gseaTrainPlot, metric="GSEA (Train set)", model=gseaParams, gseaTrainDatasetCors,npage=7)

## Predicted vs ARI plots ###
plotPredsARI(silTestPlot, ari_tbl, metric="Sil", model=silParams)
plotPredsARI(dbTestPlot, ari_tbl, metric="DB", model=dbParams)
plotPredsARI(chTestPlot, ari_tbl, metric="CH", model=chParams)
plotPredsARI(gseaTestPlot, ari_tbl, metric="GSEA", model=gseaParams)

#dev.off()
### ColData ARI Correlation heatmap ###

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
colnames(silColdataARI)[[2]] <- "SIL"

dbColdataARI <- dbARI %>%
  group_by(dataType) %>%
  group_modify(~ as.data.frame(cor(as.numeric(.$value), as.numeric(.$V1))))
colnames(dbColdataARI)[[2]] <- "DB"

chColdataARI <- chARI %>%
  group_by(dataType) %>%
  group_modify(~ as.data.frame(cor(as.numeric(.$value), as.numeric(.$V1))))
colnames(chColdataARI)[[2]] <- "CH"

gseaColdataARI <- gseaARI %>%
  group_by(dataType) %>%
  group_modify(~ as.data.frame(cor(as.numeric(.$value), as.numeric(.$V1))))
colnames(gseaColdataARI)[[2]] <- "GSEA"

corHmARI <- merge(silColdataARI, dbColdataARI, by="dataType")
corHmARI <- merge(corHmARI, chColdataARI, by="dataType")
corHmARI<- merge(corHmARI, gseaColdataARI, by="dataType")
rownames(corHmARI) <- corHmARI$dataType
corHmARI$dataType <- NULL
col_ha <- rowAnnotation("Average" = anno_barplot(rowMeans(corHmARI)))

ariHm <- ComplexHeatmap::Heatmap((as.matrix(corHmARI)), name="Correlation of\ndataset feature\nwith dataset-specific\npredictive performance\n(ARI)", left_annotation = col_ha, cluster_columns = FALSE)

### ColData ARI Cor boxplot ###
wilcoxLabelsARI <- c(wilcox.test(chARICor$V1, alternative="greater")$p.value,
                     wilcox.test(dbARICor$V1, alternative="greater")$p.value,
                     wilcox.test(gseaARICor$V1, alternative="greater")$p.value,
                     wilcox.test(silARICor$V1, alternative="greater")$p.value )
wilcoxLabelsARI <- as.data.frame(wilcoxLabelsARI)
wilcoxLabelsARI <- cbind(wilcoxLabelsARI, metric=c("CH","DB","GSEA","SIL"))
allARIdf <- mutate(allARIdf, metric=toupper(metric))
wilcoxLabelsARI$wilcoxLabelsARI <- p.adjust(wilcoxLabelsARI$wilcoxLabelsARI, method="BH")

allARIdf$metricf <- factor(allARIdf$metric, levels=c("CH", "DB", "SIL", "GSEA"))

pdf("./Documents/2022_autoMLscRNA_Cindy/PaperFigures/RFOrig_predictionsARICorBoxplot_pipelinesOnly_padj.pdf", width=5, height=5.85)
ggplot(allARIdf, aes(x=metricf, y=V1, fill=metric))+
  geom_boxplot()+
  xlab("Metric")+
  ylab("Correlation of RF predictions with ARI")+
  geom_text(
    data    = wilcoxLabelsARI,
    mapping = aes(x = metric, y = Inf, label=paste("p =", signif(wilcoxLabelsARI,3))),
    hjust   = 0.65,
    vjust   = 1.5,
    color = "black",
    size=5
  )+
  scale_fill_brewer(palette="Set3")+
  cowplot::theme_cowplot()+
  theme(axis.title.x=element_text(size=18),
        axis.title.y=element_text(size=18),
        legend.position = "none")
dev.off()

#pdf("PaperFigures/RFOrig_ariCorBar_pipelinesOnly.pdf")
ggplot(allARIdf, aes(x=tidytext::reorder_within(name, -V1, metric), y=V1))+
  geom_bar(stat="identity")+
  facet_wrap(~metric, scale="free")+
  xlab("Dataset")+
  ylab("Correlation between RF predictions and ARI on test set")+
  cowplot::theme_cowplot()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=6))+
  tidytext::scale_x_reordered()+
  facet_wrap(~ metricf, scale="free_x")+
  theme(strip.background = element_rect(fill="white"),
        strip.text = element_text(face="bold", size=20))
#dev.off()

### ColData Correlation Heatmap ###
silCol <- makeColdataDF(silTestDatasetCors)
dbCol <- makeColdataDF(dbTestDatasetCors)
chCol <- makeColdataDF(chTestDatasetCors)
gseaCol <- makeColdataDF(gseaTestDatasetCors)

alldf <- dplyr::bind_rows(list(sil=silCol, db=dbCol, ch=chCol, gsea=gseaCol), .id = 'metric')

silColdataCor <- silCol %>%
  group_by(dataType) %>%
  group_modify(~ as.data.frame(cor(as.numeric(.$value), as.numeric(.$Cor))))
colnames(silColdataCor)[[2]] <- "SIL"

dbColdataCor <- dbCol %>%
  group_by(dataType) %>%
  group_modify(~ as.data.frame(cor(as.numeric(.$value), as.numeric(.$Cor))))
colnames(dbColdataCor)[[2]] <- "DB"

chColdataCor <- chCol %>%
  group_by(dataType) %>%
  group_modify(~ as.data.frame(cor(as.numeric(.$value), as.numeric(.$Cor))))
colnames(chColdataCor)[[2]] <- "CH"

gseaColdataCor <- gseaCol %>%
  group_by(dataType) %>%
  group_modify(~ as.data.frame(cor(as.numeric(.$value), as.numeric(.$Cor))))
colnames(gseaColdataCor)[[2]] <- "GSEA"

corHmDf <- merge(silColdataCor, dbColdataCor, by="dataType")
corHmDf <- merge(corHmDf, chColdataCor, by="dataType")
corHmDf <- merge(corHmDf, gseaColdataCor, by="dataType")
rownames(corHmDf) <- corHmDf$dataType
corHmDf$dataType <- NULL
col_ha <- rowAnnotation("Average" = anno_barplot(rowMeans(corHmDf)))

corHm <- ComplexHeatmap::Heatmap((as.matrix(corHmDf)), name="Correlation of\ndataset features\nwith dataset-specific\npredictive performance", left_annotation = col_ha,cluster_columns = FALSE)

ht_opt$heatmap_row_names_gp = gpar(fontsize = 16)
ht_opt$heatmap_column_names_gp = gpar(fontsize = 16)

# pdf("PaperFigures/RFOrig_datasetFeaturePerformanceHeatmaps_pipelinesOnly.pdf")
# draw(corHm,
#      row_title="Feature", 
#      column_title="Metric", 
#      column_title_side="bottom")
# 
# draw(ariHm,
#      row_title="Feature", 
#      column_title="Metric", 
#      column_title_side="bottom")
# dev.off()

### ColData Correlation Boxplot ###
corbpdf <- merge(silTestDatasetCors, dbTestDatasetCors, by="name", suffix=c("Sil","Db"))
chgsea <- merge(chTestDatasetCors, gseaTestDatasetCors, by="name", suffix=c("Ch", "Gsea"))
corbpdf <- merge(corbpdf, chgsea, by="name")

corbpdf <- corbpdf %>%
  #as_tibble(corbpdf) %>%
  pivot_longer(cols=starts_with("Cor"), names_to="metric")

wilcoxLabels <- c(wilcox.test(chTestDatasetCors$Cor, alternative="greater")$p.value,
                  wilcox.test(dbTestDatasetCors$Cor, alternative="greater")$p.value,
                  wilcox.test(gseaTestDatasetCors$Cor, alternative="greater")$p.value,
                  wilcox.test(silTestDatasetCors$Cor, alternative="greater")$p.value )
wilcoxLabels <- as.data.frame(wilcoxLabels)
wilcoxLabels <- cbind(wilcoxLabels, metric=c("CH","DB","GSEA","SIL"))
corbpdf$metric <- sub("^Cor", "", corbpdf$metric)
wilcoxLabels$wilcoxLabels <- p.adjust(wilcoxLabels$wilcoxLabels, method="BH")

corbpdf <- corbpdf %>%
  mutate(metric=toupper(metric))

pdf("./Documents/2022_autoMLscRNA_Cindy/PaperFigures/RFOrig_predictionMetricCorBoxplot_pipelinesOnly_padj.pdf", width=6.25, height=6.1)
corbpdf$metricf <- factor(corbpdf$metric, levels=c("CH", "DB", "SIL", "GSEA"))
ggplot(corbpdf, aes(x=metricf, y=value, fill=metric))+
  geom_boxplot()+
  xlab("Metric")+
  ylab("Correlation of RF predictions\nwith corrected observed metric")+
  geom_text(
    data    = wilcoxLabels,
    mapping = aes(x = metric, y = Inf, label=paste("p =",signif(wilcoxLabels,3))),
    hjust   = 0.65,
    vjust   = 1.5,
    color = "black",
    size=5
  )+
  scale_fill_brewer(palette="Set3")+
  cowplot::theme_cowplot()+
  theme(axis.title.x=element_text(size=18),
        axis.title.y=element_text(size=18),
        legend.position = "none")
dev.off()

