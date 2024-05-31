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
})
numbers_only <- function(x) {
  sum(suppressWarnings(!is.na(as.numeric(as.character(x))))) == length(x)
}


models <- readRDS("./Documents/2022_autoMLscRNA_Cindy/data/LR_pipelinesOnly.RDS")
#filename <- "LRIntPipelinesOnly.pdf"
designMat <- readRDS("./Documents/2022_autoMLscRNA_Cindy/data/corrected_designMat_unscaled.RDS")
ari <- readRDS("./Documents/2022_autoMLscRNA_Cindy/data/ari_unscaled_cleaned.RDS")

# args <- commandArgs(trailingOnly = TRUE)
# models <- readRDS(args[[1]])
# designMat <- readRDS(args[[2]])
# ari <- readRDS(args[[3]])
# filename <- args[[4]]

#filename <- paste0("/home/campbell/cfang/automl_scrna/", filename)

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
    dplyr::rename(preds=`1`)%>%
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
    merge(nclusts, by=c("pipelines", "name")) %>%
    rename(Actual = contains("Imp"))
  return(df)
}

makeColdataDF <- function(R2df){
  dMat <- designMat %>% 
    mutate_if(numbers_only, as.numeric) %>%
    mutate_if(is.numeric, scale) %>%
    dplyr::slice(which(.$name %in% R2df$name)) %>%
    select(c("name", "detected", "sum", "ncells", "ngenes", contains("_"))) %>%
    distinct() %>%
    pivot_longer(!name, names_to="dataType") %>%
    merge(R2df, by="name") #%>%
  #mutate(PredPower = if_else(Cor >= 0.3, "R2 >= 0.1", "R2 < 0.1"))
  return(dMat)
}

createARIPlotMat <- function(df, ariMat, metric="", model=""){
  df <- df %>%
    merge(ariMat, by=c("pipelines", "name"))%>%
    dplyr::rename(ARI = value)%>%
    dplyr::rename(preds=`1`)
  ariCors <- df %>%
    group_by(name) %>%
    group_modify(~as.data.frame(cor(as.numeric(.$preds), as.numeric(.$ARI))))
  ariCors <- as.data.frame(ariCors)
  colnames(ariCors)[[2]] <- "Cor"
  return(ariCors)
}

plotPredsARI <- function(df, ariMat, metric="", npage=2){
  df <- df %>%
    merge(ariMat, by=c("pipelines", "name"))%>%
    dplyr::rename(ARI = value)
  ariCors <- df %>%
    group_by(name) %>%
    group_modify(~as.data.frame(cor(as.numeric(.$preds), as.numeric(.$ARI))))
  ariCors <- as.data.frame(ariCors)
  colnames(ariCors)[[2]] <- "V1"
  for (i in 1:npage) {
    p <- ggplot(df, aes(x=ARI, y=as.numeric(preds)))+
      geom_point(aes(colour=as.factor(res)))+
      facet_wrap(~name, scales="free")+
      labs(colour = "Clustering\nresolution")+
      xlab("Scaled ARI")+
      facet_wrap_paginate(~name, scales="free", page=i, nrow=4, ncol=2)+
      ylab(sprintf("Predicted %s from LR", metric))+
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

plotPredActuals <- function(df, metric, R2labels,npage=1){
  for (i in 1:npage){
    p <- ggplot(data=df, aes(x=as.numeric(Actual), y=as.numeric(preds)))+
      geom_point(aes(colour=as.factor(res)))+
      facet_wrap_paginate(~name, scales="free", page=i, nrow=3, ncol=3)+
      xlab(sprintf("Ground truth %s", metric))+
      ylab(sprintf("Predicted %s from LR", metric))+
      labs(colour = "Clustering\nresolution")+
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

#pdf(filename)
# ggplot(testCors, aes(x=tidytext::reorder_within(name, -Cor, metric), y=Cor, fill=metric))+
#   geom_bar(stat="identity")+
#   xlab("Dataset")+
#   ylab("Correlation between LR predictions and ground truth on test set")+
#   cowplot::theme_cowplot()+
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=6))+
#   tidytext::scale_x_reordered()+
#   facet_wrap(~ metricf, scale="free_x")+
#   theme(strip.background = element_rect(fill="white"),
#         strip.text = element_text(face="bold", size=20),
#         legend.position = "none")
# #dev.off()

#filename <- sprintf("/home/campbell/cfang/automl_scrna/results/figures/LRPlots/LRPlots.pdf")
#pdf(filename)
# Predicted vs Actual plots ###
# plotPredActuals(silTestPlot, metric="Sil (Test set)",  silTestDatasetCors, npage=3)
# plotPredActuals(silTrainPlot, metric="Sil (Train set)",  silTrainDatasetCors, npage=7)
# plotPredActuals(dbTestPlot, metric="DB (Test set)", dbTestDatasetCors, npage=3)
# plotPredActuals(dbTrainPlot, metric="DB (Train set)", dbTrainDatasetCors, npage=7)
# plotPredActuals(chTestPlot, metric="CH (Test set)",  chTestDatasetCors,npage=3)
# plotPredActuals(chTrainPlot, metric="CH (Train set)",  chTrainDatasetCors,npage=7)
# plotPredActuals(gseaTestPlot, metric="GSEA (Test set)", gseaTestDatasetCors,npage=3)
# plotPredActuals(gseaTrainPlot, metric="GSEA (Train set)",  gseaTrainDatasetCors,npage=7)
# 
# ## Predicted vs ARI plots ###
# plotPredsARI(silTestPlot, ari_tbl, metric="Sil")
# plotPredsARI(dbTestPlot, ari_tbl, metric="DB")
# plotPredsARI(chTestPlot, ari_tbl, metric="CH")
# plotPredsARI(gseaTestPlot, ari_tbl, metric="GSEA")

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
  group_modify(~ as.data.frame(cbind("Cor"=cor.test(as.numeric(.$value), as.numeric(.$Cor))$estimate, "pvalue"=cor.test(as.numeric(.$value), as.numeric(.$Cor))$p.value)))
colnames(silColdataARI)[[2]] <- "SIL"

silColdataARIPvals <- as.data.frame(cbind("dataType"=silColdataARI$dataType, "SIL"=silColdataARI$pvalue))
silColdataARI <- as.data.frame(cbind("dataType"=silColdataARI$dataType, "SIL"=silColdataARI$SIL))


dbColdataARI <- dbARI %>%
  group_by(dataType) %>%
  group_modify(~ as.data.frame(cbind("Cor"=cor.test(as.numeric(.$value), as.numeric(.$Cor))$estimate, "pvalue"=cor.test(as.numeric(.$value), as.numeric(.$Cor))$p.value)))
colnames(dbColdataARI)[[2]] <- "DB"

dbColdataARIPvals <- as.data.frame(cbind("dataType"=dbColdataARI$dataType, "DB"=dbColdataARI$pvalue))
dbColdataARI <- as.data.frame(cbind("dataType"=dbColdataARI$dataType, "DB"=dbColdataARI$DB))


chColdataARI <- chARI %>%
  group_by(dataType) %>%
  group_modify(~ as.data.frame(cbind("Cor"=cor.test(as.numeric(.$value), as.numeric(.$Cor))$estimate, "pvalue"=cor.test(as.numeric(.$value), as.numeric(.$Cor))$p.value)))
colnames(chColdataARI)[[2]] <- "CH"

chColdataARIPvals <- as.data.frame(cbind("dataType"=chColdataARI$dataType, "CH"=chColdataARI$pvalue))
chColdataARI <- as.data.frame(cbind("dataType"=chColdataARI$dataType, "CH"=chColdataARI$CH))


gseaColdataARI <- gseaARI %>%
  group_by(dataType) %>%
  group_modify(~ as.data.frame(cbind("Cor"=cor.test(as.numeric(.$value), as.numeric(.$Cor))$estimate, "pvalue"=cor.test(as.numeric(.$value), as.numeric(.$Cor))$p.value)))
colnames(gseaColdataARI)[[2]] <- "GSEA"

gseaColdataARIPvals <- as.data.frame(cbind("dataType"=gseaColdataARI$dataType, "GSEA"=gseaColdataARI$pvalue))
gseaColdataARI <- as.data.frame(cbind("dataType"=gseaColdataARI$dataType, "GSEA"=gseaColdataARI$GSEA))


corHmARI <- merge(silColdataARI, dbColdataARI, by="dataType")
corHmARI <- merge(corHmARI, chColdataARI, by="dataType")
corHmARI<- merge(corHmARI, gseaColdataARI, by="dataType")
corHmARInames <- corHmARI$dataType
corHmARI$dataType <- NULL
corHmARI <- sapply(corHmARI, as.numeric)
col_ha <- rowAnnotation("Average" = anno_barplot(rowMeans(corHmARI)))
rownames(corHmARI) <- corHmARInames

corHmARIpvals <- merge(silColdataARIPvals, dbColdataARIPvals, by="dataType")
corHmARIpvals <- merge(corHmARIpvals, chColdataARIPvals, by="dataType")
corHmARIpvals<- merge(corHmARIpvals, gseaColdataARIPvals, by="dataType")
rownames(corHmARIpvals) <- corHmARIpvals$dataType
corHmARIpvals$dataType <- NULL

corHmARIpvalsSig <- corHmARIpvals %>%
  filter(SIL < 0.05 | DB < 0.05 | CH < 0.05 | GSEA < 0.05)

print(dim(corHmARIpvalsSig))
print(corHmARI)
print(dim(corHmARI))
corHmARISig <- corHmARI[rownames(corHmARIpvalsSig),]
print(dim(corHmARISig))
corHmARISig<- sapply(as.data.frame(corHmARISig), as.numeric)
col_ha <- rowAnnotation("Average" = anno_barplot(rowMeans(corHmARISig)))
print(corHmARISig)
print(dim(corHmARISig))
print(dim(corHmARIpvalsSig))
#rownames(corHmARISig) <- names(corHmARIpvalsSig)

roundcorHmpvalsARISig <- as.data.frame(sapply(corHmARIpvalsSig, as.numeric))
roundcorHmpvalsARISig<- as.data.frame(sapply(roundcorHmpvalsARISig, signif,3))
roundcorHmpvalsARISig <- as.data.frame(sapply(roundcorHmpvalsARISig, as.character))


ariHm <- ComplexHeatmap::Heatmap((as.matrix(corHmARISig)), 
                                 name="Correlation of dataset feature with dataset-specific predictive performance (ARI)", 
                                 left_annotation = col_ha, 
                                 cluster_columns = FALSE,
                                 heatmap_legend_param = list(
                                   legend_direction = "horizontal"))
### ColData ARI Cor boxplot ###
wilcoxLabelsARI <- c(wilcox.test(chARICor$Cor, alternative="greater")$p.value,
                     wilcox.test(dbARICor$Cor, alternative="greater")$p.value,
                     wilcox.test(gseaARICor$Cor, alternative="greater")$p.value,
                     wilcox.test(silARICor$Cor, alternative="greater")$p.value )
wilcoxLabelsARI <- as.data.frame(wilcoxLabelsARI)
wilcoxLabelsARI <- cbind(wilcoxLabelsARI, metric=c("CH","DB","GSEA","SIL"))
allARIdf <- mutate(allARIdf, metric=toupper(metric))
wilcoxLabelsARI$wilcoxLabelsARI <- p.adjust(wilcoxLabelsARI$wilcoxLabelsARI, method="BH")

pdf("./Documents/2022_autoMLscRNA_Cindy/PaperFigures/predictionsARICorBoxplotLRpipelinesOnly_padj.pdf", width=5.25, height=6)
allARIdf$metricf <- factor(allARIdf$metric, levels=c("CH", "DB", "SIL", "GSEA"))
ggplot(as.data.frame(allARIdf), aes(x=metricf, y=Cor, fill=metric))+
  geom_boxplot()+
  xlab("Metric")+
  ylab("Correlation of LR predictions with ARI")+
  geom_text(
    data    = wilcoxLabelsARI,
    mapping = aes(x = metric, y = Inf, label= paste("p =", signif(wilcoxLabelsARI,3))),
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

# pdf("PaperFigures/ariCorBarLRIntOrig.pdf")
# allARIdf$metricf <- factor(allARIdf$metric, levels=c("CH", "DB", "SIL", "GSEA"))
# ggplot(allARIdf, aes(x=tidytext::reorder_within(name, -Cor, metric), y=Cor, fill=metric))+
#   geom_bar(stat="identity")+
#   facet_wrap(~metric, scale="free")+
#   xlab("Dataset")+
#   ylab("Correlation between LR predictions and ARI on test set")+
#   cowplot::theme_cowplot()+
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=6))+
#   tidytext::scale_x_reordered()+
#   facet_wrap(~ metricf, scale="free_x")+
#   theme(strip.background = element_rect(fill="white"),
#         strip.text = element_text(face="bold", size=20),
#         legend.position = "none")
# 
# dev.off()
### ColData Correlation Heatmap ###
silCol <- makeColdataDF(silTestDatasetCors)
dbCol <- makeColdataDF(dbTestDatasetCors)
chCol <- makeColdataDF(chTestDatasetCors)
gseaCol <- makeColdataDF(gseaTestDatasetCors)

print(silCol)
print(silTestDatasetCors)

alldf <- dplyr::bind_rows(list(sil=silCol, db=dbCol, ch=chCol, gsea=gseaCol), .id = 'metric')
alldf <- alldf #%>%
#mutate(xLabel = paste(metric, PredPower, sep=" "))

silColdataCor <- silCol %>%
  group_by(dataType) %>%
  group_modify(~ as.data.frame(cbind("cor"=cor.test(as.numeric(.$value), as.numeric(.$Cor))$estimate, "pvalue"=cor.test(as.numeric(.$value), as.numeric(.$Cor))$p.value)))
colnames(silColdataCor)[[2]] <- "SIL"

silColdataCorPvals <- as.data.frame(cbind("dataType"=silColdataCor$dataType, "SIL"=silColdataCor$pvalue))
silColdataCor <- as.data.frame(cbind("dataType"=silColdataCor$dataType, "SIL"=silColdataCor$SIL))


dbColdataCor <- dbCol %>%
  group_by(dataType) %>%
  group_modify(~ as.data.frame(cbind("cor"=cor.test(as.numeric(.$value), as.numeric(.$Cor))$estimate, "pvalue"=cor.test(as.numeric(.$value), as.numeric(.$Cor))$p.value)))
colnames(dbColdataCor)[[2]] <- "DB"

dbColdataCorPvals <- as.data.frame(cbind("dataType"=dbColdataCor$dataType, "DB"=dbColdataCor$pvalue))
dbColdataCor <- as.data.frame(cbind("dataType"=dbColdataCor$dataType, "DB"=dbColdataCor$DB))


chColdataCor <- chCol %>%
  group_by(dataType) %>%
  group_modify(~ as.data.frame(cbind("cor"=cor.test(as.numeric(.$value), as.numeric(.$Cor))$estimate, "pvalue"=cor.test(as.numeric(.$value), as.numeric(.$Cor))$p.value)))
colnames(chColdataCor)[[2]] <- "CH"

chColdataCorPvals <- as.data.frame(cbind("dataType"=chColdataCor$dataType, "CH"=chColdataCor$pvalue))
chColdataCor <- as.data.frame(cbind("dataType"=chColdataCor$dataType, "CH"=chColdataCor$CH))


gseaColdataCor <- gseaCol %>%
  group_by(dataType) %>%
  group_modify(~ as.data.frame(cbind("cor"=cor.test(as.numeric(.$value), as.numeric(.$Cor))$estimate, "pvalue"=cor.test(as.numeric(.$value), as.numeric(.$Cor))$p.value)))
colnames(gseaColdataCor)[[2]] <- "GSEA"

gseaColdataCorPvals <- as.data.frame(cbind("dataType"=gseaColdataCor$dataType, "GSEA"=gseaColdataCor$pvalue))
gseaColdataCor <- as.data.frame(cbind("dataType"=gseaColdataCor$dataType, "GSEA"=gseaColdataCor$GSEA))

print(silColdataCor)
print(dbColdataCor)
corHmDf <- merge(silColdataCor, dbColdataCor, by="dataType")
print(corHmDf)
corHmDf <- merge(corHmDf, chColdataCor, by="dataType")
print(corHmDf)
corHmDf <- merge(corHmDf, gseaColdataCor, by="dataType")
print(corHmDf)
corHmDfNames <- corHmDf$dataType
corHmDf$dataType <- NULL
corHmDf <- sapply(as.data.frame(corHmDf), as.numeric)
col_ha <- rowAnnotation("Average" = anno_barplot(rowMeans(corHmDf)))
rownames(corHmDf) <- corHmDfNames

corHmpvals <- merge(silColdataCorPvals, dbColdataCorPvals, by="dataType")
corHmpvals <- merge(corHmpvals, chColdataCorPvals, by="dataType")
corHmpvals <- merge(corHmpvals, gseaColdataCorPvals, by="dataType")
rownames(corHmpvals) <- corHmpvals$dataType
corHmpvals$dataType <- NULL

corHmpvalsSig <- corHmpvals %>%
  filter(SIL < 0.05 | DB < 0.05 | CH < 0.05 | GSEA < 0.05 )

corHmDfSig <- corHmDf[rownames(corHmpvalsSig),]
corHmDfSig <- sapply(as.data.frame(corHmDfSig), as.numeric)
col_ha <- rowAnnotation("Average" = anno_barplot(rowMeans(corHmDfSig)))
rownames(corHmDfSig) <- rownames(corHmpvalsSig)

roundcorHmpvalsSig <- as.data.frame(sapply(corHmpvalsSig, as.numeric))
roundcorHmpvalsSig <- as.data.frame(sapply(roundcorHmpvalsSig, signif,3))
roundcorHmpvalsSig <- as.data.frame(sapply(roundcorHmpvalsSig, as.character))

# corHm <- ComplexHeatmap::Heatmap(as.matrix(corHmDfSig), 
#                                  name="Correlation of dataset features with dataset-specific predictive performance", 
#                                  left_annotation = col_ha,
#                                  cluster_columns = FALSE,
#                                  heatmap_legend_param = list(
#                                    legend_direction = "horizontal"))
# 
# ht_opt$heatmap_row_names_gp = gpar(fontsize = 10)
# ht_opt$heatmap_column_names_gp = gpar(fontsize = 10)
# 
# pdf("PaperFigures/datasetFeaturePerformanceHeatmapsLRIntOrig.pdf")
# draw(corHm,
#      row_title="Feature", 
#      column_title="Metric", 
#      column_title_side="bottom",
#      heatmap_legend_side="top")
# 
# draw(ariHm,
#      row_title="Feature", 
#      column_title="Metric", 
#      column_title_side="bottom",
#      heatmap_legend_side="top")
# dev.off()

### ColData Correlation Boxplot ###
corbpdf <- merge(silTestDatasetCors, dbTestDatasetCors, by="name", suffix=c("Sil","Db"))
chgsea <- merge(chTestDatasetCors, gseaTestDatasetCors, by="name", suffix=c("Ch", "Gsea"))
corbpdf <- merge(corbpdf, chgsea, by="name")

corbpdf <- corbpdf %>%
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

corbpdf$metricf <- factor(corbpdf$metric, levels=c("CH", "DB", "SIL", "GSEA"))

pdf("./Documents/2022_autoMLscRNA_Cindy/PaperFigures/predictionMetricCorBoxplotLRpipelinesOnly_padj.pdf", width=6.25, height=6.1)
ggplot(corbpdf, aes(x=metricf, y=value, fill=metric))+
  geom_boxplot()+
  xlab("Metric")+
  ylab("Correlation of LR predictions\nwith corrected observed metric")+
  geom_text(
    data    = wilcoxLabels,
    mapping = aes(x = metric, y = Inf, label=paste("p =", signif(wilcoxLabels,3))),
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

corbpdfMeans <- corbpdf %>%
  ungroup()%>%
  group_by(metric) %>%
  group_modify(~as.data.frame(cbind(median=median(.$value))))



