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
  library(tidytext)
})
numbers_only <- function(x) {
  sum(suppressWarnings(!is.na(as.numeric(as.character(x))))) == length(x)
}


# models <- readRDS("RFPreds_NumericParams.RDS")
# filename <- "allFeaturesRFPlots.pdf"
# designMat <- readRDS("corrected_designMat_unscaled.RDS")
# ari <- readRDS("ari_unscaled_cleaned.RDS")
# RFType <- "Numeric"

args <- commandArgs(trailingOnly=TRUE)

models <- readRDS(args[[1]])
designMat <- readRDS(args[[2]])
ari <- readRDS(args[[3]])

filename <- args[[4]]

args <- commandArgs(trailingOnly=TRUE)
print(args)

models <- readRDS(args[[1]])
designMat <- readRDS(args[[2]])
ari <- readRDS(args[[3]])

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
    group_modify(~ as.data.frame(cbind("Cor"=cor.test(as.numeric(.$preds), as.numeric(.$Actual))$estimate, "p.value"=cor.test(as.numeric(.$preds), as.numeric(.$Actual))$p.value)))
  return(datasetCors)
}

createARIPlotMat <- function(df, ariMat, metric="", model=""){
  df <- df %>%
    merge(ariMat, by=c("pipelines", "name"))%>%
    dplyr::rename(ARI = value)
  ariCors <- df %>%
    group_by(name) %>%
    group_modify(~as.data.frame(cbind("Cor"=cor.test(as.numeric(.$preds), as.numeric(.$ARI))$estimate,"p.value"=cor.test(as.numeric(.$preds), as.numeric(.$ARI))$p.value)))
  ariCors <- as.data.frame(ariCors)
  colnames(ariCors)[[2]] <- "Cor"
  return(ariCors)
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
    group_modify(~as.data.frame(cbind("Cor"=cor.test(as.numeric(.$preds), as.numeric(.$ARI))$estimate,"p.value"=cor.test(as.numeric(.$preds), as.numeric(.$ARI))$p.value)))
  ariCors <- as.data.frame(ariCors)
  colnames(ariCors)[[2]] <- "Cor"
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

testCorsSignif <- testCors
pdf(filename)

ggplot(testCors, aes(x=tidytext::reorder_within(name, -Cor, metric), y=Cor, fill=metric))+
  geom_bar(stat="identity")+
  xlab("Dataset")+
  ylab("Correlation between RF predictions\nand ground truth on test set (all features)")+
  cowplot::theme_cowplot()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=6))+
  tidytext::scale_x_reordered()+
  facet_wrap(~ metricf, scale="free_x")+
  theme(strip.background = element_rect(fill="white"),
        strip.text = element_text(face="bold", size=20),
        legend.position = "none")+
  geom_text(
    data    = testCorsSignif,
    mapping = aes(x=tidytext::reorder_within(name, -Cor, metric), y = 0, label=paste("p =", signif(testCorsSignif$p.value,3))),
    hjust   = 0,
    vjust   = 0.1,
    color = "black",
    size=3,
    angle=90
  )

silParams = sprintf("RF (mtry=%s, ntree=%s)", silRF$mtry, silRF$ntree)
dbParams = sprintf("RF (mtry=%s, ntree=%s)", dbRF$mtry, dbRF$ntree)
chParams = sprintf("RF (mtry=%s, ntree=%s)", chRF$mtry, chRF$ntree)
gseaParams = sprintf("RF (mtry=%s, ntree=%s)", gseaRF$mtry, gseaRF$ntree)


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
corHmARI <- sapply(as.data.frame(corHmARI), as.numeric)
col_ha <- rowAnnotation("Average" = anno_barplot(rowMeans(corHmARI)))
print(dim(corHmARI))
print(length(corHmARInames))
rownames(corHmARI) <- corHmARInames

corHmARIpvals <- merge(silColdataARIPvals, dbColdataARIPvals, by="dataType")
corHmARIpvals <- merge(corHmARIpvals, chColdataARIPvals, by="dataType")
corHmARIpvals<- merge(corHmARIpvals, gseaColdataARIPvals, by="dataType")
rownames(corHmARIpvals) <- corHmARIpvals$dataType
corHmARIpvals$dataType <- NULL

corHmARIpvalsSig <- corHmARIpvals %>%
  filter(SIL < 0.05 | DB < 0.05 | CH < 0.05 | GSEA < 0.05)

corHmARISig <- corHmARI[rownames(corHmARIpvalsSig),]
corHmARISig<- sapply(as.data.frame(corHmARISig), as.numeric)
col_ha <- rowAnnotation("Average" = anno_barplot(rowMeans(corHmARISig)))
rownames(corHmARISig) <- rownames(corHmARIpvalsSig)

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

allARIdf$metricf <- factor(allARIdf$metric, levels=c("CH", "DB", "SIL", "GSEA"))

ggplot(allARIdf, aes(x=metricf, y=Cor, fill=metric))+
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

ggplot(allARIdf, aes(x=tidytext::reorder_within(name, -Cor, metric), y=Cor, fill=metric))+
  geom_bar(stat="identity")+
  #facet_wrap(~metric, scale="free")+
  xlab("Dataset")+
  ylab("Correlation between RF predictions\nand ARI on test set (all features)")+
  cowplot::theme_cowplot()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=6))+
  tidytext::scale_x_reordered()+
  facet_wrap(~ metricf, scale="free_x")+
  theme(strip.background = element_rect(fill="white"),
        strip.text = element_text(face="bold", size=20),
        legend.position = "none")#+
  # geom_text(
  #   data    = ARIsignif,
  #   mapping = aes(x=tidytext::reorder_within(name, -Cor, metric), y = 0, label=paste("p =", signif(ARIsignif$p.value,3))),
  #   hjust   = 0,
  #   vjust   = 0.1,
  #   color = "black",
  #   size=3,
  #   angle=90
  # )

### ColData Correlation Heatmap ###
silCol <- makeColdataDF(silTestDatasetCors)
dbCol <- makeColdataDF(dbTestDatasetCors)
chCol <- makeColdataDF(chTestDatasetCors)
gseaCol <- makeColdataDF(gseaTestDatasetCors)

alldf <- dplyr::bind_rows(list(sil=silCol, db=dbCol, ch=chCol, gsea=gseaCol), .id = 'metric')

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

corHmDf <- merge(silColdataCor, dbColdataCor, by="dataType")
corHmDf <- merge(corHmDf, chColdataCor, by="dataType")
corHmDf <- merge(corHmDf, gseaColdataCor, by="dataType")
corHmDfNames <- corHmDf$dataType
corHmDf$dataType <- NULL
corHmDf <- sapply(as.data.frame(corHmDf), as.numeric)
col_ha <- rowAnnotation("Average" = anno_barplot(rowMeans(corHmDf)))
rownames(corHmDf) <- corHmDfNames

corHmpvals <- merge(silColdataCorPvals, dbColdataCorPvals, by="dataType")
corHmpvals <- merge(corHmpvals, chColdataCorPvals, by="dataType")
corHmpvals <- merge(corHmpvals, gseaColdataCorPvals, by="dataType")
rownames(corHmpvals) <- corHmpvals$dataType
print(corHmpvals)
corHmpvals$dataType <- NULL

corHmpvalsSig <- corHmpvals %>%
  filter(SIL < 0.05 | DB < 0.05 | CH < 0.05 | GSEA < 0.05 )

corHmDfSig <- corHmDf[rownames(corHmpvalsSig),]
corHmDfSig <- sapply(as.data.frame(corHmDfSig), as.numeric)
col_ha <- rowAnnotation("Average" = anno_barplot(rowMeans(corHmDfSig)))
print(corHmDfSig)
print(corHmpvalsSig)
rownames(corHmDfSig) <- rownames(corHmpvalsSig)


roundcorHmpvalsSig <- as.data.frame(sapply(corHmpvalsSig, as.numeric))
roundcorHmpvalsSig <- as.data.frame(sapply(roundcorHmpvalsSig, signif,3))
roundcorHmpvalsSig <- as.data.frame(sapply(roundcorHmpvalsSig, as.character))

corHm <- ComplexHeatmap::Heatmap(as.matrix(corHmDfSig), 
                                 name="Correlation of dataset features with dataset-specific predictive performance", 
                                 left_annotation = col_ha,
                                 cluster_columns = FALSE,
                                 heatmap_legend_param = list(
                                   legend_direction = "horizontal"))

# cell_fun = function(j, i, x, y, width, height, fill) {
#   grid.text(sprintf(roundcorHmpvalsSig[i, j]), x, y, gp = gpar(fontsize = 10))
# }

ht_opt$heatmap_row_names_gp = gpar(fontsize = 10)
ht_opt$heatmap_column_names_gp = gpar(fontsize = 10)


draw(corHm,
     row_title="Feature", 
     column_title="Metric", 
     column_title_side="bottom",
     heatmap_legend_side="top")

draw(ariHm,
     row_title="Feature", 
     column_title="Metric", 
     column_title_side="bottom",
     heatmap_legend_side="top")
#dev.off()

### ColData Correlation Boxplot ###
corbpdf <- merge(silTestDatasetCors, dbTestDatasetCors, by="name", suffix=c("Sil","Db"))
chgsea <- merge(chTestDatasetCors, gseaTestDatasetCors, by="name", suffix=c("Ch", "Gsea"))
corbpdf <- merge(corbpdf, chgsea, by="name")

corbpdf <- corbpdf %>%
  as_tibble(corbpdf) %>%
  pivot_longer(cols=starts_with("Cor"), names_to="metric")

wilcoxLabels <- c(wilcox.test(chTestDatasetCors$Cor, alternative="greater")$p.value,
                  wilcox.test(dbTestDatasetCors$Cor, alternative="greater")$p.value,
                  wilcox.test(gseaTestDatasetCors$Cor, alternative="greater")$p.value,
                  wilcox.test(silTestDatasetCors$Cor, alternative="greater")$p.value )
wilcoxLabels <- as.data.frame(wilcoxLabels)
wilcoxLabels <- cbind(wilcoxLabels, metric=c("CH","DB","GSEA","SIL"))
corbpdf$metric <- sub("^Cor", "", corbpdf$metric)

corbpdf <- corbpdf %>%
  mutate(metric=toupper(metric))


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

corbpdfMeans <- corbpdf %>%
  ungroup()%>%
  group_by(metric) %>%
  group_modify(~as.data.frame(cbind(mean=mean(.$value))))
