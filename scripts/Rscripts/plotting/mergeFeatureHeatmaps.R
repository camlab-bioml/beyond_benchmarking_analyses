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


models <- readRDS("LR_int_allTerms.RDS")
pipeLRModels <- readRDS("LR_pipelinesOnly.RDS")

designMat <- readRDS("corrected_designMat_unscaled.RDS")
ari <- readRDS("ari_unscaled_cleaned.RDS")

RFIntmodels <- readRDS("RFPreds_NumericParams.RDS")
RFPipemodels <- readRDS("pipelinesOnlyRFPreds_Params.RDS")

nclusts <- designMat %>%
  select(pipelines, name, nclusts)


silLR <- models$sil
dbLR <- models$db
chLR <- models$ch
gseaLR <- models$gsea

silRF <- RFIntmodels$sil
dbRF <- RFIntmodels$db
chRF <- RFIntmodels$ch
gseaRF <- RFIntmodels$gsea

silPipeLR <- pipeLRModels$sil
dbPipeLR <- pipeLRModels$db
chPipeLR <- pipeLRModels$ch
gseaPipeLR <- pipeLRModels$gsea

silPipeRF <- RFPipemodels$sil
dbPipeRF <- RFPipemodels$db
chPipeRF <- RFPipemodels$ch
gseaPipeRF <- RFPipemodels$gsea

datasetCors <- function(x){
  datasetCors <- x %>%
    group_by(name) %>%
    group_modify(~ as.data.frame(cor(as.numeric(.$preds), as.numeric(.$Actual))))
  colnames(datasetCors)[[2]] <- "Cor"
  return(datasetCors)
}

createPlottingDfs <- function(preds, dataMat, RF=TRUE){
  if(RF){
    df <- cbind(dataMat, preds=preds) %>%
      merge(nclusts, by=c("pipelines", "name")) %>%
      dplyr::rename(Actual = contains("Imp"))
  }else{
    df <- cbind(dataMat, preds=preds[,1]) %>%
      merge(nclusts, by=c("pipelines", "name")) %>%
      dplyr::rename(Actual = contains("Imp"))
  }
  return(df)
}

makeColdataDF <- function(R2df){
  dMat <- designMat %>% 
    as_tibble() %>% 
    mutate_if(numbers_only, as.numeric) %>%
    mutate_if(is.numeric, scale) %>%
    dplyr::slice(which(.$name %in% R2df$name)) %>%
    dplyr::select(c("name", "detected", "sum", "ncells", "ngenes", contains("_"))) %>%
    distinct() %>%
    pivot_longer(!name, names_to="dataType") %>%
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
    group_modify(~as.data.frame(cor(as.numeric(.$preds), as.numeric(.$ARI))))
  ariCors <- as.data.frame(ariCors)
  colnames(ariCors)[[2]] <- "Cor"
  return(ariCors)
}



ari_tbl <- ari %>%
  as_tibble() %>%
  mutate_if(is.numeric, scale) %>%
  pivot_longer(starts_with("E"))
colnames(ari_tbl)[[3]] <- "value"

makeARIHm <- function(ch, db, sil, gsea, RF=TRUE){
  silTestPlot <- createPlottingDfs(sil$testPreds, sil$testData, RF)
  dbTestPlot <- createPlottingDfs(db$testPreds, db$testData, RF)
  chTestPlot <- createPlottingDfs(ch$testPreds, ch$testData, RF)
  gseaTestPlot <- createPlottingDfs(gsea$testPreds, gsea$testData, RF)
  
  silTestDatasetCors <- datasetCors(silTestPlot)
  dbTestDatasetCors <- datasetCors(dbTestPlot)
  chTestDatasetCors <- datasetCors(chTestPlot)
  gseaTestDatasetCors <- datasetCors(gseaTestPlot)
  
  testCors <- dplyr::bind_rows(list(sil=silTestDatasetCors, db=dbTestDatasetCors, ch=chTestDatasetCors, gsea=gseaTestDatasetCors), .id="metric")%>%
    mutate(metric = toupper(metric))
  
  
  testCors$metricf <- factor(testCors$metric, levels=c("CH", "DB", "SIL", "GSEA"))
  
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
  
  ariHmDfs <- list(ariHm = corHmARI, corARIPvals = corHmARIpvals)
  
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
  corHmpvals$dataType <- NULL
  
  corHmDfs <- list(corHm = corHmDf, corPvals = corHmpvals)
  

  return(list(ari=ariHmDfs, cor=corHmDfs))
}
LRIntHms <- makeARIHm(chLR, dbLR, silLR, gseaLR, RF=FALSE)
RFIntHms <- makeARIHm(chRF, dbRF, silRF, gseaRF, RF=TRUE)

LRPipeHms <- makeARIHm(chPipeLR, dbPipeLR, silPipeLR, gseaPipeLR, RF=FALSE)
RFPipeHms <- makeARIHm(chPipeRF, dbPipeRF, silPipeRF, gseaPipeRF, RF=TRUE)

findSigCors <- function(LRPvals, RFPvals) {
  LRSig <- LRPvals %>%
     filter(SIL < 0.05 | DB < 0.05 | CH < 0.05 | GSEA < 0.05)
  
  RFSig <- RFPvals %>%
    filter(SIL < 0.05 | DB < 0.05 | CH < 0.05 | GSEA < 0.05)
  
  allSig <- c(rownames(LRSig), rownames(RFSig)) %>%
    as.data.frame() %>%
    distinct()
  return(allSig)
}

ariSigFeatures <- findSigCors(LRIntHms$ari$corARIPvals, RFIntHms$ari$corARIPvals)
corSigFeatures <- findSigCors(LRIntHms$cor$corPvals, RFIntHms$cor$corPvals)

pipeAriSigFeatures <- findSigCors(LRPipeHms$ari$corARIPvals, RFPipeHms$ari$corARIPvals)
pipeCorSigFeatures <- findSigCors(LRPipeHms$cor$corPvals, RFPipeHms$cor$corPvals)


col_fun = circlize::colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))

makeHeatmaps <- function(sigFeatures, corDf, modelType, corType) {
  corHmDfSig <- corDf[sigFeatures,]
  metrics <- names(corHmDfSig)
  corHmDfSig <- sapply(as.data.frame(corHmDfSig), as.numeric)
  col_ha <- rowAnnotation("Average" = anno_barplot(rowMeans(corHmDfSig)))
  if(length(sigFeatures)==1){
    print(corHmDfSig)
    corHmDfSig <- t(corHmDfSig)
    rownames(corHmDfSig) <- sigFeatures
    colnames(corHmDfSig) <- metrics
    col_ha <- rowAnnotation("Average" = anno_barplot(rowMeans(corHmDfSig)))
    corHm <- ComplexHeatmap::Heatmap(as.matrix(corHmDfSig),
                                     name=sprintf("Correlation of dataset features with dataset-specific predictive performance", corType),
                                     left_annotation = col_ha,
                                     cluster_columns = FALSE,
                                     heatmap_legend_param = list(
                                       legend_direction = "horizontal"),
                                     column_title = modelType,
                                     row_title=corType,
                                     col=col_fun)
  }else{
    rownames(corHmDfSig) <- sigFeatures
    corHm <- ComplexHeatmap::Heatmap(as.matrix(corHmDfSig),
                                     name=sprintf("Correlation of dataset features with dataset-specific predictive performance", corType),
                                     left_annotation = col_ha,
                                     cluster_columns = FALSE,
                                     heatmap_legend_param = list(
                                       legend_direction = "horizontal"),
                                     column_title = modelType,
                                     row_title=corType,
                                     col=col_fun)
  }
  
  return(corHm)
}
LRIntARIHm <- makeHeatmaps(ariSigFeatures[,1], LRIntHms$ari$ariHm, modelType = "LR with interactions", corType = "ARI")
RFIntARIHm <- makeHeatmaps(ariSigFeatures[,1], RFIntHms$ari$ariHm, modelType="RF with interactions", corType = "ARI")

LRCorHm <- makeHeatmaps(corSigFeatures[,1], LRIntHms$cor$corHm, modelType="LR with interactions", corType = "Metric")
RFCorHm <- makeHeatmaps(corSigFeatures[,1], RFIntHms$cor$corHm, modelType="RF with interactions", corType = "Metric")

LRPipeARIHm <- makeHeatmaps(pipeAriSigFeatures[,1], LRPipeHms$ari$ariHm, modelType = "Pipelines only LR", corType="ARI")
RFPipeARIHm <- makeHeatmaps(pipeAriSigFeatures[,1], RFPipeHms$ari$ariHm, modelType = "Pipelines only RF", corType="ARI")

LRPipeCorHm <- makeHeatmaps(pipeCorSigFeatures[,1], LRPipeHms$cor$corHm, modelType = "Pipelines only LR", corType="Metric")
RFPipeCorHm <- makeHeatmaps(pipeCorSigFeatures[,1], RFPipeHms$cor$corHm, modelType = "Pipelines only RF", corType="Metric")


pdf("PaperFigures/pipelinesOnlyAllCorrelationHeatmaps.pdf", height=9.5, width=7.5)
draw(LRPipeARIHm %v% LRPipeCorHm,
     row_title="Feature", 
     column_title="Metric", 
     column_title_side="bottom",
     row_title_side="right",
      heatmap_legend_side="top",
     ht_gap = unit(1.3, "cm")) 

draw(RFPipeARIHm %v% RFPipeCorHm,
     row_title="Feature", 
     column_title="Metric", 
     column_title_side="bottom",
     heatmap_legend_side="top",
     ht_gap = unit(1.3, "cm")) 
dev.off()


pdf("PaperFigures/interactionsAllCorrelationHeatmaps.pdf", height=9.5, width=7.5)
draw(LRIntARIHm %v% LRCorHm,
     row_title="Feature", 
     column_title="Metric", 
     column_title_side="bottom",
     row_title_side="right",
     heatmap_legend_side="top",
     ht_gap = unit(1.5, "cm")) 

draw(RFIntARIHm %v% RFCorHm,
     row_title="Feature", 
     column_title="Metric", 
     column_title_side="bottom",
     heatmap_legend_side="top",
     ht_gap = unit(1.5, "cm")) 
dev.off()


pdf("PaperFigures/LRpipelinesOnlyModelsMergedFeatureCorrelationHeatmapsARI.pdf",
    height=4.5, width=7.5)
draw(LRPipeARIHm %v% LRPipeCorHm ,
     row_title="Feature", 
     column_title="Metric", 
     column_title_side="bottom",
     heatmap_legend_side="top")
dev.off()
pdf("PaperFigures/RFpipelinesOnlyModelsMergedFeatureCorrelationHeatmaps.pdf",
    height=3.5, width=7.5)
draw(RFPipeARIHm %v% RFPipeCorHm,
     row_title="Feature", 
     column_title="Metric", 
     column_title_side="bottom",
     heatmap_legend_side="top")
dev.off()

pdf("PaperFigures/interactionsModelsMergedFeatureCorrelationHeatmapsARI.pdf",
    height=2, width=7.5)
draw(LRIntARIHm %v% LRCorHm ,
     row_title="Feature", 
     column_title="Metric", 
     column_title_side="bottom",
     heatmap_legend_side="top")
dev.off()

pdf("PaperFigures/interactionsModelsMergedFeatureCorrelationHeatmaps.pdf",
    height=3.5, width=7.5)
draw(RFIntARIHm %v% RFCorHm,
     row_title="Feature", 
     column_title="Metric", 
     column_title_side="bottom",
     heatmap_legend_side="top")
dev.off()