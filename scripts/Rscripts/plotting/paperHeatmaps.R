suppressPackageStartupMessages({
  library(tidyverse) # CRAN
  library(here) # CRAN
  library(DT) # CRAN
  library(pheatmap) # CRAN
  library(here) # CRAN
  library(Matrix) # CRAN
  library(ggplot2)
  library(purrr)
  library(ComplexHeatmap)
  library(circlize)
  library(viridis)
})
designMatPath <- "uncorrected_designMatrix_scaled.RDS"
numClustersPath <- "num_clusters_cleaned.RDS"
designMatCor <- readRDS("corrected_designMat_unscaled.RDS")

# args <- commandArgs(trailingOnly = TRUE)
# designMatPath <- args[[1]]
# numClustersPath <- args[[2]]
# designMatCor <- readRDS(args[[3]])

nclusts <- readRDS(numClustersPath)

designMat <- as_tibble(readRDS(designMatPath))
designMat$pipelines <- gsub(pattern="\"", replacement="", designMat$pipelines)
designMat <- designMat %>%
  dplyr::rename(filt = V2)%>%
  dplyr::rename(dims = V3)%>%
  dplyr::rename(norm = V4)%>%
  dplyr::rename(res = V5)

# Merge nclusts with designMat
nclusts <- as_tibble(nclusts)
#nclusts <- pivot_longer(nclusts, cols=which(grepl("^E-*", colnames(nclusts))))
designClusts <- merge(designMat, nclusts, by=c("pipelines", "name"))
colnames(designClusts)[which(colnames(designClusts)=="value")] <- "nclusts"
designMat <- designClusts


metrics <- designMat %>%
  select(c(name, pipelines, sil, db, ch, filt, dims, norm, res))

silMat <- metrics %>%
  select(c(name, pipelines, sil,  filt, dims, norm, res)) %>%
  pivot_wider(names_from=name, values_from=sil)%>%
  column_to_rownames("pipelines")

dbMat <- metrics %>%
  select(c(name, pipelines, db, filt, dims, norm, res)) %>%
  pivot_wider(names_from=name, values_from=db)%>%
  column_to_rownames("pipelines")

chMat <- metrics %>%
  select(c(name, pipelines, ch, filt, dims, norm, res)) %>%
  pivot_wider(names_from=name, values_from=ch)%>%
  column_to_rownames("pipelines")

resdf <- silMat %>%
  as_tibble()%>%
  select(res)

normdf <- silMat %>%
  as_tibble()%>%
  select(norm)

dimsdf <- silMat %>%
  as_tibble()%>%
  select(dims)

filtdf <- silMat %>%
  as_tibble()%>%
  select(filt)
  # mutate(filt=case_when(filt=="lenient" ~ 1,
  #                       filt=="default" ~ 2,
  #                       filt=="stringent"~ 3))

dbMat <- dbMat %>%
  as_tibble() %>%
  select(-c(filt, dims, norm, res))

silMat <- silMat %>%
  as_tibble() %>%
  select(-c(filt, dims, norm, res))

chMat <- chMat %>%
  as_tibble() %>%
  select(-c(filt, dims, norm, res))

designMatCor <- designMatCor %>%
  as_tibble() %>%
  merge(nclusts, by=c("name", "pipelines"))

gseaMat <- designMatCor %>%
  select(c(name, pipelines, gseaScaledImp, filt, dims, norm, res)) %>%
  pivot_wider(names_from=name, values_from=gseaScaledImp)%>%
  column_to_rownames("pipelines")%>%
  as_tibble()%>%
  select(-c(filt, dims, norm, res))%>%
  mutate_if(is.character, as.numeric)


res_col = circlize::colorRamp2(as.numeric(unique(resdf$res)), hcl.colors(8, "Teal"))
dims_col = circlize::colorRamp2(as.numeric(unique(dimsdf$dims)), hcl.colors(4, "YlOrRd"))
#filt_col = circlize::colorRamp2(c(3,2,1), hcl.colors(3, "Blues"))

n=nrow(chMat)*ncol(chMat)
colFunCH = circlize::colorRamp2(seq(-2, 2, length = n), hcl.colors(n,"viridis"))
colFunDB = circlize::colorRamp2(seq(-2, 2, length = n), hcl.colors(n,"viridis"))
colFunSil = circlize::colorRamp2(seq(-2, 2, length = n), hcl.colors(n,"viridis"))
colFunGsea = circlize::colorRamp2(seq(-2, 2, length = n), hcl.colors(n,"viridis"))


sil_row_ha = rowAnnotation("Filtering\nstrategy" = filtdf$filt, 
                           "Normalization\nstrategy" = normdf$norm, 
                            "Clustering\nresolution" = as.numeric(resdf$res),
                           "Dimensions"=as.numeric(dimsdf$dims),
                           col = list("Filtering\nstrategy" = c("lenient" = "#F4FAFEFF", "default" = "#7FABD3FF", "stringent" = "#273871FF"),
                                      "Normalization\nstrategy" = c("seurat" = "#089392", "sctransform" = "#EAE29C", "scran" = "#CF597E"), 
                                      "Clustering\nresolution"= res_col,
                                      "Dimensions"=dims_col),
                           show_annotation_name=F,
                           width=unit(1, "cm"),
                           simple_anno_size = unit(0.27, "cm"),
                           annotation_legend_param=list(
                             "Filtering\nstrategy" = list(nrow=3, labels_gp = gpar(fontsize = 7),
                                                          title_gp = gpar(fontsize=8, fontface="bold")),
                             "Normalization\nstrategy" = list(nrow=3,labels_gp = gpar(fontsize = 7),
                                                              title_gp = gpar(fontsize=8, fontface="bold")),
                             "Clustering\nresolution"= list(direction="vertical",labels_gp = gpar(fontsize = 7), 
                                                           legend_height=unit(1.6, "cm"),
                                                           title_gp = gpar(fontsize=8, fontface="bold")),
                             "Dimensions"=list(direction="vertical",labels_gp = gpar(fontsize = 7),
                                               legend_height=unit(1.6, "cm"),
                                               title_gp = gpar(fontsize=8, fontface="bold"))))

chHM <- ComplexHeatmap::Heatmap(as.matrix(chMat), col=colFunCH, name="Metric value", 
                                column_names_gp = grid::gpar(fontsize = 0),
                                row_names_gp = gpar(fontsize=0), 
                                width=unit(3.1, "cm"), 
                                column_title="CH",
                                height = unit(10.1, "cm"),
                                heatmap_legend_param = list(direction="vertical",
                                                            labels_gp = gpar(fontsize = 7),
                                                            legend_height=unit(1.6,"cm"),
                                                            title_gp = gpar(fontsize=8, fontface="bold")))
dbHM <- ComplexHeatmap::Heatmap(as.matrix(dbMat),col=colFunDB, name="Metric value", 
                                column_names_gp = grid::gpar(fontsize = 0),
                                row_names_gp = gpar(fontsize=0), 
                                width=unit(3.1, "cm"),
                                column_title="DB",
                                height = unit(10.1, "cm"))
silHM <- ComplexHeatmap::Heatmap(as.matrix(silMat), col=colFunSil, name="Metric value", 
                                  column_names_gp = grid::gpar(fontsize = 0),
                                  row_names_gp = gpar(fontsize=0), 
                                  width=unit(3.1, "cm"),
                                  column_title="SIL",
                                 height = unit(10.1, "cm"))
gseaHM <- ComplexHeatmap::Heatmap(as.matrix(gseaMat), col=colFunGsea, name="Metric value",
                                 right_annotation=sil_row_ha,
                                 column_names_gp = grid::gpar(fontsize = 0),
                                 row_names_gp = gpar(fontsize=0), 
                                 width=unit(3.1, "cm"),
                                 column_title="GSEA",
                                 height = unit(10.1, "cm"),
                                 heatmap_legend_param = list(direction = "horizontal",
                                                             legend_gp = gpar(fontsize = 6)))
pdf("PaperFigures/MetricsHeatmaps.pdf")
#pdf("/home/campbell/cfang/automl_scrna/results/figures/uncorrectedMetricsHeatmaps.pdf")
draw(chHM + dbHM + silHM + gseaHM, 
     row_title="Pipeline", 
     column_title="Dataset",
     column_title_side="bottom",
     heatmap_legend_side="right",
     annotation_legend_side="right",
     ht_gap=unit(0.055, "cm"))
dev.off()


RFTrain <- readRDS("corrected_RF_trainMatNumeric_unscaled.RDS")
RFTest <- readRDS("corrected_RF_testMatNumeric_unscaled.RDS")

trainMat <- RFTrain$data
testMat <- RFTest$data 

rfMat <- rbind(trainMat, testMat) %>%
  as_tibble()%>%
  select(-c(pipelines, filt, dims, norm ,res))%>%
  distinct()%>%
  rename_with( ~ gsub("V", "PC", .x, fixed = TRUE))
# 
# for (i in 1:length(names(rfMat))){
#   if (grepl("^PC", names(rfMat)[[i]])){
#     num <- as.numeric(strsplit(names(rfMat)[[i]], split="PC")[[1]][[2]])-1
#     names(rfMat)[[i]] <- paste0("PC", num)
#   }
# }

pdf("PaperFigures/datasetFeaturesHeatmap.pdf")
#pdf("/home/campbell/cfang/automl_scrna/results/figures/datasetFeaturesHeatmap.pdf")
colfun <- circlize::colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
ComplexHeatmap::Heatmap(t(as.matrix(rfMat)),
                        col=colfun, 
                        name="Scaled\nvalue",
                        row_title="Feature",
                        column_title="Dataset",
                        column_title_side = "bottom",
                        row_names_gp = grid::gpar(fontsize = 7),
                        heatmap_width = unit(16, "cm"), 
                        heatmap_height = unit(11, "cm"))
dev.off()


designMat <- designMatCor%>%
  select(!contains(c("ch","db","sil", "gsea")))%>%
  select(!nclusts)%>%
  select(!pipelines)%>%
  select(!c("filt", "dims", "norm", "res"))%>%
  distinct()

designMat_long <- designMat %>%
  select(name, subsets_Mt_percent, subsets_ribosomal_percent, 
         PC4, PC6, PC11, PC16, PC18, PC17, PC20,
         subsets_Mt_sum, subsets_ribosomal_sum, subsets_coding_sum) %>%
  tidyr::pivot_longer(!name, names_to="metric", values_to="value")

designMat_sum <- designMat_long %>%
  filter(metric %in% c("subsets_coding_sum", "subsets_Mt_sum", "subsets_ribosomal_sum"))%>%
  filter(value <= 1000)

ggplot(designMat_sum, aes(x=value))+
  geom_histogram()+
  facet_wrap(~metric, scales="free")
