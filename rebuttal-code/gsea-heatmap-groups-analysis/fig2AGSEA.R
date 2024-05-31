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
designMatPath <- here("data", "uncorrected_designMatrix_scaled.RDS")
numClustersPath <- here("data", "num_clusters_cleaned.RDS")
designMatCor <- readRDS(here("data", "corrected_designMat_unscaled.RDS"))

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
#pdf("PaperFigures/MetricsHeatmaps.pdf")
#pdf("/home/campbell/cfang/automl_scrna/results/figures/uncorrectedMetricsHeatmaps.pdf")
gseaHM_all <- draw(chHM + dbHM + silHM + gseaHM,
     row_title="Pipeline",
     column_title="Dataset",
     column_title_side="bottom",
     heatmap_legend_side="right",
     annotation_legend_side="right",
     ht_gap=unit(0.055, "cm"))
#dev.off()

#gseaHM_all = draw(chHM + dbHM + silHM + gseaHM)
dend <- ComplexHeatmap::column_dend(gseaHM_all)
data_order <- column_order(gseaHM_all)

col_order <- colnames(gseaMat)[data_order[[4]]]
dend[[1]]
dend[[2]]

firstHalf <- col_order[1:35]
secondHalf <- col_order[36:86]

des <- designMat %>%
  mutate(split=case_when(name %in% firstHalf ~ "Positive",
                         name %in% secondHalf ~ "Negative")) %>%
  select(!c(sil, ch, db, nclusts, pipelines, filt, dims, norm, res, sum.1, total.1, detected.1,
            pct_counts_in_top_50_features, pct_counts_top_20_features,pct_counts_Mt, pct_coding,
            total_features, total_counts, sum,pct_ribosomal, name)) %>%
  distinct()



des_long <- des %>%
  pivot_longer(cols=!c("split"))


ggsave(here("results", "figures", "supplementalFigs", "gsea-analyses", "allGSEAboxplots.pdf"), height=15, width=25, dpi=300)
ggplot(des_long, aes(x=split, y=as.numeric(value)))+
  geom_boxplot()+
  geom_jitter()+
  facet_wrap(~name, scales="free")+
  cowplot::theme_cowplot()+
  theme(strip.background=element_blank(),
        axis.title=element_text(size=40),
        title=element_text(size=30),
        strip.text.x=element_text(size=19, hjust=0),
        axis.text.x=element_text(size=22), axis.text.y=element_text(size=22))+
  xlab("Split")+
  ylab("Value")+
  ggtitle("Dataset summary statistics of each group")
dev.off()

pdf(here("results", "figures", "supplementalFigs", "gsea-analyses", "gsea-ncells-boxplots.pdf"))
ggplot(des, aes(x=split, y=as.numeric(ncells)))+
      geom_boxplot()+
  geom_jitter()+
  xlab("Split")+
  ylab("Number of cells")+
  cowplot::theme_cowplot()+
  ggtitle("Number of cells per\ndataset in each group")+
  theme(axis.title.x=element_text(size=30),
        axis.title.y=element_text(size=30),
        title=element_text(size=32))
dev.off()

# plist <- list()
# for (i in 1:length(cols_plt)){
#   col <- cols_plt[[i]]
#   print(head(des[,col]))
#   p <- ggplot(des, aes(x=split, y=as.numeric(des[,col])), environment=environment())+
#     geom_boxplot()+
#     ylab(col)
#   p
#   plist[[i]] <- p
# }
# do.call(gridExtra::grid.arrange, c(plist, ncol=4))

# firstHalfDesign <- designMat %>%
#   filter(name %in% firstHalf) %>%
#   select(!c(db, sil, ch, nclusts, pipelines, filt, dims, norm, res))%>%
#   data.matrix()
#
#
# secondHalfDesign <- designMat %>%
#   filter(name %in% secondHalf)%>%
#   select(!c(db, sil, ch, nclusts, pipelines, filt, dims, norm, res))%>%
#   data.matrix()
#
# design<- designMat %>%
#   select(!c(db, sil, ch, nclusts, pipelines, filt, dims, norm, res))%>%
#   dplyr::distinct()
#
# rownames(design) <- design$name
#
#
# HM <- ComplexHeatmap::Heatmap(data.matrix(design), cluster_rows=TRUE)
# HM = draw(HM)
#
# HM_order <- row_order(HM)
# HM_dend <- ComplexHeatmap::row_dend(HM)
#
# HM_rownames <- rownames(design)[HM_order]
#
# HM_dend[[1]]
# HM_dend[[2]]
# firstGroup <- HM_rownames[1:36]
# secondGroup <- HM_rownames[37:86]
