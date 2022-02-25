suppressPackageStartupMessages({
  library(SingleCellExperiment) # BioConductor
  library(DropletUtils) # BioConductor
  library(tidyverse) # CRAN
  library(here) # CRAN
  library(DT) # CRAN
  library(here) # CRAN
  library(Matrix) # CRAN
  library(DESeq2) # BioConductor
  library(pipeComp) # Github
  library(scran) # BioConductor
  library(fgsea) # BioConductor
  library(org.Hs.eg.db) # BioConductor
  library(AnnotationDbi) # BioConductor
  library(purrr)
  library(scuttle)
})

# Reading in the SCE and clusters
args = commandArgs(trailingOnly = TRUE)
print(args)
sce_name <- args[[1]]
sce <- readRDS(sce_name)
clusters_name <- args[[2]]
clusters <- readRDS(clusters_name)
whichClust <- args[[3]]
pathways.hallmark <- gmtPathways("/home/campbell/cfang/automl_scrna/scripts/Rscripts/c5.all.v7.4.symbols.gmt")
l <- lengths(pathways.hallmark)
pathways.hallmark <- pathways.hallmark[which(l > 10 & l < 500)]
# pathways.hallmark <- gmtPathways("c5.all.v7.4.symbols.gmt")
# sce <- readRDS("E-CURD-10-SCE.RDS")
# clusters <- readRDS("E-CURD-10res.dataset1.endOutputs.rds")
# whichClust <- 150

exp.name <- strsplit(sce_name, "/")[[1]][[8]]
print(exp.name)
whichClust <- as.numeric(whichClust)
print(clusters[[whichClust]])
print(length(unique(clusters[[whichClust]])))
if (length(unique(clusters[[whichClust]])) <= 1){
  run.name <- whichClust
  print("single clust")
  filename <- paste(run.name, "-gsea.csv", sep="")
  filename <- paste(exp.name, filename, sep="-")
  gsea_filename <- paste("/home/campbell/cfang/automl_scrna/results/gsea_results", filename, sep="/")
  print(gsea_filename)
  write.csv("Single cluster", gsea_filename)
} else {
  cl <- clusters[[whichClust]]
  filt_sce <- sce[,rownames(as.matrix(cl))]
  filt_sce <- logNormCounts(filt_sce)
  fms <- findMarkers(filt_sce, groups = cl)
  logFCs <- lapply(fms,
                   function(fm){
                     print(fm)
                     logfc <- fm[[5]]
                     names(logfc) <- rownames(fm)
                     logfc
                   })
  mapped_logFCs <- lapply(logFCs, function(logfc){
    ens2symbol <- AnnotationDbi::select(org.Hs.eg.db,
                                        key=names(logfc),
                                        columns="SYMBOL",
                                        keytype="ENSEMBL")
    ens2symbol <- as_tibble(ens2symbol)
    res <- inner_join(as_tibble(logfc, rownames="row"), ens2symbol, by=c("row"="ENSEMBL"))
    res2 <- res %>%
      dplyr::select(SYMBOL, value) %>%
      na.omit() %>%
      distinct() %>%
      group_by(SYMBOL)
    ranks <- deframe(res2)
    ranks
  })
  
  gseas <- lapply(mapped_logFCs, function(mapped_logfc){
    fgseaRes <- fgsea(mapped_logfc, pathways = pathways.hallmark)
    fgseaResTidy <- fgseaRes %>%
      as_tibble() %>%
      arrange(desc(NES))
    mean(abs(fgseaResTidy$ES))
  })
  
  print(gseas)
  gseas_df <- as.data.frame(gseas)
  run.name <- whichClust
  filename <- paste(run.name, "-gsea.csv", sep="")
  filename <- paste(exp.name, filename, sep="-")
  gsea_filename <- paste("/home/campbell/cfang/automl_scrna/results/gsea_results", filename, sep="/")
  write.csv(gseas_df, gsea_filename)
}
