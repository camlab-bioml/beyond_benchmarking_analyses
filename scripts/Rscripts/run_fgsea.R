suppressPackageStartupMessages({
  #library(scater) # BioConductor
  library(SingleCellExperiment) # BioConductor
  library(DropletUtils) # BioConductor
  library(tidyverse) # CRAN
  library(here) # CRAN
  library(DT) # CRAN
  #library(pheatmap) # CRAN
  library(here) # CRAN
  library(Matrix) # CRAN
  library(DESeq2) # BioConductor
  library(pipeComp) # Github
  library(scran) # BioConductor
  library(fgsea) # BioConductor
  library(org.Hs.eg.db) # BioConductor
  library(AnnotationDbi) # BioConductor
  #library(ggplot2)
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
pathways.hallmark <- gmtPathways("/home/campbell/cfang/automl_scrna/scripts/Rscripts/h.all.v7.2.symbols.gmt")
#sce <- readRDS("E-CURD-10-SCE.RDS")
#clusters <- readRDS("E-CURD-10res.dataset1.endOutputs.rds")[[1]]

exp.name <- strsplit(sce_name, "/")[[1]][[8]]
print(exp.name)
removeInds <- c()
for (i in 1:length(clusters)){
  if (length(unique(clusters[[i]])) <= 1){
    run.name <- i
    filename <- paste(run.name, "-gsea.csv", sep="")
    filename <- paste(exp.name, filename, sep="-")
    gsea_filename <- paste("/home/campbell/cfang/automl_scrna/results/gsea_results", filename, sep="/")
    print(gsea_filename)
    write.csv("Single cluster", gsea_filename)
    removeInds <- c(removeInds, i)
  }
}
print(removeInds)
# fms <- lapply(clusters, function(cl) {
#   for (i in 1:length(clusters)){
#     if (i %in% removeInds){
#       print(i)
#       return(1)
#     }
#   }
#   filt_sce <- sce[,rownames(as.matrix(cl))]
#   filt_sce <- logNormCounts(filt_sce)
#   findMarkers(filt_sce, groups = cl)
# }
# )
fms <- c()
for (i in 1:length(clusters)){
  if (i %in% removeInds){
    print(i)
    fms[[i]] <- 1
  }
  else {
  cl <- clusters[[i]]
  filt_sce <- sce[,rownames(as.matrix(cl))]
  filt_sce <- logNormCounts(filt_sce)
  fms[[i]] <- findMarkers(filt_sce, groups = cl)
  }
}
fms.name <- paste("/home/campbell/cfang/automl_scrna/results/gsea_results/", exp.name)
fms.name <- paste(fms.name, "-markers.RDS")
saveRDS(fms, fms.name)

print(fms)
for (i in 1:length(fms)){
  if (i %in% removeInds){
    print(i)
  }
  else{
  print(fms[[i]])
  logFCs <- lapply(fms[[i]],
                   function(fm){
                     print(fm)
                     logfc <- fm[[5]]
                     names(logfc) <- rownames(fm)
                     logfc
                   }
  )
  
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
  ## mean(abs(enrichment_scores))
  gseas_df <- as.data.frame(gseas)
  print(names(fms[[i]]))
  rownames(gseas_df) <- names(fms)[[i]]
  run.name <- i
  filename <- paste(run.name, "-gsea.csv", sep="")
  filename <- paste(exp.name, filename, sep="-")
  gsea_filename <- paste("/home/campbell/cfang/automl_scrna/results/gsea_results", filename, sep="/")
  write.csv(gseas_df, gsea_filename)
  }
  
}



#Save fgsea results

