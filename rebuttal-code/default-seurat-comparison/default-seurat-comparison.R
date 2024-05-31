library(here)
library(Seurat)
library(SingleCellExperiment)
library(biomaRt)


sce <- readRDS(here("data", "experiments", "E-MTAB-8410-SCE.RDS"))
sce <- scuttle::logNormCounts(sce)
rowData(sce)$ensembl_gene_id <- rownames(sce)


mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <-  rownames(sce)
gene_IDs <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),
                  values = genes, mart= mart)


rowData(sce) <- merge(rowData(sce), gene_IDs, by=c("ensembl_gene_id"), all=TRUE)

rownames(sce) <- scuttle::uniquifyFeatureNames(ID=rowData(sce)$ensembl_gene_id, names=rowData(sce)$hgnc_symbol)



seu <- as.Seurat(sce)

#### QC ####
seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu)
seu <- ScaleData(seu)
seu <- RunPCA(seu)
seu <- FindNeighbors(seu)
seu <- FindClusters(seu)
seu <- RunUMAP(seu, dims=1:30)
DimPlot(seu,reduction = "umap")

saveRDS(seu, here("data", "experiments", "E-MTAB-8410-SEU.RDS"))


