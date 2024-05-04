library(Seurat)
library(SeuratData)
library(magrittr)

adata <- lapply(adata, function(x) as.Seurat(x, counts='X', data=NULL))
int_feats <- SelectIntegrationFeatures(adata)
int_anchors <- FindIntegrationAnchors(object.list = adata, anchor.features = int_feats)
cca <- IntegrateData(anchorset = int_anchors)
remove(int_anchors, int_feats)
cca <- cca %>%
  RunPCA(verbose = FALSE) %>%
  FindNeighbors(dims = 1:30) %>%
  RunUMAP(dims = 1:30) %>%
  FindClusters()
resolution.range <- seq(from = 0, to = 1, by = 0.1)
cca <- FindClusters(cca, resolution = resolution.range)
cca <- as.SingleCellExperiment(cca)