library(Seurat)
library(SeuratData)
library(magrittr)

adata <- lapply(adata, function(x) as.Seurat(x, counts='X', data=NULL))
int_feats <- SelectIntegrationFeatures(adata)
int_anchors <- FindIntegrationAnchors(object.list = adata,
                                      anchor.features = int_feats,
                                      reduction = "rpca")
rpca <- IntegrateData(anchorset = int_anchors)
remove(int_anchors, int_feats)
rpca <- rpca %>%
  RunPCA(verbose = FALSE) %>%
  FindNeighbors(dims = 1:30) %>%
  RunUMAP(dims = 1:30) %>%
  FindClusters()
resolution.range <- seq(from = 0, to = 1, by = 0.1)
rpca <- FindClusters(rpca, resolution = resolution.range)
rpca <- as.SingleCellExperiment(rpca)