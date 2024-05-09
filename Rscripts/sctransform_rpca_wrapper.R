library(Seurat)
library(SeuratData)
library(magrittr)
library(sctransform)
library(glmGamPoi)

# convert anndata objects to Seurat
adata <- lapply(adata, function(x) as.Seurat(x, counts='X', data=NULL))
# apply sctransform
adata <- lapply(adata, function(x) {
  x <- x %>%
    SCTransform(method = "glmGamPoi", assay="originalexp") %>%
    FindVariableFeatures() %>%
    ScaleData() %>%
    RunPCA(npcs = 30) %>%
    FindNeighbors(k.parm = 30) %>%
    FindClusters() %>%
    RunUMAP(dims = 1:30)
})

int_feats <- SelectIntegrationFeatures(adata)
int_list <- PrepSCTIntegration(object.list = adata,
                               anchor.features = int_feats)
int_anchors <- FindIntegrationAnchors(object.list = int_list,
                                      normalization.method = "SCT",
                                      anchor.features = int_feats,
                                      reduction = "rpca")
rpca <- IntegrateData(anchorset = int_anchors,
                      normalization.method = "SCT")
remove(int_anchors, int_list, int_feats)
rpca <- rpca %>%
  RunPCA(verbose = FALSE) %>%
  FindNeighbors(dims = 1:30) %>%
  RunUMAP(dims = 1:30) %>%
  FindClusters()
# resolution.range <- seq(from = 0, to = 1, by = 0.1)
# rpca <- FindClusters(rpca, resolution = resolution.range)
rpca <- FindClusters(rpca, resolution = 0.3)

# convert back to anndata object
rpca <- as.SingleCellExperiment(rpca)