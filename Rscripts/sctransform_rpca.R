# Enter commands in R (or R studio, if installed)
library(SeuratObject)
library(Seurat)
library(remotes)
library(sctransform)
library(glmGamPoi)
library(SeuratData)
library(ggplot2)
library(clustree)

path = "/mnt/shared/nationwide/cell_type_datasets/immune.rds"
data = readRDS(path)
data.list = SplitObject(data, split.by="development_stage")

target_data = list(data.list$`31-year-old human stage`,
                   data.list$`30-year-old human stage`)

# normalization
target_data <- lapply(X = target_data, FUN = SCTransform, method = "glmGamPoi")

features <- SelectIntegrationFeatures(object.list = target_data, nfeatures = 3000)
target_data <- PrepSCTIntegration(object.list = target_data, anchor.features = features)

# npcs must be less than the number of cells in the assay
target_data <- lapply(X = target_data, FUN = RunPCA, features = features)

# Very difficult to integrate across datasets with varying cells
anchors <- FindIntegrationAnchors(object.list = target_data, normalization.method = "SCT",
                                  anchor.features = features, dims = 1:30, reduction = "rpca", k.anchor = 20)

combined.sct <- IntegrateData(anchorset = anchors, normalization.method = "SCT", dims = 1:30)
combined.sct <- RunPCA(combined.sct)
resolution.range <- seq(from = 0, to = 1, by = 0.1)
combined.sct <- FindNeighbors(combined.sct, dims=1:30)
combined.sct  <- RunUMAP(combined.sct, reduction = "pca", dims = 1:30)
combined.sct <- FindClusters(combined.sct, resolution = resolution.range)

tree <- clustree(combined.sct)
ggsave(tree,
       filename = paste0("/Users/damonlin/Bioinformatics/nationwide/rpca_clustree.jpeg"),
       width = 8,
       height = 10,
       units = "in")
combined.sct <- FindClusters(combined.sct, resolution = 0.3)

# Adjust the contrast in the plot
DimPlot(object = combined.sct, reduction = "umap")
