library(Seurat)
library(SeuratData)
library(SeuratDisk)
library(sctransform)
library(magrittr)
library(glmGamPoi)

adata <- as.Seurat(adata, counts="X", data=NULL)

adata <- adata %>%
          SCTransform(method = "glmGamPoi", assay="originalexp") %>%
          FindVariableFeatures() %>%
          ScaleData() %>%
          RunPCA(npcs = 30) %>%
          FindNeighbors(k.parm = 30) %>%
          FindClusters() %>%
          RunUMAP(dims = 1:30)


adata <- as.SingleCellExperiment(adata)