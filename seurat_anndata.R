library(Seurat)
library(magrittr)
library(SeuratData)
library(sctransform)
library(anndata)

data <- read_h5ad("nsg_bus_0_data.h5ad")
