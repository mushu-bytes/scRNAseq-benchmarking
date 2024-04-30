library(Seurat)
library(magrittr)
library(SeuratData)
library(sctransform)
library(KernSmooth)
library(parallel)
library(DoubletFinder)
library(glmGamPoi)
library(fields)
library(clustree)
library(ROCR)
library(dplyr)

# Read in samples 
nsg_bus_1.data <- Read10X("/mnt/shared/nationwide/Counts/NSG_BUS_1/outs/filtered_feature_bc_matrix") 
nsg_bus_2.data <- Read10X("/mnt/shared/nationwide/Counts/NSG_BUS_2/outs/filtered_feature_bc_matrix") 
nsg_bus_3.data <- Read10X("/mnt/shared/nationwide/Counts/NSG_BUS_3/outs/filtered_feature_bc_matrix") 

# Create Seurat objects 
nsg_bus_1 <- CreateSeuratObject(counts = nsg_bus_1.data,
                                project = "nsg_bus_1",
                                min.cells = 3)
nsg_bus_2 <- CreateSeuratObject(counts = nsg_bus_2.data,
                                project = "nsg_bus_2",
                                min.cells = 3)
nsg_bus_3 <- CreateSeuratObject(counts = nsg_bus_3.data,
                                project = "nsg_bus_3",
                                min.cells = 3)

# make list of seurat objects and save
objs <- list(nsg_bus_1 = nsg_bus_1,
             nsg_bus_2 = nsg_bus_2,
             nsg_bus_3 = nsg_bus_3)

# calculating QC metrics
for (sample in names(objs)) {
  data <- objs[[sample]]
  data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^mt-")
  objs[[sample]] <- data
}

# subset data based on QC metrics
for (sample in names(objs)){
  data <- objs[[sample]]
  data <- subset(data,
                 subset = nFeature_RNA > 1000 &
                   nFeature_RNA < 6000 &
                   nCount_RNA < 25000 &
                   percent.mt < 5)
  objs[[sample]] <- data
}

# run pre-processing on each object individually
for (sample in names(objs)){
  data <- objs[[sample]]
  data <- data %>%
    SCTransform(method = "glmGamPoi") %>%
    FindVariableFeatures() %>%
    ScaleData() %>%
    RunPCA(npcs = 30) %>%
    FindNeighbors(k.parm = 30) %>%
    FindClusters() %>%
    RunUMAP(dims = 1:30)
  objs[[sample]] <- data
}

# run doublet finder to get subset with just singlets 
for (sample in names(objs)){
  data <- objs[[sample]]
  # pK Identification (no ground-truth)
  sweep.res <- paramSweep(data, PCs = 1:30, sct = T)
  sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  # store pK as a vector
  pk <- as.numeric(as.vector(top_n(bcmvn,1,BCmetric)[["pK"]]))
  # homotypic doublet proportion estimate
  nExp_poi <- round(0.075*nrow(data@meta.data))
  # run DoubletFinder
  data <- doubletFinder(data, PCs = 1:30, pN = 0.25, pK = pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = T)
  # label doublets
  Idents(data) <- data@meta.data[[paste0("DF.classifications_0.25_",pk,"_",nExp_poi)]]
  # plot number of doublets and singlets
  df <- matrix(data = NA, nrow = 2, ncol = 1) %>%
    as.data.frame()
  colnames(df) <- c("count")
  rownames(df) <- c("Singlet", "Doublet")
  df$DoubletFinder <- rownames(df)
  df["Singlet", "count"] <- sum(data@meta.data$DF.classifications_0.25_0.05_727 == "Singlet")
  df["Doublet", "count"] <- sum(data@meta.data$DF.classifications_0.25_0.05_727 == "Doublet")
  # remove doublets
  data <- subset(data, idents = "Singlet")
  # store back in list
  objs[[sample]] <- data
}

#Add metadata 
objs[["nsg_bus_1"]]$condition <- "busulfan"
objs[["nsg_bus_2"]]$condition <- "busulfan"
objs[["nsg_bus_3"]]$condition <- "busulfan"
objs[["nsg_control_1"]]$condition <- "control"
objs[["nsg_control_2"]]$condition <- "control"
objs[["nsg_control_3"]]$condition <- "control"

objs[["nsg_bus_1"]]$pool <- "3F"
objs[["nsg_bus_2"]]$pool <- "2F2M"
objs[["nsg_bus_3"]]$pool <- "2M"
objs[["nsg_control_1"]]$pool <- "3F"
objs[["nsg_control_2"]]$pool <- "2F2M"
objs[["nsg_control_3"]]$pool <- "2M"

# perform cca integration
seurat_objs <- objs[1:3]

int_feats <- SelectIntegrationFeatures(seurat_objs)
int_list <- PrepSCTIntegration(object.list = seurat_objs,
                               anchor.features = int_feats)
int_anchors <- FindIntegrationAnchors(object.list = int_list,
                                      normalization.method = "SCT",
                                      anchor.features = int_feats)
cca <- IntegrateData(anchorset = int_anchors,
                     normalization.method = "SCT")
remove(int_anchors, int_list, int_feats)
cca <- cca %>%
  RunPCA(verbose = FALSE) %>%
  FindNeighbors(dims = 1:30) %>%
  RunUMAP(dims = 1:30) %>%
  FindClusters()
resolution.range <- seq(from = 0, to = 1, by = 0.1)
cca <- FindClusters(cca, resolution = resolution.range)
tree <- clustree(cca)
cca <- FindClusters(cca, resolution = 0.3)

# perform rpca integration
int_feats <- SelectIntegrationFeatures(seurat_objs)
int_list <- PrepSCTIntegration(object.list = seurat_objs,
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
resolution.range <- seq(from = 0, to = 1, by = 0.1)
rpca <- FindClusters(rpca, resolution = resolution.range)
tree <- clustree(rpca)
rpca <- FindClusters(rpca, resolution = 0.3)

