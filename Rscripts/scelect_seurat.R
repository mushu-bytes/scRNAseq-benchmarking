library(SeuratObject)
library(Seurat)
library(remotes)
library(sctransform)
library(glmGamPoi)
library(ggplot2)
library(glue)
library(cluster)
library(dplyr)
library(SeuratDisk)

# args are input data, output path, dataset key to split object by
args <- commandArgs(trailingOnly = TRUE)

input_path <- args[[1]]
output_path <- args[[2]]
dataset_key <- args[[3]]

read_single_cell_data <- function(filepath) {
  tryCatch({
    # Try to read as 10X Genomics HDF5
    return(Read10X_h5(filepath))
  }, error = function(e) {
    # Handle error silently and try the next format
  })
  
  tryCatch({
    # Try to read as h5seurat
    return(LoadH5Seurat(filepath))
  }, error = function(e) {
  })
  
  tryCatch({
    # Try to read as rds
    return(readRDS(filepath))
  }, error = function(e) {
  })
  
  tryCatch({
    # Try to read as H5AD
    Convert(filepath, dest = glue("tmp/data.h5seurat"), overwrite = TRUE)
    return(LoadH5Seurat("tmp/data.h5seurat"))
  }, error = function(e) {
  })
  
  stop("Unsupported file format or file could not be read")
}
  
data <- read_single_cell_data(input_path)
data <- SplitObject(data, split.by = dataset_key)[1:5]
normalized_data <- lapply(X = data, FUN = SCTransform, method = "glmGamPoi")

reduction <- list("rpca", "cca")
integrations <- list()

# do the integrations
for (x in reduction) {
    features <- SelectIntegrationFeatures(object.list = normalized_data,
                                          nfeatures = 3000)
    pre_integration <- PrepSCTIntegration(object.list = normalized_data,
                                          anchor.features = features)

    # npcs must be less than the number of cells in the assay
    pre_integration <- lapply(X = pre_integration,
                              FUN = RunPCA,
                              features = features)

    # Very difficult to integrate across datasets with varying cells
    print(x)
    anchors <- FindIntegrationAnchors(object.list = pre_integration,
                                      normalization.method = "SCT",
                                      anchor.features = features,
                                      dims = 1:30,
                                      reduction = x,
                                      k.anchor = 20)

    integrated <- IntegrateData(anchorset = anchors,
                                normalization.method = "SCT",
                                dims = 1:30)
    integrated <- RunPCA(integrated)
    resolution_range <- seq(from = 0, to = 1, by = 0.1)
    integrated <- FindNeighbors(integrated, dims = 1:30)
    integrated  <- RunUMAP(integrated, reduction = "pca", dims = 1:30)
    integrated <- FindClusters(integrated, resolution = resolution_range)
    integrations <- append(integrations, integrated)

    # Adjust the contrast in the plot
    dim_plot <- DimPlot(object = integrated, reduction = "umap", group.by = dataset_key)

    # make sure to add string interpolcation here
    ggsave(filename = glue("{output_path}/{x}_umap_plot.png"), plot = dim_plot)
}

compute_silhouette_avg <- function(integration) {
  clusters <- Idents(integration)
  coords <- Embeddings(integration, reduction = "pca")
  silhouette_scores <- silhouette(as.numeric(clusters), dist(coords))
  sil_df <- data.frame(
    cluster = factor(clusters),
    silhouette_width = silhouette_scores[, 3]
  )
}

# compute silhouette scores and visualize
results1 <- compute_silhouette_avg(integrations[[1]])
results2 <- compute_silhouette_avg(integrations[[2]])
avg_sil_score1 <- mean(results1$silhouette_width)
avg_sil_score2 <- mean(results2$silhouette_width)

results1$dataset <- "Dataset 1"
results2$dataset <- "Dataset 2"

combined_sil_scores <- bind_rows(results1, results2)

p <- ggplot(combined_sil_scores, aes(x = cluster, y = silhouette_width, fill = dataset)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_hline(aes(yintercept = avg_sil_score1, color = "Dataset 1 Avg"), linetype = "dashed", show.legend = TRUE) +
  geom_hline(aes(yintercept = avg_sil_score2, color = "Dataset 2 Avg"), linetype = "dashed", show.legend = TRUE) +
  theme_minimal() +
  labs(title = "Silhouette Scores by Cluster",
       x = "Cluster",
       y = "Silhouette Width") +
  theme(legend.position = "bottom")

ggsave(glue("{output_path}/seurat_report.png"), plot = p, width = 10, height = 6)

# write to disk
if (avg_sil_score1 > avg_sil_score2) {
  SaveH5Seurat(integrations[[1]], glue("{output_path}/seurat_integration.h5Seurat"), overwrite=TRUE)
}
if (avg_sil_score1 < avg_sil_score2) {
  SaveH5Seurat(integrations[[2]], glue("{output_path}/seurat_integration.h5Seurat"), overwrite=TRUE)
}







