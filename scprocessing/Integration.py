import scanpy as sc
import anndata as ad
from typing import List
from anndata import AnnData
import scanpy.external as sce
from .PipelineStep import PipelineStep


class Integration(PipelineStep):
    def __init__(
        self, method: str = "scanorama", key: str = "Trial", resolution: int = 0.3
    ):
        """
        Method: Integration Method
        Key: Key differentiating different datasets
        Resolution: Clustering resolution
        will add additional parameters for normalization
        """
        self.method = method
        self.key = key
        self.resolution = resolution

    def apply(self, datasets: List[AnnData]) -> AnnData:
        """
        Parameters:
            Datasets: list of datasets to concatenate

        Return Value:
            AnnData object containing PCA, Integration, Nearest Neighbor, UMAP, and Clustering results
        """
        print("Integrating Datasets")
        dataset = ad.concat(datasets)
        sc.pp.pca(dataset)
        if self.method == "scanorama":
            sce.pp.scanorama_integrate(dataset, self.key, verbose=1)
            sc.pp.neighbors(dataset, use_rep="X_scanorama")
        elif self.method == "harmony":
            sce.pp.harmony_integrate(dataset, self.key, verbose=1)
            sc.pp.neighbors(dataset, use_rep="X_pca_harmony")
        elif self.method == "merge":
            sc.pp.neighbors(dataset, use_rep="X_pca")
        else:
            raise ValueError("Invalid Integration Method")
        sc.tl.umap(dataset)
        sc.tl.leiden(
            dataset,
            key_added="clusters",
            n_iterations=2,
            directed=False,
            resolution=self.resolution,
        )
        sc.pl.umap(dataset, color=["clusters"], palette=sc.pl.palettes.default_20)
        return dataset
