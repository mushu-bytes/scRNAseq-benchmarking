import scanpy as sc
import anndata as ad
from typing import List, Union, Callable
from anndata import AnnData
import scanpy.external as sce
from .PipelineStep import PipelineStep


class Integration(PipelineStep):
    def __init__(
        self,
        method: Union[str, Callable] = "scanorama",
        key: str = "Type",
        resolution: int = 0.3,
    ):
        """
        Method: Integration Method
        Key: Key differentiating different datasets
        Resolution: Clustering resolution
        will add additional parameters for normalization
        """
        methods = {
            "scanorama": sce.pp.scanorama_integrate,
            "harmony": sce.pp.harmony_integrate,
            "merge": None,  # just do nothing if merge
        }
        if type(method) == str:
            if method in methods:
                self.integration = methods[method]
            else:
                raise ValueError("Invalid Integration Method")
        else:
            self.integration = method

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

        # integrate, or if merge don't
        if self.integration:
            self.integration(dataset, self.key, adjusted_basis="X_integration")
            sc.pp.neighbors(dataset, use_rep="X_integration")
        else:
            sc.pp.neighbors(dataset)

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
