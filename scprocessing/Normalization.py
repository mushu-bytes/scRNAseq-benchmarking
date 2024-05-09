import scanpy as sc
import anndata as ad
from typing import List
from anndata import AnnData
from .PipelineStep import PipelineStep


class Normalization(PipelineStep):
    def __init__(self, method: str = "zheng17"):
        """
        method: A string indicating how the data should be normalized
        """
        self.method = method

    def apply(self, datasets: List[AnnData]) -> List[AnnData]:
        """
        Parameters:
            dataset: AnnData object
        Return Value: A normalized AnnData object
        """
        print("Normalizing Datasets")
        Normalized_data = []
        for data in datasets:
            if self.method == "zheng17":
                Normalized_data.append(sc.pp.recipe_zheng17(data, copy=True))
            elif self.method == "seurat":
                Normalized_data.append(sc.pp.recipe_seurat(data, copy=True))
            elif self.method == "weinreb17":
                Normalized_data.append(sc.pp.recipe_weinreb17(data, copy=True))

        return Normalized_data
