import scanpy as sc
import anndata as ad
from typing import List, Union, Callable
from anndata import AnnData
from .PipelineStep import PipelineStep


class Normalization(PipelineStep):
    def __init__(self, method: Union[str, Callable] = "zheng17"):
        """
        method: A string indicating how the data should be normalized
        """
        methods = {
            "zheng17": sc.pp.recipe_zheng17,
            "seurat": sc.pp.recipe_seurat,
            "weinreb17": sc.pp.recipe_weinreb17,
        }
        if type(method) == str:
            if method in methods:
                self.normalization = methods[method]
            else:
                raise ValueError("Invalid Normalization Method")
        else:
            self.normalization = method

    def apply(self, datasets: List[AnnData]) -> List[AnnData]:
        """
        Parameters:
            dataset: AnnData object
        Return Value: A normalized AnnData object
        """
        print("Normalizing Datasets")
        Normalized_data = []
        for data in datasets:
            Normalized_data.append(self.normalization(data, copy=True))
        return Normalized_data
