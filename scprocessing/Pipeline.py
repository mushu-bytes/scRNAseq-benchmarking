import scanpy as sc
import anndata as ad
from typing import List
from anndata import AnnData
import scanpy.external as sce
from .PipelineStep import PipelineStep

class Pipeline():
    def __init__(self):
        self.steps = []
    def add_step(self, step: PipelineStep) -> None:
        self.steps.append(step)
    def execute(self, data: List[AnnData]) -> AnnData:
        """
        Parameters:
            List of AnnData Objects
        Return Value: An AnnData Object
        """
        processed_data = [adata.copy() for adata in data]
        for step in self.steps:
            processed_data = step.apply(processed_data)
        return processed_data