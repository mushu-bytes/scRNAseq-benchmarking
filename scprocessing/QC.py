import scanpy as sc
import anndata as ad
from typing import List
from anndata import AnnData
from .PipelineStep import PipelineStep


class QC(PipelineStep):
    def __init__(self, method=""):
        self.method = method

    def apply(self, datasets: List[AnnData]) -> List[AnnData]:
        """
        Parameters:
            datasets: List of AnnData objects

        Return Value: A List of QC'd AnnData objects
        """
        print("Quality Control")
        QC_data = []
        for data in datasets:
            sc.pp.calculate_qc_metrics(data)
            sc.pp.filter_cells(data, min_genes=100)
            sc.pp.filter_genes(data, min_cells=3)
            sc.pp.scrublet(data, batch_key="Trial")
            QC_data.append(data[data.obs["predicted_doublet"] == False])
        return QC_data
