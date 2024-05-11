from .Pipeline import Pipeline
import scanpy as sc
import anndata as ad
from typing import List, Tuple
from anndata import AnnData
from pandas import DataFrame
import warnings
import time
import pandas as pd

from .Normalization import Normalization
from .Integration import Integration
from .QC import QC
from .metrics import evaluate


class SelectPipeline:
    def __init__(
        self,
        qc: List[str] = [],
        normalization: List[str] = ["zheng17", "seurat", "weinreb17"],
        integration: List[str] = ["merge", "harmony", "scanorama"],
        metrics: List[str] = ["jaccard", "silhouette", "davies", "calinski"],
        store_all_clusters: bool = True,
    ) -> None:
        self.qc = qc
        self.normalization = normalization
        self.integration = integration
        self.metrics = metrics
        self.store_all_clusters = store_all_clusters
        self.clusters = {}

    def search(
        self,
        data: List[AnnData],
        key_metric: str = "jaccard",
    ) -> Tuple[AnnData, DataFrame, Pipeline]:
        """
        Parameters:
            data, a list of AnnData objects
            key_metric, metric to determine best method
        Returns:
            A Tuple containing:
                Best Integrated Data
                DataFrame containing results
                a formalized pipeline object
        """
        report = {}
        # TODO: Add more qc steps and clean this code up
        if self.qc:
            qc_data = QC().apply(datasets=data)  # skipping qc for now
        else:
            qc_data = data

        for norm in self.normalization:
            try:
                # wrapping in a try catch in case failure occurs
                start_time = time.time()
                norm_data = Normalization(method=norm).apply(qc_data)
                end_time = time.time()
                norm_time = end_time - start_time
            except Exception as e:
                warnings.warn(f"Error occured while using {norm} normalization: {e}")
                warnings.warn(f"Skipping {norm} normalization")
                continue

            for integrate in self.integration:
                try:
                    start_time = time.time()
                    integration_data = Integration(method=integrate).apply(norm_data)
                    end_time = time.time()
                    integrate_time = end_time - start_time
                except Exception as e:
                    warnings.warn(
                        f"Error occured while using {norm} normalization and {integrate} integration: {e}"
                    )
                    warnings.warn(
                        f"Skipping {norm} normalization and {integrate} integration"
                    )
                    continue

                # storing clusters inside the pipeline selection:
                if self.store_all_clusters:
                    self.clusters[(norm, integrate)] = integration_data

                # adding runtime into the report
                report[(norm, integrate)] = evaluate(
                    integration_data, metrics=self.metrics
                ) + [norm_time + integrate_time]

        # generate report
        report_df = pd.DataFrame(report)
        report_df.index = self.metrics + [
            "Runtime"
        ]  # adding an additional element into the index for times
        report_df = (
            report_df.transpose()
        )  # df should be norm/integration method on index, measurements as columns
        pipeline_steps = report_df[
            key_metric
        ].idxmax()  # retrieves index of highest key_metric

        best_pipeline = Pipeline(steps=list(pipeline_steps))  # create pipeline
        return self.clusters[pipeline_steps], report_df, best_pipeline
