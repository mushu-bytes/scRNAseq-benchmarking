import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd
import warnings
from typing import List
from anndata import AnnData
from sklearn.metrics import (
    silhouette_score,
    davies_bouldin_score,
    calinski_harabasz_score,
    jaccard_score,
)
from scib.metrics import nmi, cell_cycle, hvg_overlap, ari
from typing import List


def silhouette(dataset: AnnData, key: str = "clusters") -> float:
    """
    Parameters
        dataset: AnnData object with clusters precalculated
        key: where the clusters are located within ann data object
        random_state: for reproducibility
    Return Value: silhouette
    """
    sil = silhouette_score(dataset.X, dataset.obs[key])
    return sil


def davies(dataset: AnnData, key: str = "clusters") -> float:
    """
    Parameters
        dataset: AnnData object with clusters precalculated
        key: where the clusters are located within ann data object
        random_state: for reproducibility
    Return Value: davies
    """
    davies = davies_bouldin_score(dataset.X, dataset.obs[key])
    return davies


def calinski(dataset: AnnData, key: str = "clusters") -> float:
    """
    Parameters
        dataset: AnnData object with clusters precalculated
        key: where the clusters are located within ann data object
        random_state: for reproducibility
    Return Value: calinski
    """
    calinski = calinski_harabasz_score(dataset.X, dataset.obs[key])
    return calinski


def cell_type_stats(
    data: AnnData, raw: AnnData, data_key: str = "Type", label: str = "cell_type"
) -> List[int]:
    """
    Parameters:
        data: Integrated dataset
        raw: Unintegrated dataset. Preproccessed data concatenated
        key: column delineating datasets
        label: column containing cell type labels
    Return Value: List of scores
    """
    ari_score = ari(data, data_key, label)
    nmi_score = nmi(data, data_key, label)
    hvg_overlap_score = hvg_overlap(raw, data, data_key)
    cell_cycle_score = cell_cycle(raw, data, data_key, organism="human")
    return [ari_score, nmi_score, hvg_overlap_score, cell_cycle_score]


def jaccard(dataset: AnnData, key: str = "clusters", num_genes: int = 100) -> float:
    """
    Parameters
        dataset: AnnData object with clusters precalculated
        key: where the clusters are located within ann data object
        num_genes: number of genes to compare
    Return Value: average jaccard score across pairs of clusters
    """
    clusters = dataset.obs[key].unique()  # get clusters
    top_genes_per_cluster = {}
    for c in clusters:
        cluster_indices = np.where(dataset.obs["clusters"] == c)[
            0
        ]  # get cells of the cluster
        top_umis = (
            dataset.X[cluster_indices].mean(axis=0).argsort()[-num_genes:]
        )  # extract the top 20 genes
        top_genes_per_cluster[c] = dataset.var.iloc[top_umis].index.tolist()

    # aggregate scores
    distances = []
    for c_1 in clusters:
        for c_2 in clusters:
            if c_1 != c_2:
                # distance is 1 - jaccard
                distances.append(
                    1
                    - jaccard_score(
                        top_genes_per_cluster[c_1],
                        top_genes_per_cluster[c_2],
                        average="macro",
                    )
                )
    return np.mean(distances)


def evaluate(
    dataset: AnnData,
    key: str = "clusters",
    metrics: List[str] = ["jaccard", "silhouette", "davies", "calinski"],
    num_genes: int = 100,
) -> List[float]:
    """
    Parameters
        dataset: AnnData object with clusters precalculated
        key: where the clusters are located within ann data object
        metrics: what metrics to use
        num_genes: how many genes to compare
    Return Value: List of scores from metrics chosen
    """
    scores = []
    for m in metrics:
        if m == "jaccard":
            scores.append(jaccard(dataset, key, num_genes))
        elif m == "silhouette":
            scores.append(silhouette(dataset, key))
        elif m == "davies":
            scores.append(davies(dataset, key))
        elif m == "calinski":
            scores.append(calinski(dataset, key))
        else:
            warnings.warn(f"Invalid Metric: {m}")
    return scores
