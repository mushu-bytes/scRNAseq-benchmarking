import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd
from typing import List
from anndata import AnnData
from sklearn.metrics import silhouette_score, davies_bouldin_score, calinski_harabasz_score, jaccard_score
from typing import List

def silhouette(dataset: AnnData, key: str = "clusters", random_state: int = 42) -> float:
    """
    Parameters
        dataset: AnnData object with clusters precalculated
        key: where the clusters are located within ann data object
        random_state: for reproducibility 
    Return Value: silhouette
    """
    sil = silhouette_score(dataset.X, dataset.obs[key], random_state=random_state)
    return sil

def davies(dataset: AnnData, key: str = "clusters", random_state: int = 42) -> float:
    """
    Parameters
        dataset: AnnData object with clusters precalculated
        key: where the clusters are located within ann data object
        random_state: for reproducibility 
    Return Value: davies
    """
    davies = davies_bouldin_score(dataset.X, dataset.obs[key], random_state=random_state)
    return davies

def calinski(dataset: AnnData, key: str = "clusters", random_state: int = 42) -> float:
    """
    Parameters
        dataset: AnnData object with clusters precalculated
        key: where the clusters are located within ann data object
        random_state: for reproducibility 
    Return Value: calinski
    """
    calinski = calinski_harabasz_score(dataset.X, dataset.obs[key], random_state=random_state)
    return calinski

def jaccard(dataset: AnnData, key: str = "clusters") -> float:
    """
    Parameters
        dataset: AnnData object with clusters precalculated
        key: where the clusters are located within ann data object
    Return Value: calinski
    """
    clusters = dataset.obs[key].unique() # get clusters
    top_genes_per_cluster = {}
    for c in clusters:
        cluster_indices = np.where(dataset.obs["clusters"] == c)[0] # get cells of the cluster
        top_umis = np.absolute(dataset.X[cluster_indices]).sum(axis=0).argsort()[-20:] # extract the top 20 genes
        top_genes_per_cluster[c] = dataset.var.iloc[top_umis].index.tolist()

    # aggregate scores
    distances = []
    for c_1 in clusters:
        for c_2 in clusters: 
            if c_1 != c_2:
                # distance is 1 - jaccard
                distances.append(1 - jaccard_score(top_genes_per_cluster[c_1], top_genes_per_cluster[c_2], average="macro"))
    return np.mean(distances)
