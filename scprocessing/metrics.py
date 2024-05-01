from sklearn.metrics import silhouette_score, davies_bouldin_score, calinski_harabasz_score

def silhouette(dataset: AnnData, key: str = "clusters", random_state: int = 42) -> tuple:
    """
    Parameters
        dataset: AnnData object with clusters precalculated
        key: where the clusters are located within ann data object
        random_state: for reproducibility 
    Return Value: silhouette
    """
    sil = silhouette_score(dataset.X, dataset.obs[key], random_state=random_state)
    return sil

def davies(dataset: AnnData, key: str = "clusters", random_state: int = 42) -> tuple:
    """
    Parameters
        dataset: AnnData object with clusters precalculated
        key: where the clusters are located within ann data object
        random_state: for reproducibility 
    Return Value: davies
    """
    davies = davies_bouldin_score(dataset.X, dataset.obs[key], random_state=random_state)
    return davies

def calinksi(dataset: AnnData, key: str = "clusters", random_state: int = 42) -> tuple:
    """
    Parameters
        dataset: AnnData object with clusters precalculated
        key: where the clusters are located within ann data object
        random_state: for reproducibility 
    Return Value: calinski
    """
    calinski = calinski_harabasz_score(dataset.X, dataset.obs[key], random_state=random_state)
    return calinski

