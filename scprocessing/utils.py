from typing import List
from anndata import AnnData
import scanpy as sc


def splitAD(dataset: AnnData, key: str) -> List[AnnData]:
    """
    Parameters:
        dataset: AnnData object
        key: observation field to split the dataset into
    Return:
        List of anndata split by key
    """
    unique_values = dataset.obs[key].unique()
    subsets = []
    
    # Iterate over unique values and create subsets
    for value in unique_values:
        subsets.append(dataset[dataset.obs[key] == value].copy())
    return subsets

def read_single_cell_data(filepath):
    try:
        # Try to read as 10X Genomics HDF5
        return sc.read_10x_h5(filepath)
    except:
        pass
    
    try:
        # Try to read as H5AD
        return sc.read_h5ad(filepath)
    except:
        pass


    
    raise ValueError("Unsupported file format or file could not be read")