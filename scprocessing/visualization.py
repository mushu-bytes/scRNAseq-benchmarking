import pandas as pd
from pandas import DataFrame
import seaborn as sns
import matplotlib.pyplot as plt
import scanpy as sc
import anndata as ad
from anndata import AnnData

def visualize_report(report: DataFrame, data: AnnData, path: str) -> None:
    """
    Parameters:
        Report: Dataframe of results generated from Select Pipeline
        data: AnnData object generated from Select Pipeline
        path: Where to save figure
    Returns:
        None
    """

    res = report.reset_index().rename(columns={"level_0": "Normalization", "level_1": "Integration"})
    
    fig, axs = plt.subplots(2, 3, figsize=(48, 24))
    # Remove the last subplot
    axes = [axs[0, 0], axs[0, 1], axs[0, 2], axs[1, 0], axs[1, 1]]
    for i in range(len(axes)):
        # Create grouped bar chart using Seaborn
        sns.barplot(data=res, x='Normalization', y=report.columns[i], hue='Integration', palette='flare', ax=axes[i])
        
        # Add labels and title
        axes[i].set_xlabel('Normalization+Integration')
        axes[i].set_ylabel(report.columns[i])
        axes[i].legend()
    sc.pl.umap(data, color=["clusters"], palette=sc.pl.palettes.default_20, ax=axs[1, 2], legend_loc="None")
    fig.savefig(path)
