import pandas as pd
from pandas import DataFrame
import seaborn as sns
import matplotlib.pyplot as plt
import scanpy as sc
import anndata as ad
from anndata import AnnData


def visualize_report(
    report: DataFrame,
    data: AnnData,
    path: str,
    labels=[
        "Jaccard",
        "Adjusted Rand Index",
        "Normalized Mutual Information",
        "Runtime",
        "Number of Clusters",
    ],
) -> None:
    """
    Parameters:
        Report: Dataframe of results generated from Select Pipeline
        data: AnnData object generated from Select Pipeline
        path: Where to save figure
    Returns:
        None
    """

    res = report.reset_index().rename(
        columns={"level_0": "Normalization", "level_1": "Integration"}
    )
    fig, axs = plt.subplots(2, 3, figsize=(48, 32))
    # Remove the last subplot
    axes = [axs[0, 0], axs[0, 1], axs[0, 2], axs[1, 0], axs[1, 1]]
    for i in range(len(axes)):
        # Create grouped bar chart using Seaborn
        bar = sns.barplot(
            data=res,
            x="Normalization",
            y=report.columns[i],
            hue="Integration",
            palette="flare",
            ax=axes[i],
        )

        # Add labels
        axes[i].set_xlabel("Normalization and Integration", fontsize=28)
        axes[i].set_ylabel(labels[i], fontsize=28)
        axes[i].legend(fontsize=24)
        axes[i].tick_params(axis="both", which="major", labelsize=24)

    # adding the clustering
    umap = sc.pl.umap(
        data,
        color=["clusters"],
        palette=sc.pl.palettes.default_20,
        ax=axs[1, 2],
        legend_loc="None",
        title=f"{data.uns['Pipeline Steps']}",
    )
    fig.savefig(path)
