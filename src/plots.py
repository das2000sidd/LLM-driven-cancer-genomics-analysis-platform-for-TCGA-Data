import numpy as np
import matplotlib.pyplot as plt


def create_volcano_plot(
    results,
    output_file
):

    plot_data = results.copy()
    
    effect_column = "difference"

    plot_data["minus_log10_FDR"] = (
        -np.log10(
            plot_data["FDR"].clip(
                lower=1e-300
            )
        )
    )

    significant = (
        (plot_data["FDR"] < 0.05)
        &
        (abs(plot_data["difference"]) > 1)
    )

    plt.figure(
        figsize=(9, 7)
    )

    plt.scatter(
        plot_data[effect_column],
        plot_data["minus_log10_FDR"],
        alpha=0.6
    )

    plt.scatter(
        plot_data.loc[significant, effect_column],
        plot_data.loc[significant, "minus_log10_FDR"]
    )

    plt.axvline(
        x=1,
        linestyle="--"
    )

    plt.axvline(
        x=-1,
        linestyle="--"
    )

    plt.axhline(
        y=-np.log10(0.05),
        linestyle="--"
    )

    plt.xlabel("Difference in VSD expression")
    
    plt.ylabel("-log10 FDR")
    
    plt.title("ER-positive vs ER-negative")

    plt.tight_layout()

    plt.savefig(
        output_file,
        dpi=300
    )

    plt.close()
    
import matplotlib.pyplot as plt
from pathlib import Path


def plot_pca(
    pca_data,
    variable_of_interest,
    explained_variance,
    output_path
):
    """
    Create and save a PCA scatter plot.

    Samples are colored/grouped according to
    variable_of_interest.
    """

    output_path = Path(output_path)

    # Create output directory if necessary
    output_path.parent.mkdir(
        parents=True,
        exist_ok=True
    )

    fig, ax = plt.subplots(
        figsize=(8, 6)
    )

    # Plot each group separately
    for group in pca_data[
        variable_of_interest
    ].dropna().unique():

        subset = pca_data[
            pca_data[variable_of_interest] == group
        ]

        ax.scatter(
            subset["PC1"],
            subset["PC2"],
            label=str(group),
            alpha=0.7,
            s=25
        )

    # Explained variance
    pc1 = explained_variance[0] * 100
    pc2 = explained_variance[1] * 100

    ax.set_xlabel(
        f"PC1 ({pc1:.2f}% variance)"
    )

    ax.set_ylabel(
        f"PC2 ({pc2:.2f}% variance)"
    )

    ax.set_title(
        f"PCA colored by {variable_of_interest}"
    )

    ax.legend(
        title=variable_of_interest
    )

    plt.tight_layout()

    fig.savefig(
        output_path,
        dpi=300,
        bbox_inches="tight"
    )

    plt.close(fig)

    return str(output_path)