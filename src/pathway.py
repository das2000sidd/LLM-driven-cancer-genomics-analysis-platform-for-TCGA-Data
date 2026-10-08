import pandas as pd
import gseapy as gp


def prepare_ranked_gene_list(
    results,
    gene_column="symbol",
    score_column="difference"
):
    """
    Prepare ranked gene list for GSEA.

    Positive scores:
        higher expression in ER-positive tumors

    Negative scores:
        higher expression in ER-negative tumors
    """

    ranked = results[
        [gene_column, score_column]
    ].dropna()

    ranked = ranked.drop_duplicates(
        subset=gene_column
    )

    ranked = ranked.sort_values(
        score_column,
        ascending=False
    )

    ranked = ranked.set_index(
        gene_column
    )[score_column]

    return ranked


def run_gsea(
    ranked_genes,
    gene_sets="MSigDB_Hallmark_2020",
    output_directory="../results/gsea"
):
    """
    Run preranked GSEA.
    """

    ranking = ranked_genes.reset_index()

    ranking.columns = [
        "gene",
        "score"
    ]

    pre_res = gp.prerank(
        rnk=ranking,
        gene_sets=gene_sets,
        min_size=10,
        max_size=500,
        permutation_num=1000,
        outdir=output_directory,
        seed=42,
        verbose=False
    )

    results = pre_res.res2d.copy()

    return results