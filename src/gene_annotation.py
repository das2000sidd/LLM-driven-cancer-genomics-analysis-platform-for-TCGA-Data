import pandas as pd
import mygene


def prepare_gene_ids(results):
    """
    Add version-free Ensembl IDs to DE results.
    """

    results = results.copy()

    results["ensembl_id"] = (
        results["gene"]
        .astype(str)
        .str.split(".")
        .str[0]
    )

    return results


def map_ensembl_to_symbol(results):
    """
    Map Ensembl gene IDs to HGNC gene symbols.
    """

    results = results.copy()

    mg = mygene.MyGeneInfo()

    ensembl_ids = results["ensembl_id"].dropna().unique().tolist()

    annotations = mg.querymany(
        ensembl_ids,
        scopes="ensembl.gene",
        fields="symbol",
        species="human",
        as_dataframe=True,
        returnall=False
    )

    annotations = annotations.reset_index()

    # Depending on MyGene response, the Ensembl ID may be
    # stored in different columns.
    if "query" in annotations.columns:
        annotations = annotations.rename(
            columns={"query": "ensembl_id"}
        )

    annotations = annotations[
        ["ensembl_id", "symbol"]
    ].dropna()

    annotations = annotations.drop_duplicates(
        subset="ensembl_id"
    )

    results = results.merge(
        annotations,
        on="ensembl_id",
        how="left"
    )

    return results