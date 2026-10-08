import numpy as np
import pandas as pd
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests


def run_differential_expression(
    expression,
    metadata,
    group_column="ER_status",
    sample_column="sample_id"
):
    """
    Perform differential expression analysis between two groups.

    Parameters
    ----------
    expression : pandas.DataFrame
        Expression matrix with samples as rows and genes as columns.
        Must contain sample_id.

    metadata : pandas.DataFrame
        Clinical metadata containing sample_id and group_column.

    group_column : str
        Column defining the two groups.

    sample_column : str
        Column containing sample identifiers.
    """

    # Make sure required columns exist
    if sample_column not in expression.columns:
        raise ValueError(
            f"{sample_column} not found in expression data"
        )

    if sample_column not in metadata.columns:
        raise ValueError(
            f"{sample_column} not found in metadata"
        )

    if group_column not in metadata.columns:
        raise ValueError(
            f"{group_column} not found in metadata"
        )

    # Keep only samples with definitive ER status
    valid_groups = ["Positive", "Negative"]
    
    metadata = metadata[ metadata[group_column].isin(valid_groups) ].copy()
    
    print("ER groups after filtering:")
    
    print(metadata[group_column].value_counts())
    
    group1 = "Positive"
    
    group2 = "Negative"
    
    # Get sample IDs for each group
    group1_ids = metadata.loc[
        metadata[group_column] == group1,
        sample_column
    ]

    group2_ids = metadata.loc[
        metadata[group_column] == group2,
        sample_column
    ]

    # Restrict expression matrix to samples in each group
    group1_expression = expression[
        expression[sample_column].isin(group1_ids)
    ].copy()

    group2_expression = expression[
        expression[sample_column].isin(group2_ids)
    ].copy()

    print(f"{group1} samples: {len(group1_expression)}")
    print(f"{group2} samples: {len(group2_expression)}")

    # Identify gene columns
    excluded_columns = {sample_column}

    genes = [
        column
        for column in expression.columns
        if column not in excluded_columns
        and pd.api.types.is_numeric_dtype(expression[column])
    ]

    results = []

    for gene in genes:

        values1 = group1_expression[gene].dropna()
        values2 = group2_expression[gene].dropna()

        if len(values1) < 3 or len(values2) < 3:
            continue

        statistic, pvalue = ttest_ind(
            values1,
            values2,
            equal_var=False
        )

        mean1 = values1.mean()
        mean2 = values2.mean()

        # Expression is already transformed,
        # so calculate the difference in means directly.
        difference = mean1 - mean2

        results.append({
    "gene": gene,
    "positive_mean": mean1,
    "negative_mean": mean2,
    "difference": difference,
    "statistic": statistic,
    "pvalue": pvalue
    })

    results = pd.DataFrame(results)

    if results.empty:
        return results

    # Multiple-testing correction
    results["FDR"] = multipletests(
        results["pvalue"],
        method="fdr_bh"
    )[1]

    results = results.sort_values("FDR")

    return results