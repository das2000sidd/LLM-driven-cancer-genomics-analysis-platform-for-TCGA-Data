import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler


def run_pca(
    expression,
    metadata,
    variable_of_interest="ER_status",
    sample_column="sample_id",
    n_components=2
):
    """
    Run PCA on gene-expression data and annotate samples
    using a clinical variable of interest.
    """

    # Check inputs
    if sample_column not in expression.columns:
        raise ValueError(f"{sample_column} not found in expression data")

    if sample_column not in metadata.columns:
        raise ValueError(f"{sample_column} not found in metadata")

    if variable_of_interest not in metadata.columns:
        raise ValueError(
            f"{variable_of_interest} not found in metadata"
        )

    # --------------------------------------------------
    # 1. Identify gene-expression columns
    # --------------------------------------------------

    genes = [
        column
        for column in expression.columns
        if column != sample_column
        and pd.api.types.is_numeric_dtype(expression[column])
    ]

    expression_matrix = expression[
        [sample_column] + genes
    ].copy()

    # --------------------------------------------------
    # 2. Merge expression with metadata
    # --------------------------------------------------

    data = metadata[
        [sample_column, variable_of_interest]
    ].merge(
        expression_matrix,
        on=sample_column,
        how="inner"
    )

    # Remove samples without the variable of interest
    data = data.dropna(subset=[variable_of_interest])

    # --------------------------------------------------
    # 3. Extract expression matrix
    # --------------------------------------------------

    X = data[genes].copy()

    # Remove genes with missing values
    valid_genes = X.columns[X.notna().all(axis=0)]
    X = X[valid_genes]

    # --------------------------------------------------
    # 4. Standardize genes
    # --------------------------------------------------

    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    # --------------------------------------------------
    # 5. PCA
    # --------------------------------------------------

    pca = PCA(n_components=n_components)

    principal_components = pca.fit_transform(X_scaled)

    # --------------------------------------------------
    # 6. Build PCA result table
    # --------------------------------------------------

    pca_data = pd.DataFrame(
        principal_components,
        columns=[
            f"PC{i + 1}"
            for i in range(n_components)
        ]
    )

    pca_data[sample_column] = data[sample_column].values

    pca_data[variable_of_interest] = (
        data[variable_of_interest].values
    )

    # --------------------------------------------------
    # 7. Explained variance
    # --------------------------------------------------

    explained_variance = (
        pca.explained_variance_ratio_
    )

    return {
        "pca_data": pca_data,
        "explained_variance": explained_variance.tolist(),
        "n_samples": int(len(pca_data)),
        "n_genes": int(len(valid_genes)),
        "variable_of_interest": variable_of_interest
    }