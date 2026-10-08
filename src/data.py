import pandas as pd


def load_clinical(path):
    """Load and standardize TCGA clinical metadata."""

    clinical = pd.read_csv(path)

    required = ["bcr_patient_barcode", "er_status_by_ihc"]

    missing = [x for x in required if x not in clinical.columns]

    if missing:
        raise ValueError(
            f"Missing required clinical columns: {missing}"
        )

    # Keep only the columns we currently need
    clinical = clinical.copy()

    # Rename to standardized names used by OncoAgent
    clinical = clinical.rename(
        columns={
            "bcr_patient_barcode": "sample_id",
            "er_status_by_ihc": "ER_status"
        }
    )

    return clinical


def load_expression(path):
    """Load TCGA gene-expression matrix."""

    expression = pd.read_csv(path, index_col=0)

    print("Original expression shape:", expression.shape)

    # TCGA expression matrix is expected to be:
    # rows = genes
    # columns = patients
    #
    # Transpose so that:
    # rows = patients
    # columns = genes

    expression = expression.T

    # Convert index into sample_id column
    expression.index.name = "sample_id"
    expression = expression.reset_index()

    print("Transposed expression shape:", expression.shape)

    return expression


def merge_data(clinical, expression):
    """Merge clinical metadata with gene expression."""

    data = clinical.merge(
        expression,
        on="sample_id",
        how="inner"
    )

    print("Merged data shape:", data.shape)

    return data
    
def cohort_summary(clinical):

    summary = {
        "n_samples": len(clinical),
        "er_positive": (
            clinical["ER_status"] == "Positive"
        ).sum(),
        "er_negative": (
            clinical["ER_status"] == "Negative"
        ).sum(),
    }

    return summary