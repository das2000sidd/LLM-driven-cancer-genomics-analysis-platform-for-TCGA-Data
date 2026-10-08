import numpy as np
import pandas as pd
import pytest


@pytest.fixture
def clinical():
    return pd.DataFrame({
        "sample_id": [
            "S1", "S2", "S3", "S4",
            "S5", "S6", "S7", "S8"
        ],
        "ER_status": [
            "Positive", "Positive", "Positive", "Positive",
            "Negative", "Negative", "Negative", "Negative"
        ],
        "vital_status": [
            "Alive", "Dead", "Alive", "Dead",
            "Dead", "Dead", "Alive", "Dead"
        ],
        "death_days_to": [
            np.nan, 200, np.nan, 400,
            50, 150, np.nan, 350
        ],
        "last_contact_days_to": [
            100, np.nan, 300, np.nan,
            np.nan, np.nan, 250, np.nan
        ]
    })


@pytest.fixture
def expression():
    return pd.DataFrame({
        "sample_id": [
            "S1", "S2", "S3", "S4",
            "S5", "S6", "S7", "S8"
        ],
        "GENE1": [5.0, 5.2, 5.1, 5.3, 2.0, 2.2, 2.1, 2.3],
        "GENE2": [3.0, 3.1, 3.2, 3.3, 4.0, 4.1, 4.2, 4.3],
        "GENE3": [6.0, 6.1, 6.2, 6.3, 6.0, 6.1, 6.2, 6.3],
        "GENE4": [1.0, 1.2, 1.1, 1.3, 2.0, 2.2, 2.1, 2.3],
        "GENE5": [7.0, 7.1, 7.2, 7.3, 5.0, 5.1, 5.2, 5.3],
    })


@pytest.fixture
def de_results(expression, clinical):
    from src.differential_expression import run_differential_expression

    results = run_differential_expression(
        expression,
        clinical
    )

    results["symbol"] = results["gene"]

    return results
