from pathlib import Path

from src.pca import run_pca
from src.plots import plot_pca


def test_pca_plot_is_created(tmp_path, expression, clinical):

    result = run_pca(
        expression,
        clinical,
        variable_of_interest="ER_status"
    )

    output = tmp_path / "pca_test.png"

    plot_pca(
        result["pca_data"],
        "ER_status",
        result["explained_variance"],
        output
    )

    assert output.exists()