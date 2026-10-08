from src.pca import run_pca


def test_pca_returns_two_components(expression, clinical):

    result = run_pca(
        expression,
        clinical,
        variable_of_interest="ER_status"
    )

    assert "pca_data" in result
    assert "explained_variance" in result

    assert result["pca_data"].shape[1] == 4
    assert len(result["explained_variance"]) == 2
    
    
def test_pca_has_samples(expression, clinical):

    result = run_pca(
        expression,
        clinical,
        variable_of_interest="ER_status"
    )

    assert result["n_samples"] > 0
    assert result["n_genes"] > 0