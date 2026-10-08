from src.differential_expression import run_differential_expression


def test_de_returns_expected_columns(expression, clinical):

    result = run_differential_expression(
        expression,
        clinical
    )

    expected = [
        "gene",
        "positive_mean",
        "negative_mean",
        "difference",
        "statistic",
        "pvalue",
        "FDR"
    ]

    for column in expected:
        assert column in result.columns