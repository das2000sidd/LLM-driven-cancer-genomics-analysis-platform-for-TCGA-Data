from src.data import load_clinical, load_expression


def test_clinical_has_sample_id():
    clinical = load_clinical("tests/test_data/clinical.csv")

    assert "sample_id" in clinical.columns


def test_expression_has_sample_id():
    expression = load_expression("tests/test_data/expression.csv")

    assert "sample_id" in expression.columns