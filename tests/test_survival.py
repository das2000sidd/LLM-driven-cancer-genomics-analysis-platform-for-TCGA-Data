from src.survival import run_survival_analysis


def test_survival_returns_result(clinical):

    result = run_survival_analysis(clinical)

    assert result is not None