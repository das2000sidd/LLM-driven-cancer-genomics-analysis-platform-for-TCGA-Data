import pandas as pd

from src import llm_agent


def test_initialize_data(expression, clinical):
    llm_agent.initialize_data(expression, clinical)

    assert llm_agent._expression is expression
    assert llm_agent._clinical is clinical


def test_execute_differential_expression(expression, clinical, monkeypatch):
    llm_agent.initialize_data(expression, clinical)

    def fake_annotation(results):
        results = results.copy()
        results["symbol"] = results["gene"]
        return results

    monkeypatch.setattr(
        llm_agent,
        "map_ensembl_to_symbol",
        fake_annotation
    )

    result = llm_agent.execute_tool(
        "run_differential_expression",
        {}
    )

    assert isinstance(result, dict)
    assert result["analysis"] == "ER-positive vs ER-negative differential expression"
    assert result["n_er_positive"] == 4
    assert result["n_er_negative"] == 4
    assert result["n_genes_tested"] == 5

def test_execute_gsea(expression, clinical, monkeypatch):
    llm_agent.initialize_data(expression, clinical)

    def fake_gsea(
        ranked_genes,
        gene_sets="MSigDB_Hallmark_2020",
        output_directory=None
    ):
        return pd.DataFrame({
            "Term": ["Hallmark test pathway"],
            "NES": [1.5],
            "NOM p-val": [0.01],
            "FDR q-val": [0.05]
        })

    monkeypatch.setattr(
        llm_agent,
        "run_gsea",
        fake_gsea
    )

    def fake_annotation(results):
        results = results.copy()
        results["symbol"] = results["gene"]
        return results

    monkeypatch.setattr(
        llm_agent,
        "map_ensembl_to_symbol",
        fake_annotation
    )

    result = llm_agent.execute_tool(
        "run_gsea",
        {}
    )

    assert isinstance(result, dict)
    assert result["analysis"] == "Hallmark GSEA"


def test_execute_survival(clinical):
    llm_agent.initialize_data(
        pd.DataFrame(),
        clinical
    )

    result = llm_agent.execute_tool(
        "run_survival_analysis",
        {}
    )

    assert isinstance(result, dict)
    assert "Overall survival" in result["analysis"]


def test_execute_pca(expression, clinical):
    llm_agent.initialize_data(expression, clinical)

    result = llm_agent.execute_tool(
        "run_pca",
        {"variable_of_interest": "ER_status"}
    )

    assert isinstance(result, dict)
    assert result["analysis"] == "Principal component analysis"
    assert result["n_samples"] == 8
    assert result["n_genes"] == 5
    assert "PC1_variance" in result
    assert "PC2_variance" in result
    assert "plot_path" in result