from types import SimpleNamespace
from unittest.mock import Mock
import pandas as pd
from src import llm_agent


def test_agent_executes_requested_tool(monkeypatch):
    # Simulate the LLM requesting GSEA.
    tool_call_response = SimpleNamespace(
        id="resp_1",
        output=[
            SimpleNamespace(
                type="function_call",
                name="run_gsea",
                arguments="{}",
                call_id="call_123"
            )
        ],
        output_text=""
    )

    # Simulate the LLM's final response after receiving GSEA results.
    final_response = SimpleNamespace(
        id="resp_2",
        output=[],
        output_text=(
            "Hallmark pathway enrichment completed. "
            "Interpret the results using the returned statistics."
        )
    )

    # Return the simulated responses in sequence.
    mock_create = Mock(
        side_effect=[tool_call_response, final_response]
    )
    monkeypatch.setattr(
        llm_agent.client.responses,
        "create",
        mock_create
    )

    # Simulate the deterministic GSEA tool.
    mock_execute = Mock(return_value={
        "analysis": "Hallmark GSEA",
        "n_pathways": 1
    })
    monkeypatch.setattr(
        llm_agent,
        "execute_tool",
        mock_execute
    )

    answer = llm_agent.ask_agent(
        "Which pathways are enriched in ER-positive compared with ER-negative breast cancer?"
    )

    # Check that the correct tool was executed.
    mock_execute.assert_called_once_with("run_gsea", {})

    # Check that the LLM was called first to select a tool,
    # then again to produce the final response.
    assert mock_create.call_count == 2

    # Check that the agent returned the final response.
    assert "pathway enrichment completed" in answer.lower()


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