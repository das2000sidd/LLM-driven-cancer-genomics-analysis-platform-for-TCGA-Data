import json
import os
import numpy as np

from openai import OpenAI

from src.differential_expression import run_differential_expression
from src.pathway import prepare_ranked_gene_list, run_gsea
from src.gene_annotation import prepare_gene_ids, map_ensembl_to_symbol
from src.survival import run_survival_analysis
from src.pca import run_pca
from src.plots import plot_pca

client = OpenAI(
    api_key=os.environ["OPENAI_API_KEY"]
)


SYSTEM_PROMPT = """
You are an expert cancer bioinformatics scientist.

You are an analysis agent working with TCGA-BRCA gene-expression data.

Your job is to answer scientific questions by using the available
bioinformatics tools rather than inventing results.

Important rules:

1. Use tools whenever the question requires analysis of the TCGA data.
2. Never invent numerical results.
3. Treat the output of the Python analysis tools as the source of truth.
4. Clearly distinguish ER-positive from ER-negative biology.
5. Do not claim causality from observational TCGA data.
6. Explain statistical results in biologically meaningful terms.
7. Mention important limitations when appropriate.

8. For GSEA:
   - positive NES means enrichment toward ER-positive tumors
   - negative NES means enrichment toward ER-negative tumors.

9. For survival analyses:
   - Use the survival analysis tool when the question concerns
     overall survival, mortality, prognosis, or survival.
   - Report the sample size and number of observed deaths.
   - Clearly identify the reference group for the hazard ratio.
   - Interpret hazard ratios in the correct direction.
   - Treat survival results as associations, not causal effects.
   - Do not make individual patient mortality predictions.
   - Do not claim that ER status determines survival.
   - Mention important limitations of observational TCGA survival data.

10. If the available tools cannot answer a question:
    - Say explicitly that the analysis is not currently available.
    - Do not invent results or perform calculations mentally.
    - Explain what analysis or data would be required.

11. When combining results from multiple tools:
    - Treat each tool's reported cohort size as authoritative
      for that analysis.
    - Do not reuse or substitute sample sizes from one analysis
      when describing another analysis.
    - Clearly state when different analyses use different numbers
      of patients because of data availability or filtering.
"""


# ============================================================
# Global data
# ============================================================

_expression = None
_clinical = None

# Stores structured results from tools used during the
# most recent agent request.
_last_tool_results = {}


def initialize_data(expression, clinical):
    """
    Store the loaded TCGA data for use by the analysis tools.
    """

    global _expression
    global _clinical

    _expression = expression
    _clinical = clinical

# ============================================================
# Tool 1: PCA plot
# ============================================================
def pca_tool(variable_of_interest="ER_status"):

    if _expression is None or _clinical is None:
        raise RuntimeError(
            "Data have not been initialized."
        )

    # -----------------------------------------
    # Run PCA
    # -----------------------------------------

    result = run_pca(
        _expression,
        _clinical,
        variable_of_interest=variable_of_interest
    )

    # -----------------------------------------
    # Create output path
    # -----------------------------------------

    output_path = (
        f"results/pca_"
        f"{variable_of_interest}.png"
    )

    # -----------------------------------------
    # Create PCA plot
    # -----------------------------------------

    plot_path = plot_pca(
        result["pca_data"],
        variable_of_interest,
        result["explained_variance"],
        output_path
    )

    # -----------------------------------------
    # Prepare structured result for LLM
    # -----------------------------------------

    return {
        "analysis": "Principal component analysis",

        "variable_of_interest":
            variable_of_interest,

        "n_samples":
            int(result["n_samples"]),

        "n_genes":
            int(result["n_genes"]),

        "PC1_variance":
            float(
                result["explained_variance"][0]
            ),

        "PC2_variance":
            float(
                result["explained_variance"][1]
            ),

        "plot_path":
            plot_path
    }


# ============================================================
# Tool 2: Differential expression
# ============================================================

def differential_expression_tool():
    """
    Run differential expression comparing ER-positive
    and ER-negative TCGA-BRCA tumors.

    Returns statistical results for the LLM together
    with structured data for Streamlit visualization.
    """

    if _expression is None or _clinical is None:
        raise RuntimeError(
            "Data have not been initialized."
        )

    # --------------------------------------------------------
    # Run differential expression
    # --------------------------------------------------------

    results = run_differential_expression(
        _expression,
        _clinical
    )

    # --------------------------------------------------------
    # Determine actual expression-matched cohort sizes
    # --------------------------------------------------------

    valid_clinical = _clinical[
        _clinical["ER_status"].isin(
            ["Positive", "Negative"]
        )
    ].copy()

    er_positive_ids = set(
        valid_clinical.loc[
            valid_clinical["ER_status"] == "Positive",
            "sample_id"
        ]
    )

    er_negative_ids = set(
        valid_clinical.loc[
            valid_clinical["ER_status"] == "Negative",
            "sample_id"
        ]
    )

    expression_ids = set(
        _expression["sample_id"]
    )

    n_er_positive = len(
        er_positive_ids.intersection(expression_ids)
    )

    n_er_negative = len(
        er_negative_ids.intersection(expression_ids)
    )

    # --------------------------------------------------------
    # Debug information
    # --------------------------------------------------------

    print("\nDEBUG: differential expression results")
    print("Type:", type(results))
    print("Shape:", results.shape)
    print("Columns:", results.columns.tolist())
    print(results.head())

    print("\nDEBUG: cohort")
    print("ER-positive:", n_er_positive)
    print("ER-negative:", n_er_negative)
    print("Genes tested:", len(results))

    # --------------------------------------------------------
    # Add Ensembl IDs and gene symbols
    # --------------------------------------------------------

    results = prepare_gene_ids(results)
    results = map_ensembl_to_symbol(results)

    # --------------------------------------------------------
    # Prepare plotting data
    # --------------------------------------------------------

    plot_results = results.copy()

    plot_results["label"] = (
        plot_results["symbol"]
        .fillna(plot_results["gene"])
        .astype(str)
    )

    plot_results["neg_log10_fdr"] = -np.log10(
        plot_results["FDR"].clip(lower=1e-300)
    )

    volcano_data = plot_results[
        [
            "label",
            "difference",
            "FDR",
            "neg_log10_fdr"
        ]
    ].to_dict(orient="records")

    # --------------------------------------------------------
    # Select top genes
    # --------------------------------------------------------

    top_by_fdr = (
        results
        .sort_values("FDR")
        .head(15)
    )

    top_er_positive = (
        results
        .sort_values(
            "difference",
            ascending=False
        )
        .head(15)
    )

    top_er_negative = (
        results
        .sort_values(
            "difference",
            ascending=True
        )
        .head(15)
    )

    output_columns = [
        "symbol",
        "gene",
        "ensembl_id",
        "positive_mean",
        "negative_mean",
        "difference",
        "statistic",
        "pvalue",
        "FDR"
    ]

    # --------------------------------------------------------
    # Return structured result
    # --------------------------------------------------------

    return {
        "analysis": (
            "ER-positive vs ER-negative "
            "differential expression"
        ),

        "n_er_positive": int(
            n_er_positive
        ),

        "n_er_negative": int(
            n_er_negative
        ),

        "n_genes_tested": int(
            len(results)
        ),

        "top_by_fdr": (
            top_by_fdr[output_columns]
            .to_dict(orient="records")
        ),

        "top_er_positive_by_difference": (
            top_er_positive[output_columns]
            .to_dict(orient="records")
        ),

        "top_er_negative_by_difference": (
            top_er_negative[output_columns]
            .to_dict(orient="records")
        ),

        "volcano_data": volcano_data
    }

# ============================================================
# Tool 3: GSEA
# ============================================================

def gsea_tool():
    """
    Run differential expression, gene annotation,
    and Hallmark GSEA.
    """

    if _expression is None or _clinical is None:
        raise RuntimeError(
            "Data have not been initialized."
        )

    # --------------------------------------------------------
    # Differential expression
    # --------------------------------------------------------

    results = run_differential_expression(
        _expression,
        _clinical
    )

    # --------------------------------------------------------
    # Ensembl ID preparation
    # --------------------------------------------------------

    results = prepare_gene_ids(results)

    # --------------------------------------------------------
    # Ensembl -> gene symbol
    # --------------------------------------------------------

    results = map_ensembl_to_symbol(results)

    # --------------------------------------------------------
    # Create ranked gene list
    # --------------------------------------------------------

    ranked_genes = prepare_ranked_gene_list(
        results,
        gene_column="symbol",
        score_column="statistic"
    )

    # --------------------------------------------------------
    # Run Hallmark GSEA
    # --------------------------------------------------------

    gsea_results = run_gsea(
        ranked_genes,
        gene_sets="MSigDB_Hallmark_2020",
        output_directory="results/gsea"
    )

    # --------------------------------------------------------
    # Select structured output
    # --------------------------------------------------------

    output = gsea_results[
        [
            "Term",
            "NES",
            "NOM p-val",
            "FDR q-val"
        ]
    ].copy()

    return {
        "analysis": "Hallmark GSEA",

        "comparison": (
            "ER-positive vs ER-negative"
        ),

        "ranking": (
            "Welch t-statistic"
        ),

        "n_genes_ranked": int(
            len(ranked_genes)
        ),

        "pathways": output.to_dict(
            orient="records"
        )
    }

# ============================================================
# Tool 4: Survival analysis
# ============================================================

def survival_analysis_tool():
    """
    Run overall survival analysis comparing
    ER-positive and ER-negative patients.
    """

    if _clinical is None:
        raise RuntimeError(
            "Clinical data have not been initialized."
        )

    result = run_survival_analysis(
        clinical=_clinical
    )

    return result


# ============================================================
# Tool definitions supplied to the LLM
# ============================================================

TOOLS = [

    {
        "type": "function",
        "name": "run_differential_expression",
        "description": (
            "Run differential expression comparing "
            "ER-positive and ER-negative TCGA-BRCA "
            "tumors. Use this when the user asks "
            "which genes differ between the groups."
        ),
        "parameters": {
            "type": "object",
            "properties": {},
            "required": []
        }
    },

    {
        "type": "function",
        "name": "run_gsea",
        "description": (
            "Run Hallmark gene-set enrichment analysis "
            "comparing ER-positive and ER-negative "
            "TCGA-BRCA tumors. Use this when the user "
            "asks about biological pathways, signaling "
            "programs, or pathway enrichment."
        ),
        "parameters": {
            "type": "object",
            "properties": {},
            "required": []
        }
    },

    {
        "type": "function",
        "name": "run_survival_analysis",
        "description": (
            "Run an unadjusted overall survival analysis "
            "comparing ER-positive and ER-negative "
            "TCGA-BRCA patients. Returns sample sizes, "
            "number of deaths, Kaplan-Meier median survival, "
            "and a univariable Cox proportional hazards "
            "model including hazard ratio, 95% confidence "
            "interval, and p-value. The ER-negative group "
            "is the reference group for the Cox model. "
            "Use this when the user asks whether ER status "
            "is associated with overall survival, mortality, "
            "prognosis, or survival. This analysis is "
            "unadjusted and should not be interpreted as "
            "an independent prognostic effect or as causal."
        ),
        "parameters": {
            "type": "object",
            "properties": {},
            "required": []
        }
    },

    {
        "type": "function",
        "name": "run_pca",
        "description": (
            "Run principal component analysis on "
            "TCGA-BRCA gene-expression data and create "
            "a PCA plot annotated by a clinical variable. "
            "Use this when the user asks for PCA, sample "
            "clustering, expression-space separation, or "
            "a PCA visualization."
        ),
        "parameters": {
            "type": "object",
            "properties": {
                "variable_of_interest": {
                    "type": "string",
                    "description": (
                        "Clinical variable used to group "
                        "samples in the PCA plot, such as "
                        "ER_status."
                    )
                }
            },
            "required": [
                "variable_of_interest"
            ]
        }
    }

]

# ============================================================
# Tool dispatcher
# ============================================================

def execute_tool(tool_name, arguments):
    """
    Execute the requested bioinformatics tool.
    """

    if tool_name == "run_differential_expression":
        return differential_expression_tool()

    elif tool_name == "run_gsea":
        return gsea_tool()

    elif tool_name == "run_survival_analysis":
        return survival_analysis_tool()

    elif tool_name == "run_pca":
        return pca_tool(
            variable_of_interest=arguments.get(
                "variable_of_interest",
                "ER_status"
            )
        )

    else:
        raise ValueError(f"Unknown tool: {tool_name}")

# ============================================================
# LLM Agent
# ============================================================

def ask_agent(question, tool_callback=None):
    """
    Ask the LLM a question and allow it to call
    bioinformatics tools.

    Parameters
    ----------
    question : str
        Natural-language cancer genomics question.

    tool_callback : callable, optional
        Function called whenever the agent invokes
        a bioinformatics tool.

    Returns
    -------
    str
        Final LLM-generated answer.
    """

    # Clear results from the previous question.
    _last_tool_results.clear()

    # --------------------------------------------------------
    # Initial LLM request
    # --------------------------------------------------------

    response = client.responses.create(
        model="gpt-5.6",

        instructions=SYSTEM_PROMPT,

        tools=TOOLS,

        input=[
            {
                "role": "user",
                "content": question
            }
        ]
    )

    # ========================================================
    # Tool-calling loop
    # ========================================================

    while True:

        tool_calls = [
            item
            for item in response.output
            if item.type == "function_call"
        ]

        # No more tools required
        if not tool_calls:
            break

        tool_outputs = []

        # ====================================================
        # Execute requested tools
        # ====================================================

        for tool_call in tool_calls:

            if tool_callback is not None:

                tool_callback(
                    tool_call.name
                )

            arguments = json.loads(
                tool_call.arguments
            )

            result = execute_tool(
                tool_call.name,
                arguments
            )

            # Store structured result so that Streamlit
            # can use it for tables and visualizations.
            _last_tool_results[
                tool_call.name
            ] = result

            tool_outputs.append(
                {
                    "type": "function_call_output",

                    "call_id": (
                        tool_call.call_id
                    ),

                    "output": json.dumps(
                        result
                    )
                }
            )

        # ====================================================
        # Send tool results back to the LLM
        # ====================================================

        response = client.responses.create(
            model="gpt-5.6",

            instructions=SYSTEM_PROMPT,

            tools=TOOLS,

            previous_response_id=response.id,

            input=tool_outputs
        )

    # ========================================================
    # Return final biological interpretation
    # ========================================================

    return response.output_text


# ============================================================
# Structured results accessor
# ============================================================

def get_last_tool_results():
    """
    Return the structured results produced by the
    most recent agent request.
    """

    return _last_tool_results