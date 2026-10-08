import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.express as px
import streamlit as st



# Project paths
PROJECT_ROOT = Path(__file__).resolve().parent.parent
SRC_DIR = PROJECT_ROOT / "src"

if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))



# Imports from project


from data import load_clinical, load_expression

from llm_agent import (
    initialize_data,
    ask_agent,
    get_last_tool_results
)

# Data paths
CLINICAL_FILE = os.getenv(
    "ONCOAGENT_CLINICAL_FILE",
    str(
        PROJECT_ROOT
        / "data"
        / "TCGA_BRCA_clinical.csv"
    )
)

EXPRESSION_FILE = os.getenv(
    "ONCOAGENT_EXPRESSION_FILE",
    str(
        PROJECT_ROOT
        / "data"
        / "TCGA_BRCA_expression.csv"
    )
)



# Page configuration


st.set_page_config(
    page_title="OncoAgent",
    page_icon="🧬",
    layout="wide"
)



# Page title
st.title("OncoAgent")

st.markdown(
    """
    **LLM-driven cancer genomics analysis of TCGA-BRCA data**

    Ask a natural-language question about gene expression,
    biological pathways, or survival. OncoAgent selects the
    appropriate bioinformatics analysis and returns the
    statistical results together with a biological interpretation.
    """
)



# Load data
@st.cache_data
def load_data():

    clinical = load_clinical(
        CLINICAL_FILE
    )

    expression = load_expression(
        EXPRESSION_FILE
    )

    return clinical, expression


try:

    clinical, expression = load_data()

    initialize_data(
        expression,
        clinical
    )

except Exception as e:

    st.error(
        f"Failed to load TCGA data: {e}"
    )

    st.stop()



# Sidebar
with st.sidebar:

    st.header("Dataset")

    st.metric(
        "Clinical samples",
        len(clinical)
    )

    st.metric(
        "Expression samples",
        expression.shape[0]
    )

    st.metric(
        "Genes",
        expression.shape[1] - 1
    )

    st.divider()

    st.subheader("Example questions")

    st.markdown(
        """
        **Differential expression**

        • Which genes differ between ER-positive
        and ER-negative tumors?

        **Pathways**

        • What biological pathways are associated
        with ER status?

        **Survival**

        • Is ER status associated with overall survival?
        """
    )



# User question
st.subheader("Ask a cancer genomics question")

question = st.text_area(
    "Question",
    placeholder=(
        "For example: Which genes differ between "
        "ER-positive and ER-negative tumors?"
    ),
    height=100
)



# Run analysis


if st.button(
    "Analyze",
    type="primary"
):

    if not question.strip():

        st.warning(
            "Please enter a question."
        )

        st.stop()

    # Tool callback
    def show_tool(tool_name):

        tool_labels = {
            "run_differential_expression":
                "Running differential expression...",

            "run_gsea":
                "Running Hallmark pathway enrichment...",

            "run_survival_analysis":
                "Running overall survival analysis..."
        }

        message = tool_labels.get(
            tool_name,
            f"Running {tool_name}..."
        )

        st.info(message)

    # Run LLM agent
    try:

        with st.spinner(
            "OncoAgent is analyzing the question..."
        ):

            answer = ask_agent(
                question,
                tool_callback=show_tool
            )

            # Retrieve the structured results generated
            # by the tools.
            tool_results = get_last_tool_results()

    except Exception as e:

        st.error(
            f"Analysis failed: {e}"
        )

        st.stop()

    # LLM interpretation
    st.markdown(
        "### Biological interpretation"
    )

    st.markdown(answer)

    # Differential expression results
    de_result = tool_results.get(
        "run_differential_expression"
    )

    if de_result is not None:

        st.markdown(
            "### Differential expression"
        )

        # Cohort summary
        col1, col2, col3 = st.columns(3)

        with col1:

            st.metric(
                "ER-positive",
                de_result["n_er_positive"]
            )

        with col2:

            st.metric(
                "ER-negative",
                de_result["n_er_negative"]
            )

        with col3:

            st.metric(
                "Genes tested",
                de_result["n_genes_tested"]
            )

        # Top ER-positive genes
        st.markdown(
            "#### Genes associated with ER-positive tumors"
        )

        positive_genes = pd.DataFrame(
            de_result[
                "top_er_positive_by_difference"
            ]
        )

        if not positive_genes.empty:

            st.dataframe(
                positive_genes,
                use_container_width=True,
                hide_index=True
            )

        # ----------------------------------------------------
        # Top ER-negative genes
        # ----------------------------------------------------

        st.markdown(
            "#### Genes associated with ER-negative tumors"
        )

        negative_genes = pd.DataFrame(
            de_result[
                "top_er_negative_by_difference"
            ]
        )

        if not negative_genes.empty:

            st.dataframe(
                negative_genes,
                use_container_width=True,
                hide_index=True
            )

        # Volcano plot
        st.markdown(
            "### Differential expression volcano plot"
        )

        volcano = pd.DataFrame(
            de_result["volcano_data"]
        )

        if not volcano.empty:

            # Define significance
            volcano["significant"] = (
                (volcano["FDR"] < 0.05)
                &
                (volcano["difference"].abs() >= 1)
            )

            # Create plot
            fig = px.scatter(
                volcano,

                x="difference",

                y="neg_log10_fdr",

                hover_name="label",

                color="significant",

                labels={
                    "difference":
                        "Mean expression difference "
                        "(ER+ − ER−)",

                    "neg_log10_fdr":
                        "−log10(FDR)",

                    "significant":
                        "Significant"
                },

                title=(
                    "ER-positive vs ER-negative "
                    "differential expression"
                )
            )

            # Add thresholds
            fig.add_vline(
                x=1,
                line_dash="dash"
            )

            fig.add_vline(
                x=-1,
                line_dash="dash"
            )

            fig.add_hline(
                y=-np.log10(0.05),
                line_dash="dash"
            )

            # Display plot
            st.plotly_chart(
                fig,
                use_container_width=True
            )

        else:

            st.warning(
                "No volcano plot data were returned."
            )

    # GSEA results
    gsea_result = tool_results.get(
        "run_gsea"
    )

    if gsea_result is not None:

        st.markdown(
            "### Hallmark pathway enrichment"
        )

        pathways = pd.DataFrame(
            gsea_result["pathways"]
        )

        if not pathways.empty:

            # Select the most informative pathways
            pathways["FDR q-val"] = pd.to_numeric(
                pathways["FDR q-val"],
                errors="coerce"
            )

            pathways["NES"] = pd.to_numeric(
                pathways["NES"],
                errors="coerce"
            )

            pathways = (
                pathways
                .sort_values("FDR q-val")
                .head(20)
                .sort_values("NES")
            )

            # Plot
            fig = px.bar(
                pathways,

                x="NES",

                y="Term",

                orientation="h",

                hover_data=[
                    "FDR q-val"
                ],

                title=(
                    "Hallmark pathway enrichment "
                    "by ER status"
                )
            )

            fig.add_vline(
                x=0,
                line_dash="dash"
            )

            st.plotly_chart(
                fig,
                use_container_width=True
            )

            # Table
            st.markdown(
                "#### Pathway enrichment results"
            )

            st.dataframe(
                pathways,
                use_container_width=True,
                hide_index=True
            )

    # Survival results
    survival_result = tool_results.get(
        "run_survival_analysis"
    )

    if survival_result is not None:

        st.markdown(
            "### Overall survival"
        )

        # Display the structured survival result.
        # This section can be expanded once the exact
        # survival.py output structure is finalized.

        st.json(
            survival_result
        )