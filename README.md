# OncoAgent

## LLM-driven cancer genomics analysis platform for TCGA-BRCA

OncoAgent is an LLM-powered cancer genomics analysis application that converts natural-language questions into reproducible bioinformatics analyses of **TCGA breast cancer (TCGA-BRCA) data**.

The LLM acts as an **analysis orchestrator and biological interpreter**, while statistical analyses are performed by deterministic Python functions. Analysis results are returned as structured data that can be used both for biological interpretation and interactive visualization.

## Architecture

```text
Natural-language question
          ↓
    LLM analysis agent
          ↓
      Tool selection
          ↓
 Deterministic Python analysis
          ↓
   Structured results
       ↙       ↘
Visualization   LLM biological
                interpretation
```

This architecture separates **LLM reasoning and interpretation from statistical computation**, reducing the risk of the LLM inventing numerical results.

## Analyses

OncoAgent currently supports:

* **Differential expression** — ER-positive vs ER-negative TCGA-BRCA tumors
* **Hallmark pathway enrichment (GSEA)** — biological pathways associated with ER status
* **Overall survival analysis** — ER status and survival using Kaplan-Meier and Cox regression
* **Interactive visualization** — differential-expression volcano plots and pathway enrichment plots
* **Natural-language biological interpretation** of statistical results

## Example question

> Which genes differ between ER-positive and ER-negative breast tumors?

OncoAgent:

1. Identifies the appropriate differential-expression analysis.
2. Executes the deterministic Python analysis.
3. Applies multiple-testing correction.
4. Annotates genes with Ensembl IDs and gene symbols.
5. Returns structured statistical results.
6. Generates an interactive volcano plot.
7. Provides a biological interpretation based on the computed results.

### Example output

The differential-expression analysis compares:

* **807 ER-positive tumors**
* **237 ER-negative tumors**
* **3,796 genes**

The application returns ranked differentially expressed genes together with effect sizes, test statistics, p-values and FDR-adjusted p-values.

## Technology

* **Python**
* **OpenAI API**
* **Pandas**
* **NumPy**
* **SciPy**
* **Statsmodels**
* **GSEApy**
* **Lifelines**
* **Plotly**
* **Streamlit**

## Scientific safeguards

The application is designed to distinguish statistical computation from generative interpretation.

* Statistical analyses are performed by deterministic Python functions.
* Numerical results returned by the analysis tools are treated as the source of truth.
* Multiple-testing correction is applied to differential-expression results.
* Cohort sizes are explicitly reported for each analysis.
* Survival results report hazard ratios and confidence intervals with an explicit reference group.
* TCGA observational associations are not interpreted as causal relationships.
* The application does not make individual patient mortality predictions.

## Data

The project uses processed **TCGA-BRCA gene-expression and clinical data**.

The underlying data files are not included in this repository. See `data/README.md` for information on the required input format and data preparation.

## Project structure

```text
OncoAgent/
│
├── src/
│   ├── app.py
│   ├── llm_agent.py
│   ├── data.py
│   ├── differential_expression.py
│   ├── pathway.py
│   ├── survival.py
│   └── gene_annotation.py
│
├── data/
│   └── README.md
│
├── results/
│
├── tests/
│
├── requirements.txt
└── README.md
```

## Future development

Planned extensions include:

* Kaplan-Meier visualization directly in the Streamlit interface
* Expanded survival-analysis visualization
* Integrated multi-tool questions combining gene expression, pathways and survival
* Additional cancer genomics analyses
* Automated testing and validation of analysis tools
* Deployment as an interactive web application

## Project motivation

OncoAgent was developed to explore how large language models can be integrated with reproducible cancer genomics workflows.

Rather than allowing an LLM to perform statistical calculations directly, the system uses the LLM to determine **what analysis is required**, delegates the computation to deterministic bioinformatics functions, and then uses the resulting evidence to generate a biological interpretation.
