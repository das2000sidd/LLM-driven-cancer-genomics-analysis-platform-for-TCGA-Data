import pandas as pd
from lifelines import KaplanMeierFitter
from lifelines import CoxPHFitter


def prepare_survival_data(
    clinical,
    group_column="ER_status",
    positive_label="Positive",
    negative_label="Negative"
):
    """
    Construct overall survival data from TCGA clinical metadata.

    Dead patients:
        OS_time = death_days_to
        OS_event = 1

    Alive patients:
        OS_time = last_contact_days_to
        OS_event = 0
    """

    required = [
        "sample_id",
        group_column,
        "vital_status",
        "death_days_to",
        "last_contact_days_to"
    ]

    missing = [
        column
        for column in required
        if column not in clinical.columns
    ]

    if missing:
        raise ValueError(
            f"Missing survival columns: {missing}"
        )

    data = clinical[
        [
            "sample_id",
            group_column,
            "vital_status",
            "death_days_to",
            "last_contact_days_to"
        ]
    ].copy()

    # Keep ER-positive and ER-negative patients only
    data = data[
        data[group_column].isin(
            [positive_label, negative_label]
        )
    ].copy()

    # Convert survival fields to numeric
    data["death_days_to"] = pd.to_numeric(
        data["death_days_to"],
        errors="coerce"
    )

    data["last_contact_days_to"] = pd.to_numeric(
        data["last_contact_days_to"],
        errors="coerce"
    )

    # --------------------------------------------------
    # Construct overall survival endpoint
    # --------------------------------------------------

    data["OS_event"] = (
        data["vital_status"]
        .astype(str)
        .str.lower()
        .eq("dead")
        .astype(int)
    )

    data["OS_time"] = data["last_contact_days_to"]

    dead = data["OS_event"] == 1

    data.loc[dead, "OS_time"] = (
        data.loc[dead, "death_days_to"]
    )

    # Remove patients without a usable survival time
    data = data.dropna(
        subset=["OS_time"]
    )

    # Remove invalid survival times
    data = data[
        data["OS_time"] >= 0
    ].copy()

    return data


def run_survival_analysis(
    clinical,
    group_column="ER_status",
    positive_label="Positive",
    negative_label="Negative"
):
    data = prepare_survival_data(
        clinical=clinical,
        group_column=group_column,
        positive_label=positive_label,
        negative_label=negative_label
    )

    # ---------------------------------------------------------
    # Sample sizes
    # ---------------------------------------------------------

    n_total = len(data)

    n_positive = (
        data[group_column] == positive_label
    ).sum()

    n_negative = (
        data[group_column] == negative_label
    ).sum()

    # ---------------------------------------------------------
    # Number of deaths
    # ---------------------------------------------------------

    n_deaths_total = int(data["OS_event"].sum())

    n_deaths_positive = int(
        data.loc[
            data[group_column] == positive_label,
            "OS_event"
        ].sum()
    )

    n_deaths_negative = int(
        data.loc[
            data[group_column] == negative_label,
            "OS_event"
        ].sum()
    )

    # ---------------------------------------------------------
    # Kaplan-Meier analysis
    # ---------------------------------------------------------

    km_positive = KaplanMeierFitter()
    km_negative = KaplanMeierFitter()

    positive = data[
        data[group_column] == positive_label
    ]

    negative = data[
        data[group_column] == negative_label
    ]

    km_positive.fit(
        durations=positive["OS_time"],
        event_observed=positive["OS_event"],
        label="ER-positive"
    )

    km_negative.fit(
        durations=negative["OS_time"],
        event_observed=negative["OS_event"],
        label="ER-negative"
    )

    median_positive = km_positive.median_survival_time_
    median_negative = km_negative.median_survival_time_

    # Convert infinite median survival to None
    if pd.isna(median_positive) or median_positive == float("inf"):
        median_positive = None
    else:
        median_positive = float(median_positive)

    if pd.isna(median_negative) or median_negative == float("inf"):
        median_negative = None
    else:
        median_negative = float(median_negative)

    # ---------------------------------------------------------
    # Univariable Cox proportional hazards model
    # ---------------------------------------------------------

    cox_data = data[
        ["OS_time", "OS_event", group_column]
    ].copy()

    cox_data["ER_positive"] = (
        cox_data[group_column] == positive_label
    ).astype(int)

    cox_data = cox_data[
        ["OS_time", "OS_event", "ER_positive"]
    ]

    cph = CoxPHFitter()

    cph.fit(
        cox_data,
        duration_col="OS_time",
        event_col="OS_event"
    )

    hazard_ratio = float(
        cph.hazard_ratios_["ER_positive"]
    )

    ci_lower = float(
        cph.summary.loc[
            "ER_positive",
            "exp(coef) lower 95%"
        ]
    )

    ci_upper = float(
        cph.summary.loc[
            "ER_positive",
            "exp(coef) upper 95%"
        ]
    )

    p_value = float(
        cph.summary.loc[
            "ER_positive",
            "p"
        ]
    )

    # ---------------------------------------------------------
    # Return structured result
    # ---------------------------------------------------------

    return {
        "analysis": (
            "Overall survival comparison between "
            "ER-positive and ER-negative "
            "TCGA-BRCA patients"
        ),

        "endpoint": {
            "name": "Overall survival",
            "time_variable": "OS_time",
            "event_variable": "OS_event",
            "death_definition": (
                "vital_status == Dead"
            ),
            "censoring_definition": (
                "vital_status == Alive; time measured "
                "using last_contact_days_to"
            )
        },

        "model": (
            "Univariable Cox proportional hazards model"
        ),

        "adjusted": False,

        "reference_group": negative_label,

        "sample_sizes": {
            "total": int(n_total),
            "ER_positive": int(n_positive),
            "ER_negative": int(n_negative)
        },

        "events": {
            "total_deaths": n_deaths_total,
            "ER_positive_deaths": n_deaths_positive,
            "ER_negative_deaths": n_deaths_negative
        },

        "kaplan_meier": {
            "median_survival_days_ER_positive": (
                median_positive
            ),
            "median_survival_days_ER_negative": (
                median_negative
            )
        },

        "cox_model": {
            "model_type": (
                "Univariable Cox proportional hazards"
            ),
            "comparison": (
                "ER-positive vs ER-negative"
            ),
            "reference_group": (
                negative_label
            ),
            "hazard_ratio": hazard_ratio,
            "ci_lower_95": ci_lower,
            "ci_upper_95": ci_upper,
            "p_value": p_value
        }
    }