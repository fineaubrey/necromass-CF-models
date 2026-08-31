#!/usr/bin/env python3
"""
Observed dependence diagnostics for Dataset 2 soil variables.

PURPOSE
-------
Summarize the empirical dependence among soil muramic acid (MS),
soil glucosamine (GS), and soil organic carbon (S) in the complete-case
Dataset 2 used to parameterize the uncertainty analysis.

This analysis is DESCRIPTIVE ONLY.

The empirical correlations reported here are not imposed on the Sobol
sampling design. Standard Sobol sensitivity analysis in the manuscript
varies uncertain inputs independently.

RUN FROM REPOSITORY ROOT
------------------------
python -m scripts.Python.soil_dependence_diagnostics

OUTPUTS
-------
results/tables/
    TableB2_soil_spearman_correlations.csv
    TableB2_soil_spearman_matrix.csv
"""

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from . import model_definitions as md


# =====================================================================
# Paths
# =====================================================================

OUTDIR = (
    md.REPO_ROOT
    / "results"
    / "tables"
)

OUTDIR.mkdir(
    parents=True,
    exist_ok=True,
)


# =====================================================================
# Canonical Dataset 2 subset
# =====================================================================

soil = (
    md.soils[
        [
            "MS",
            "GS",
            "S",
        ]
    ]
    .copy()
)

soil = (
    soil
    .replace(
        [np.inf, -np.inf],
        np.nan,
    )
    .dropna()
)

n = len(soil)

if n == 0:
    raise ValueError(
        "No complete soil observations available."
    )


print("=" * 78)
print("OBSERVED SOIL-VARIABLE DEPENDENCE DIAGNOSTICS")
print("=" * 78)

print()
print("Complete observations:", n)

print()
print(
    "NOTE: These correlations describe the empirical Dataset 2 "
    "joint structure."
)

print(
    "They are NOT imposed on the independent Sobol sampling design."
)

print()


# =====================================================================
# Variable labels
# =====================================================================

labels = {
    "MS": "M_S",
    "GS": "G_S",
    "S": "S",
}

descriptions = {
    "MS": "Soil muramic acid",
    "GS": "Soil glucosamine",
    "S": "Soil organic carbon",
}


# =====================================================================
# Pairwise Spearman correlations
# =====================================================================

pairs = [
    ("MS", "GS"),
    ("MS", "S"),
    ("GS", "S"),
]

rows = []

for x, y in pairs:

    rho, p_value = spearmanr(
        soil[x].to_numpy(dtype=float),
        soil[y].to_numpy(dtype=float),
        nan_policy="omit",
    )

    rows.append({
        "variable_1": labels[x],
        "variable_1_description": descriptions[x],
        "variable_2": labels[y],
        "variable_2_description": descriptions[y],
        "N": n,
        "spearman_rho": rho,
        "p_value": p_value,
    })


cor_long = pd.DataFrame(
    rows
)


# =====================================================================
# Formatted values for manuscript checking
# =====================================================================

cor_long["rho_rounded"] = (
    cor_long["spearman_rho"]
    .round(3)
)

cor_long["p_formatted"] = cor_long[
    "p_value"
].apply(
    lambda p: (
        "<0.001"
        if p < 0.001
        else f"{p:.3f}"
    )
)


# =====================================================================
# Spearman matrix
# =====================================================================

rho_matrix = (
    soil[
        [
            "MS",
            "GS",
            "S",
        ]
    ]
    .corr(
        method="spearman"
    )
    .rename(
        index=labels,
        columns=labels,
    )
)

rho_matrix = rho_matrix.round(
    3
)


# =====================================================================
# Save
# =====================================================================

long_path = (
    OUTDIR
    / "TableB2_soil_spearman_correlations.csv"
)

matrix_path = (
    OUTDIR
    / "TableB2_soil_spearman_matrix.csv"
)

cor_long.to_csv(
    long_path,
    index=False,
)

rho_matrix.to_csv(
    matrix_path,
)


# =====================================================================
# Console output
# =====================================================================

print("=" * 78)
print("PAIRWISE SPEARMAN CORRELATIONS")
print("=" * 78)

print(
    cor_long[
        [
            "variable_1",
            "variable_2",
            "N",
            "rho_rounded",
            "p_formatted",
        ]
    ]
    .rename(
        columns={
            "rho_rounded": "Spearman_rho",
            "p_formatted": "P",
        }
    )
    .to_string(
        index=False
    )
)

print()

print("=" * 78)
print("SPEARMAN CORRELATION MATRIX")
print("=" * 78)

print(
    rho_matrix.to_string()
)

print()

print("Saved:")
print(" ", long_path)
print(" ", matrix_path)

print()
print("Dependence diagnostic complete.")
