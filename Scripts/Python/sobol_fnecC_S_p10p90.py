#!/usr/bin/env python3
"""
Robustness analysis: restrict SOC (S) to its empirical P10-P90 range.

PURPOSE
-------
The primary fnecC Sobol analyses allow:

    MS -> empirical P10-P90
    GS -> empirical P10-P90
    S  -> full observed range

Because S is the denominator of fnecC and spans a substantially broader
range than MS and GS, this robustness analysis asks whether the main
sensitivity conclusions depend on the selected S range.

This script repeats the raw, unbounded Sobol analyses for:

1. Five-input composite fnecC
2. Ten-input trait-resolved fnecC

under two S designs:

    primary:
        S uses the same full-range marginal as the canonical analysis

    S_P10_P90:
        S is sampled from the fitted gamma distribution truncated to
        the empirical P10-P90 range of Dataset 2

All other parameter distributions, model equations, Sobol coordinates,
sample sizes, and seeds are unchanged.

The original and restricted-S analyses use the IDENTICAL Sobol unit-
hypercube coordinates within each formulation. Thus, the comparison
isolates the effect of changing the sampled range of S.

RUN FROM REPOSITORY ROOT
------------------------
python -m scripts.Python.sobol_fnecC_S_p10p90

OUTPUTS
-------
data/derived/gsa/
    sobol_fnecC5_S_range_robustness.csv
    sobol_fnecC10_S_range_robustness.csv
    sobol_fnecC_S_range_robustness_all.csv
    sobol_fnecC_S_range_robustness_compact.csv
    sobol_fnecC_S_range_robustness_diagnostics.csv
"""

from pathlib import Path

import numpy as np
import pandas as pd

from SALib.sample.sobol import sample as sobol_sample
from SALib.analyze.sobol import analyze as sobol_analyze

from . import model_definitions as md


# =====================================================================
# Settings
# =====================================================================

SEED = 1234
NUM_RESAMPLES = 100

N_BASE_COMPOSITE = 4096
N_BASE_TRAIT = 8192

GSA_DIR = (
    md.REPO_ROOT
    / "data"
    / "derived"
    / "gsa"
)

LHS_DIR = (
    md.REPO_ROOT
    / "data"
    / "derived"
    / "lhs"
)

GSA_DIR.mkdir(
    parents=True,
    exist_ok=True,
)


# =====================================================================
# Helpers
# =====================================================================

def make_invq(values):
    """
    Construct an empirical inverse-CDF mapping from finite values.

    This matches the approach used in the canonical Sobol scripts.
    """
    x = np.asarray(
        values,
        dtype=float,
    )

    x = x[
        np.isfinite(x)
    ]

    x = np.sort(x)

    if x.size < 50:
        raise ValueError(
            "Too few finite observations to construct inverse CDF: "
            f"n={x.size}"
        )

    u_grid = np.linspace(
        0.0,
        1.0,
        x.size,
        endpoint=True,
    )

    def invq(u):

        u = np.asarray(
            u,
            dtype=float,
        )

        u = np.clip(
            u,
            1e-12,
            1.0 - 1e-12,
        )

        return np.interp(
            u,
            u_grid,
            x,
        )

    return invq


def sobol_summary(
    problem,
    y,
    parameter_names,
    seed=SEED,
):
    """
    Calculate first-order and total-order Sobol indices.
    """

    Si = sobol_analyze(
        problem,
        np.asarray(
            y,
            dtype=float,
        ),
        calc_second_order=False,
        num_resamples=NUM_RESAMPLES,
        conf_level=0.95,
        print_to_console=False,
        seed=seed,
    )

    out = pd.DataFrame({
        "parameter": parameter_names,
        "S1": Si["S1"],
        "S1_conf": Si["S1_conf"],
        "ST": Si["ST"],
        "ST_conf": Si["ST_conf"],
    })

    out["S1_rank"] = out["S1"].rank(
        ascending=False,
        method="min",
    ).astype(int)

    out["ST_rank"] = out["ST"].rank(
        ascending=False,
        method="min",
    ).astype(int)

    return out


def qS_p10p90(u):
    """
    Sample SOC from the same fitted gamma distribution used in the
    canonical model, but truncate it to the EMPIRICAL Dataset 2
    P10-P90 range.
    """

    return md.q_truncated(
        u,
        md._s_dist,
        S_P10,
        S_P90,
    )


# =====================================================================
# Calculate empirical SOC bounds from Dataset 2
# =====================================================================

S_OBS = (
    md.soils["S"]
    .to_numpy(dtype=float)
)

S_P10, S_P90 = np.quantile(
    S_OBS,
    [0.10, 0.90],
)

S_P10 = float(S_P10)
S_P90 = float(S_P90)

S_MIN = float(
    np.min(S_OBS)
)

S_MAX = float(
    np.max(S_OBS)
)


print("=" * 78)
print("SOC RANGE ROBUSTNESS ANALYSIS")
print("=" * 78)

print()
print("Dataset 2 SOC bounds:")
print(
    f"  observed min-max : "
    f"{S_MIN:.6g} - {S_MAX:.6g} mg C g^-1 soil"
)

print(
    f"  empirical P10-P90: "
    f"{S_P10:.6g} - {S_P90:.6g} mg C g^-1 soil"
)

print()


# =====================================================================
# 1. COMPOSITE MODEL
# =====================================================================

def run_composite():

    print("=" * 78)
    print("COMPOSITE fnecC")
    print("=" * 78)

    lhs_path = (
        LHS_DIR
        / "fnecC5_lhs.csv"
    )

    if not lhs_path.exists():
        raise FileNotFoundError(
            f"Missing canonical LHS file: {lhs_path}"
        )

    lhs = pd.read_csv(
        lhs_path
    )

    parameter_names = [
        "MS",
        "GS",
        "S",
        "CFB",
        "CFF",
    ]

    missing = [
        x
        for x in parameter_names
        if x not in lhs.columns
    ]

    if missing:
        raise ValueError(
            f"{lhs_path} missing required columns: {missing}"
        )

    invq = {
        x: make_invq(
            lhs[x].to_numpy()
        )
        for x in parameter_names
    }

    unit_names = [
        "u_MS",
        "u_GS",
        "u_S",
        "u_CFB",
        "u_CFF",
    ]

    problem = {
        "num_vars": 5,
        "names": unit_names,
        "bounds": [[0.0, 1.0]] * 5,
    }

    U = sobol_sample(
        problem,
        N=N_BASE_COMPOSITE,
        calc_second_order=False,
        scramble=True,
        seed=SEED,
    )

    expected_rows = (
        N_BASE_COMPOSITE
        * (5 + 2)
    )

    if U.shape != (
        expected_rows,
        5,
    ):
        raise RuntimeError(
            "Unexpected composite Sobol matrix dimensions."
        )

    # ---------------------------------------------------------------
    # All parameters except the S robustness manipulation
    # ---------------------------------------------------------------

    MS = invq["MS"](
        U[:, 0]
    )

    GS = invq["GS"](
        U[:, 1]
    )

    # Primary S distribution exactly as used by canonical GSA.
    S_primary = invq["S"](
        U[:, 2]
    )

    # Robustness S distribution:
    # same fitted gamma distribution but restricted to empirical P10-P90.
    S_restricted = qS_p10p90(
        U[:, 2]
    )

    CFB = invq["CFB"](
        U[:, 3]
    )

    CFF = invq["CFF"](
        U[:, 4]
    )

    # ---------------------------------------------------------------
    # Primary raw response
    # ---------------------------------------------------------------

    y_primary = md.fnecC_composite_raw(
        MS=MS,
        GS=GS,
        S=S_primary,
        CFB=CFB,
        CFF=CFF,
    )

    # ---------------------------------------------------------------
    # Restricted-S raw response
    # ---------------------------------------------------------------

    y_restricted = md.fnecC_composite_raw(
        MS=MS,
        GS=GS,
        S=S_restricted,
        CFB=CFB,
        CFF=CFF,
    )

    # ---------------------------------------------------------------
    # Sobol analysis
    # ---------------------------------------------------------------

    primary = sobol_summary(
        problem,
        y_primary,
        parameter_names,
    )

    primary["S_design"] = (
        "Primary: observed S range"
    )

    restricted = sobol_summary(
        problem,
        y_restricted,
        parameter_names,
    )

    restricted["S_design"] = (
        "Robustness: S P10-P90"
    )

    combined = pd.concat(
        [
            primary,
            restricted,
        ],
        ignore_index=True,
    )

    combined.insert(
        0,
        "formulation",
        "Composite",
    )

    combined.to_csv(
        GSA_DIR
        / "sobol_fnecC5_S_range_robustness.csv",
        index=False,
    )

    # ---------------------------------------------------------------
    # Console output
    # ---------------------------------------------------------------

    print()
    print("Primary raw ST:")
    print(
        primary[
            [
                "parameter",
                "ST",
                "ST_conf",
                "ST_rank",
            ]
        ]
        .sort_values(
            "ST",
            ascending=False,
        )
        .to_string(
            index=False
        )
    )

    print()
    print("S restricted to P10-P90:")
    print(
        restricted[
            [
                "parameter",
                "ST",
                "ST_conf",
                "ST_rank",
            ]
        ]
        .sort_values(
            "ST",
            ascending=False,
        )
        .to_string(
            index=False
        )
    )

    diagnostics = pd.DataFrame({
        "formulation": [
            "Composite",
            "Composite",
        ],
        "S_design": [
            "Primary: observed S range",
            "Robustness: S P10-P90",
        ],
        "S_min": [
            np.min(S_primary),
            np.min(S_restricted),
        ],
        "S_max": [
            np.max(S_primary),
            np.max(S_restricted),
        ],
        "output_median": [
            np.median(y_primary),
            np.median(y_restricted),
        ],
        "output_variance": [
            np.var(
                y_primary,
                ddof=1,
            ),
            np.var(
                y_restricted,
                ddof=1,
            ),
        ],
        "prop_fnecC_gt_100": [
            np.mean(
                y_primary > 100.0
            ),
            np.mean(
                y_restricted > 100.0
            ),
        ],
    })

    return (
        primary,
        restricted,
        combined,
        diagnostics,
    )


# =====================================================================
# 2. TRAIT-RESOLVED MODEL
# =====================================================================

def run_trait_resolved():

    print()
    print("=" * 78)
    print("TRAIT-RESOLVED fnecC")
    print("=" * 78)

    lhs_path = (
        LHS_DIR
        / "fnecC10_lhs.csv"
    )

    if not lhs_path.exists():
        raise FileNotFoundError(
            f"Missing canonical LHS file: {lhs_path}"
        )

    lhs = pd.read_csv(
        lhs_path
    )

    source_names = [
        "MS",
        "GS",
        "S",
        "cB",
        "MGP",
        "MGN",
        "fGP",
        "cF",
        "GF",
        "rB_molar",
    ]

    real_names = [
        "MS",
        "GS",
        "S",
        "cB",
        "MGP",
        "MGN",
        "fGP",
        "cF",
        "GF",
        "rB",
    ]

    missing = [
        x
        for x in source_names
        if x not in lhs.columns
    ]

    if missing:
        raise ValueError(
            f"{lhs_path} missing required columns: {missing}"
        )

    invq = {
        x: make_invq(
            lhs[x].to_numpy()
        )
        for x in source_names
    }

    unit_names = [
        "u_MS",
        "u_GS",
        "u_S",
        "u_cB",
        "u_MGP",
        "u_MGN",
        "u_fGP",
        "u_cF",
        "u_GF",
        "u_rB",
    ]

    problem = {
        "num_vars": 10,
        "names": unit_names,
        "bounds": [[0.0, 1.0]] * 10,
    }

    U = sobol_sample(
        problem,
        N=N_BASE_TRAIT,
        calc_second_order=False,
        scramble=True,
        seed=SEED,
    )

    expected_rows = (
        N_BASE_TRAIT
        * (10 + 2)
    )

    if U.shape != (
        expected_rows,
        10,
    ):
        raise RuntimeError(
            "Unexpected trait-resolved Sobol matrix dimensions."
        )

    # ---------------------------------------------------------------
    # Parameters
    # ---------------------------------------------------------------

    MS = invq["MS"](
        U[:, 0]
    )

    GS = invq["GS"](
        U[:, 1]
    )

    S_primary = invq["S"](
        U[:, 2]
    )

    S_restricted = qS_p10p90(
        U[:, 2]
    )

    cB = invq["cB"](
        U[:, 3]
    )

    MGP = invq["MGP"](
        U[:, 4]
    )

    MGN = invq["MGN"](
        U[:, 5]
    )

    fGP = invq["fGP"](
        U[:, 6]
    )

    cF = invq["cF"](
        U[:, 7]
    )

    GF = invq["GF"](
        U[:, 8]
    )

    rB_molar = invq["rB_molar"](
        U[:, 9]
    )

    # ---------------------------------------------------------------
    # Primary raw response
    # ---------------------------------------------------------------

    y_primary = md.fnecC_trait_raw(
        MS=MS,
        GS=GS,
        S=S_primary,
        cB=cB,
        MGP=MGP,
        MGN=MGN,
        fGP=fGP,
        cF=cF,
        GF=GF,
        rB_molar=rB_molar,
    )

    # ---------------------------------------------------------------
    # Restricted-S raw response
    # ---------------------------------------------------------------

    y_restricted = md.fnecC_trait_raw(
        MS=MS,
        GS=GS,
        S=S_restricted,
        cB=cB,
        MGP=MGP,
        MGN=MGN,
        fGP=fGP,
        cF=cF,
        GF=GF,
        rB_molar=rB_molar,
    )

    # ---------------------------------------------------------------
    # Sobol analysis
    # ---------------------------------------------------------------

    primary = sobol_summary(
        problem,
        y_primary,
        real_names,
    )

    primary["S_design"] = (
        "Primary: observed S range"
    )

    restricted = sobol_summary(
        problem,
        y_restricted,
        real_names,
    )

    restricted["S_design"] = (
        "Robustness: S P10-P90"
    )

    combined = pd.concat(
        [
            primary,
            restricted,
        ],
        ignore_index=True,
    )

    combined.insert(
        0,
        "formulation",
        "Trait-resolved",
    )

    combined.to_csv(
        GSA_DIR
        / "sobol_fnecC10_S_range_robustness.csv",
        index=False,
    )

    # ---------------------------------------------------------------
    # Console output
    # ---------------------------------------------------------------

    print()
    print("Primary raw ST:")
    print(
        primary[
            [
                "parameter",
                "ST",
                "ST_conf",
                "ST_rank",
            ]
        ]
        .sort_values(
            "ST",
            ascending=False,
        )
        .to_string(
            index=False
        )
    )

    print()
    print("S restricted to P10-P90:")
    print(
        restricted[
            [
                "parameter",
                "ST",
                "ST_conf",
                "ST_rank",
            ]
        ]
        .sort_values(
            "ST",
            ascending=False,
        )
        .to_string(
            index=False
        )
    )

    diagnostics = pd.DataFrame({
        "formulation": [
            "Trait-resolved",
            "Trait-resolved",
        ],
        "S_design": [
            "Primary: observed S range",
            "Robustness: S P10-P90",
        ],
        "S_min": [
            np.min(S_primary),
            np.min(S_restricted),
        ],
        "S_max": [
            np.max(S_primary),
            np.max(S_restricted),
        ],
        "output_median": [
            np.median(y_primary),
            np.median(y_restricted),
        ],
        "output_variance": [
            np.var(
                y_primary,
                ddof=1,
            ),
            np.var(
                y_restricted,
                ddof=1,
            ),
        ],
        "prop_fnecC_gt_100": [
            np.mean(
                y_primary > 100.0
            ),
            np.mean(
                y_restricted > 100.0
            ),
        ],
    })

    return (
        primary,
        restricted,
        combined,
        diagnostics,
    )


# =====================================================================
# 3. Run both formulations
# =====================================================================

(
    comp_primary,
    comp_restricted,
    comp_all,
    comp_diag,
) = run_composite()

(
    trait_primary,
    trait_restricted,
    trait_all,
    trait_diag,
) = run_trait_resolved()


# =====================================================================
# 4. Full combined output
# =====================================================================

all_results = pd.concat(
    [
        comp_all,
        trait_all,
    ],
    ignore_index=True,
)

all_results.to_csv(
    GSA_DIR
    / "sobol_fnecC_S_range_robustness_all.csv",
    index=False,
)


# =====================================================================
# 5. Compact manuscript-oriented comparison table
# =====================================================================

def comparison_rows(
    primary,
    restricted,
    formulation,
    parameters,
):

    p = (
        primary[
            primary["parameter"].isin(
                parameters
            )
        ]
        [
            [
                "parameter",
                "S1",
                "S1_conf",
                "S1_rank",
                "ST",
                "ST_conf",
                "ST_rank",
            ]
        ]
        .rename(
            columns={
                "S1": "S1_primary",
                "S1_conf": "S1_conf_primary",
                "S1_rank": "S1_rank_primary",
                "ST": "ST_primary",
                "ST_conf": "ST_conf_primary",
                "ST_rank": "ST_rank_primary",
            }
        )
    )

    r = (
        restricted[
            restricted["parameter"].isin(
                parameters
            )
        ]
        [
            [
                "parameter",
                "S1",
                "S1_conf",
                "S1_rank",
                "ST",
                "ST_conf",
                "ST_rank",
            ]
        ]
        .rename(
            columns={
                "S1": "S1_S_P10_P90",
                "S1_conf": "S1_conf_S_P10_P90",
                "S1_rank": "S1_rank_S_P10_P90",
                "ST": "ST_S_P10_P90",
                "ST_conf": "ST_conf_S_P10_P90",
                "ST_rank": "ST_rank_S_P10_P90",
            }
        )
    )

    out = p.merge(
        r,
        on="parameter",
        how="inner",
    )

    out.insert(
        0,
        "formulation",
        formulation,
    )

    out["ST_change"] = (
        out["ST_S_P10_P90"]
        - out["ST_primary"]
    )

    out["S1_change"] = (
        out["S1_S_P10_P90"]
        - out["S1_primary"]
    )

    return out


compact_composite = comparison_rows(
    comp_primary,
    comp_restricted,
    formulation="Composite",
    parameters=[
        "S",
        "GS",
        "CFF",
    ],
)

compact_trait = comparison_rows(
    trait_primary,
    trait_restricted,
    formulation="Trait-resolved",
    parameters=[
        "S",
        "GS",
        "GF",
    ],
)

compact = pd.concat(
    [
        compact_composite,
        compact_trait,
    ],
    ignore_index=True,
)

compact.to_csv(
    GSA_DIR
    / "sobol_fnecC_S_range_robustness_compact.csv",
    index=False,
)


# =====================================================================
# 6. Diagnostics
# =====================================================================

diagnostics = pd.concat(
    [
        comp_diag,
        trait_diag,
    ],
    ignore_index=True,
)

diagnostics["S_empirical_P10"] = S_P10
diagnostics["S_empirical_P90"] = S_P90

diagnostics.to_csv(
    GSA_DIR
    / "sobol_fnecC_S_range_robustness_diagnostics.csv",
    index=False,
)


# =====================================================================
# 7. Final console summary
# =====================================================================

print()
print("=" * 78)
print("COMPACT ROBUSTNESS COMPARISON")
print("=" * 78)

print(
    compact[
        [
            "formulation",
            "parameter",
            "ST_primary",
            "ST_rank_primary",
            "ST_S_P10_P90",
            "ST_rank_S_P10_P90",
            "ST_change",
        ]
    ]
    .to_string(
        index=False
    )
)

print()
print("=" * 78)
print("DIAGNOSTICS")
print("=" * 78)

print(
    diagnostics[
        [
            "formulation",
            "S_design",
            "S_min",
            "S_max",
            "output_median",
            "output_variance",
            "prop_fnecC_gt_100",
        ]
    ]
    .to_string(
        index=False
    )
)

print()
print("=" * 78)
print("FILES SAVED")
print("=" * 78)

for filename in [
    "sobol_fnecC5_S_range_robustness.csv",
    "sobol_fnecC10_S_range_robustness.csv",
    "sobol_fnecC_S_range_robustness_all.csv",
    "sobol_fnecC_S_range_robustness_compact.csv",
    "sobol_fnecC_S_range_robustness_diagnostics.csv",
]:

    print(
        GSA_DIR
        / filename
    )

print()
print("Robustness analysis complete.")